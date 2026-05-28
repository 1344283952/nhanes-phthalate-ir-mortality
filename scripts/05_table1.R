# ============================================
# 009 / scripts/05_table1.R
# Table 1: 基线特征 by (a) IR binary (HOMA-IR>=2.5) (b) Σ-DEHP 四分位
# 输入: data/processed/nhanes_final.RData + nhanes_design.RData
# 输出: output/tables/table1_baseline_by_ir.csv
#       output/tables/table1_baseline_by_dehp_quartile.csv
#
# 加权统计: weighted median(IQR) for continuous / n (weighted %) for categorical
#           Kruskal-Wallis (svyranktest) / chi² (svychisq) via tableone::svyCreateTableOne
# ============================================

suppressPackageStartupMessages({
  library(survey); library(dplyr); library(tableone)
})

cat("========================================\n")
cat("009 / Table 1: 基线 by IR + by Σ-DEHP quartile\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

# ----------------------------------------------------------
# Prepare: derive DEHP quartile + label IR binary
# ----------------------------------------------------------
dehp_q_cuts <- quantile(nhanes_final$sum_dehp_mol, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
nhanes_final$dehp_q <- cut(nhanes_final$sum_dehp_mol, breaks = dehp_q_cuts,
                          include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
nhanes_final$ir_label <- factor(ifelse(nhanes_final$ir_binary == 1, "IR", "Non-IR"),
                                levels = c("Non-IR","IR"))

# Rebuild design (carries new variables)
design_main <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

# ----------------------------------------------------------
# Variables for Table 1
# ----------------------------------------------------------
cont_vars <- c("age","bmi","waist","pir","homa_ir",
               "LBXGLU","LBXIN","LBXGH",
               "sum_dehp_mol","sum_hmw","sum_lmw")

# 安全地添加 fish_freq_30d / kcal_day / hscrp_log (若存在)
for (v in c("fish_freq_30d","kcal_day","hscrp_log")) {
  if (v %in% names(nhanes_final)) cont_vars <- c(cont_vars, v)
}

cat_vars <- c("age_group","sex_male","race","education","smoke",
              "hypertension","cycle_tag")
# smoke_objective only if non-trivial
if ("smoke_objective" %in% names(nhanes_final) &&
    sum(!is.na(nhanes_final$smoke_objective)) > 100) {
  cat_vars <- c(cat_vars, "smoke_objective")
}

# Cast cat_vars to factor
for (v in cat_vars) {
  if (v %in% names(nhanes_final) && !is.factor(nhanes_final[[v]])) {
    nhanes_final[[v]] <- factor(nhanes_final[[v]])
  }
}
# rebuild design with these
design_main <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

all_vars <- c(cont_vars, cat_vars)

# ----------------------------------------------------------
# (a) Table 1: by IR binary
# ----------------------------------------------------------
cat("--- (a) Table 1: by IR binary ---\n")
tab1a <- tryCatch(
  svyCreateTableOne(vars = all_vars,
                    strata = "ir_label",
                    data = design_main,
                    factorVars = cat_vars,
                    test = TRUE,
                    smd = FALSE),
  error = function(e) { cat("svyCreateTableOne error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(tab1a)) {
  out_a <- print(tab1a,
                 nonnormal = cont_vars,        # report median (IQR) for all continuous
                 showAllLevels = TRUE,
                 catDigits = 1, contDigits = 2, pDigits = 3,
                 printToggle = FALSE,
                 noSpaces = TRUE,
                 test = TRUE)
  cat("\n--- preview ---\n"); print(out_a)
  if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
  write.csv(out_a, "output/tables/table1_baseline_by_ir.csv", row.names = TRUE)
  cat("\n[OK] output/tables/table1_baseline_by_ir.csv\n")
}

# ----------------------------------------------------------
# (b) Table 1: by Σ-DEHP quartile + P-trend
# ----------------------------------------------------------
cat("\n--- (b) Table 1: by Σ-DEHP quartile ---\n")
tab1b <- tryCatch(
  svyCreateTableOne(vars = all_vars,
                    strata = "dehp_q",
                    data = design_main,
                    factorVars = cat_vars,
                    test = TRUE,
                    smd = FALSE),
  error = function(e) { cat("svyCreateTableOne error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(tab1b)) {
  out_b <- print(tab1b,
                 nonnormal = cont_vars,
                 showAllLevels = TRUE,
                 catDigits = 1, contDigits = 2, pDigits = 3,
                 printToggle = FALSE,
                 noSpaces = TRUE,
                 test = TRUE)
  cat("\n--- preview ---\n"); print(out_b)

  # ----------------------------------------------------------
  # P-trend: 把 dehp_q 转 numeric (Q1..Q4 → 1..4) 跑 svyglm/svykm test
  # ----------------------------------------------------------
  nhanes_final$dehp_q_num <- as.integer(nhanes_final$dehp_q)
  design_trend <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                            weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

  ptrend_rows <- list()
  # cont: linear regression on quartile-ordinal
  for (v in cont_vars) {
    fml <- as.formula(sprintf("%s ~ dehp_q_num", v))
    fit <- tryCatch(svyglm(fml, design = design_trend), error = function(e) NULL)
    p_t <- if (!is.null(fit)) summary(fit)$coefficients["dehp_q_num", "Pr(>|t|)"] else NA_real_
    ptrend_rows[[length(ptrend_rows)+1]] <- data.frame(variable = v, p_trend = p_t)
  }
  # cat: chi-sq with linear test (use svyglm with binomial after dichotomizing) - simplified: use Cochran-Armitage
  for (v in cat_vars) {
    # 简单做法: svychisq trend (factor x num) - use svyranktest
    fml <- as.formula(sprintf("dehp_q_num ~ %s", v))
    p_t <- tryCatch({
      kt <- svyranktest(as.formula(sprintf("dehp_q_num ~ %s", v)),
                       design = design_trend, test = "KruskalWallis")
      kt$p.value
    }, error = function(e) NA_real_)
    ptrend_rows[[length(ptrend_rows)+1]] <- data.frame(variable = v, p_trend = p_t)
  }
  ptrend_df <- do.call(rbind, ptrend_rows)
  cat("\n--- P-trend by Σ-DEHP quartile ---\n"); print(ptrend_df)

  # Save Table 1b + p-trend side-by-side
  # 把 ptrend_df 加列到 out_b: 通过 row name 匹配比较粗暴，直接 dump 两份
  if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
  write.csv(out_b, "output/tables/table1_baseline_by_dehp_quartile.csv", row.names = TRUE)
  write.csv(ptrend_df, "output/tables/table1_dehp_quartile_ptrend.csv", row.names = FALSE)
  cat("\n[OK] output/tables/table1_baseline_by_dehp_quartile.csv\n")
  cat("[OK] output/tables/table1_dehp_quartile_ptrend.csv\n")
}

cat("\n========================================\n")
cat("Table 1 done.\n")
cat("========================================\n")
