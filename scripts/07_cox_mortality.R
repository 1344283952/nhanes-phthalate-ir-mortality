# ============================================
# 009 / scripts/07_cox_mortality.R
# Cox 全因 + CM 死亡 (Stack 3 N=2,238)
# 输入: data/processed/nhanes_final.RData + nhanes_design.RData
# 输出: output/tables/cox_mortality_allcause.csv
#       output/tables/cox_mortality_cm.csv
#
# 模型递进调整:
#   Crude: phth ~ 0
#   M1:    + age + sex + race
#   M2:    + edu + pir + bmi + smoke + (drink) + hypertension
#
# 各 Phthalate 单独跑:
#   8 metabolites: URXMEP/MBP/MIB/MZP/MHP/MHH/MOH/ECP (都用 _z continuous)
#   3 个 sum: sum_dehp_mol_z / sum_hmw_z / sum_lmw_z
#
# 输出: HR (95% CI) + p / P-trend (四分位中位数线性进入)
# 用 survey::svycoxph
# ============================================

suppressPackageStartupMessages({
  library(survey); library(dplyr); library(survival); library(broom); library(purrr)
})

cat("========================================\n")
cat("009 / Cox 全因死亡 + CM 死亡 (Stack 3)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

# ----------------------------------------------------------
# Prepare mortality cohort - use df_mort (已 ELIGSTAT==1)
# ----------------------------------------------------------
# W16 R-NHANES C2 fix:
# Use 03_clean's pre-capped permth (max 200 mo biological ceiling).
# Drop rows with NA permth (i.e., outliers > 200 mo set to NA in 03_clean).
n_before_permth <- nrow(df_mort)
df_mort <- df_mort %>% dplyr::filter(!is.na(permth), permth > 0)
cat(sprintf("[R-NHANES C2 fix] PERMTH valid rows: %d / %d (dropped %d with NA/0)\n",
            nrow(df_mort), n_before_permth, n_before_permth - nrow(df_mort)))

df_mort$followup_months <- pmax(df_mort$permth, 0.5)   # 防 0 follow-up; 已 capped at 200 in 03_clean
cat(sprintf("[verify] followup_months range: [%.1f, %.1f] mo (cap = 200 mo)\n",
            min(df_mort$followup_months, na.rm = TRUE),
            max(df_mort$followup_months, na.rm = TRUE)))

cat(sprintf("Stack 3 mortality cohort: N = %d\n", nrow(df_mort)))
cat(sprintf("  All-cause deaths: %d  (%.1f%%)\n",
            sum(df_mort$mort_allcause == 1, na.rm = TRUE),
            100 * mean(df_mort$mort_allcause == 1, na.rm = TRUE)))
cat(sprintf("  CM deaths:        %d  (%.1f%%)\n",
            sum(df_mort$mort_cm == 1, na.rm = TRUE),
            100 * mean(df_mort$mort_cm == 1, na.rm = TRUE)))
cat(sprintf("  Median follow-up (months): %.1f\n",
            median(df_mort$followup_months, na.rm = TRUE)))

# ----------------------------------------------------------
# 协变量准备
# ----------------------------------------------------------
# drink 全 NA → fall back to ALQ101 if avail
if (all(is.na(df_mort$drink)) && "ALQ101" %in% names(df_mort)) {
  df_mort$drink <- factor(case_when(
    df_mort$ALQ101 == 1 ~ "Yes",
    df_mort$ALQ101 == 2 ~ "No"
  ), levels = c("No","Yes"))
  cat(sprintf("[fallback] drink ← ALQ101: non-NA = %d\n", sum(!is.na(df_mort$drink))))
}
use_drink <- !all(is.na(df_mort$drink)) && length(unique(na.omit(df_mort$drink))) > 1

# Cast factors
for (v in c("race","education","smoke","drink","sex_male","hypertension")) {
  if (v %in% names(df_mort) && !is.factor(df_mort[[v]])) {
    df_mort[[v]] <- factor(df_mort[[v]])
  }
}

# Rebuild design with df_mort that contains follow-up
design_mort <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = ~wt_pooled, data = df_mort, nest = TRUE)

# ----------------------------------------------------------
# Exposure list
# ----------------------------------------------------------
phth_exposures <- c(
  "URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
  "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z",
  "sum_dehp_mol_z","sum_hmw_z","sum_lmw_z"
)
phth_exposures <- intersect(phth_exposures, names(df_mort))
cat(sprintf("\n%d exposures to model:\n  %s\n",
            length(phth_exposures), paste(phth_exposures, collapse=", ")))

# ----------------------------------------------------------
# Models
# ----------------------------------------------------------
cov_M1 <- "age + sex_male + race"
cov_M2_base <- "age + sex_male + race + education + pir + bmi + smoke + hypertension"
cov_M2 <- if (use_drink) paste(cov_M2_base, "+ drink") else cov_M2_base
if (!use_drink) cat("[note] drink 排除自 M2 (全 NA / no variation)\n")

# ----------------------------------------------------------
# Helper: 跑单个 phth × model_label × outcome
# ----------------------------------------------------------
run_cox <- function(design, exposure, outcome, model_label, adj) {
  fml <- as.formula(sprintf("Surv(followup_months, %s) ~ %s + %s",
                           outcome, exposure, adj))
  fit <- tryCatch(svycoxph(fml, design = design),
                  error = function(e) {
                    cat(sprintf("[err] %s/%s/%s: %s\n",
                                exposure, outcome, model_label, conditionMessage(e)))
                    NULL
                  })
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  ci <- s$conf.int
  cf <- s$coefficients
  if (!exposure %in% rownames(ci)) return(NULL)
  data.frame(
    exposure = exposure,
    outcome = outcome,
    model = model_label,
    HR = ci[exposure, "exp(coef)"],
    lo = ci[exposure, "lower .95"],
    hi = ci[exposure, "upper .95"],
    p  = cf[exposure, "Pr(>|z|)"],
    n_event = fit$nevent,
    stringsAsFactors = FALSE
  )
}

# Quartile P-trend helper: 用 phth 原值四分位 → median 替换 → linear
run_cox_qtrend <- function(design, exposure_raw, outcome, model_label, adj) {
  # exposure_raw 是 _z 对应的原始值列名 (去 _z)
  # 用 _imp 列计算四分位
  base <- sub("_z$", "_imp", exposure_raw)
  if (!base %in% names(design$variables)) base <- sub("_z$", "", exposure_raw)
  if (!base %in% names(design$variables)) return(NULL)

  x <- design$variables[[base]]
  qs <- quantile(x, probs = c(0,.25,.5,.75,1), na.rm = TRUE)
  if (any(duplicated(qs))) return(NULL)
  q_idx <- cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE)
  med_in_q <- tapply(x, q_idx, median, na.rm = TRUE)
  q_med <- as.numeric(med_in_q)[q_idx]

  design$variables[[paste0(exposure_raw,"_qmed")]] <- q_med
  design$variables[[paste0(exposure_raw,"_q")]] <- factor(q_idx, levels = 1:4,
                                                         labels = c("Q1","Q2","Q3","Q4"))

  # P-trend: continuous q_med
  fml_t <- as.formula(sprintf("Surv(followup_months, %s) ~ %s + %s",
                              outcome, paste0(exposure_raw,"_qmed"), adj))
  fit_t <- tryCatch(svycoxph(fml_t, design = design),
                    error = function(e) NULL)

  qmed_var <- paste0(exposure_raw,"_qmed")
  p_trend <- if (!is.null(fit_t) && qmed_var %in% rownames(summary(fit_t)$coefficients))
    summary(fit_t)$coefficients[qmed_var, "Pr(>|z|)"] else NA_real_

  # Quartile-level HRs (Q2/Q3/Q4 vs Q1)
  fml_q <- as.formula(sprintf("Surv(followup_months, %s) ~ %s + %s",
                              outcome, paste0(exposure_raw,"_q"), adj))
  fit_q <- tryCatch(svycoxph(fml_q, design = design),
                    error = function(e) NULL)
  q_rows <- NULL
  if (!is.null(fit_q)) {
    sq <- summary(fit_q)
    for (qlab in c("Q2","Q3","Q4")) {
      term <- paste0(exposure_raw, "_q", qlab)
      if (term %in% rownames(sq$conf.int)) {
        q_rows <- rbind(q_rows, data.frame(
          exposure = exposure_raw,
          outcome  = outcome,
          model    = model_label,
          term     = paste0("q_",qlab,"_vs_Q1"),
          HR       = sq$conf.int[term,"exp(coef)"],
          lo       = sq$conf.int[term,"lower .95"],
          hi       = sq$conf.int[term,"upper .95"],
          p        = sq$coefficients[term, "Pr(>|z|)"],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  if (!is.null(q_rows)) q_rows$p_trend <- p_trend
  q_rows
}

# ----------------------------------------------------------
# 跑所有
# ----------------------------------------------------------
outcomes <- c("mort_allcause","mort_cm")
specs <- list(
  Crude = "1",
  M1    = cov_M1,
  M2    = cov_M2
)

# 单 phth continuous HR
res_cont <- list()
for (oc in outcomes) {
  for (exp in phth_exposures) {
    for (m in names(specs)) {
      adj <- specs[[m]]
      fml_str <- if (adj == "1") sprintf("Surv(followup_months, %s) ~ %s", oc, exp) else
                                  sprintf("Surv(followup_months, %s) ~ %s + %s", oc, exp, adj)
      fit <- tryCatch(svycoxph(as.formula(fml_str), design = design_mort),
                      error = function(e) {
                        cat(sprintf("[err] %s/%s/%s: %s\n", exp, oc, m, conditionMessage(e)))
                        NULL
                      })
      if (is.null(fit)) next
      s <- summary(fit)
      if (!exp %in% rownames(s$conf.int)) next
      res_cont[[length(res_cont)+1]] <- data.frame(
        exposure = exp,
        outcome  = oc,
        model    = m,
        HR_per_SD = s$conf.int[exp,"exp(coef)"],
        lo        = s$conf.int[exp,"lower .95"],
        hi        = s$conf.int[exp,"upper .95"],
        p_per_SD  = s$coefficients[exp, "Pr(>|z|)"],
        n_event   = fit$nevent,
        stringsAsFactors = FALSE
      )
    }
  }
}
res_cont_df <- do.call(rbind, res_cont)

# Quartile + P-trend
res_q <- list()
for (oc in outcomes) {
  for (exp in phth_exposures) {
    for (m in names(specs)) {
      adj <- specs[[m]]
      qres <- run_cox_qtrend(design_mort, exp, oc, m, adj)
      if (!is.null(qres)) res_q[[length(res_q)+1]] <- qres
    }
  }
}
res_q_df <- do.call(rbind, res_q)

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------
out_all  <- res_cont_df %>% filter(outcome == "mort_allcause")
out_cm   <- res_cont_df %>% filter(outcome == "mort_cm")
out_all_q <- res_q_df   %>% filter(outcome == "mort_allcause")
out_cm_q  <- res_q_df   %>% filter(outcome == "mort_cm")

cat("\n--- All-cause (continuous, per-SD) preview ---\n")
print(out_all, row.names = FALSE)
cat("\n--- All-cause (quartile + P-trend) preview ---\n")
print(out_all_q, row.names = FALSE)
cat("\n--- CM (continuous, per-SD) preview ---\n")
print(out_cm, row.names = FALSE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(out_all,    "output/tables/cox_mortality_allcause.csv",         row.names = FALSE)
write.csv(out_all_q,  "output/tables/cox_mortality_allcause_quartile.csv", row.names = FALSE)
write.csv(out_cm,     "output/tables/cox_mortality_cm.csv",                row.names = FALSE)
write.csv(out_cm_q,   "output/tables/cox_mortality_cm_quartile.csv",       row.names = FALSE)

cat("\n[OK] output/tables/cox_mortality_allcause.csv\n")
cat("[OK] output/tables/cox_mortality_allcause_quartile.csv\n")
cat("[OK] output/tables/cox_mortality_cm.csv\n")
cat("[OK] output/tables/cox_mortality_cm_quartile.csv\n")

# ----------------------------------------------------------
# W16 R-Stats C1 fix: Schoenfeld residual PH test
# Manuscript §2.6 declares "Schoenfeld-residual proportionality tests".
# Implement cox.zph() per Grambsch-Therneau 1994.
# ----------------------------------------------------------
cat("\n--- Schoenfeld PH test (R-Stats C1 fix) ---\n")
schoen_rows <- list()

run_schoenfeld <- function(design, exposure, outcome, model_label, adj) {
  fml_str <- if (adj == "1")
    sprintf("Surv(followup_months, %s) ~ %s", outcome, exposure) else
    sprintf("Surv(followup_months, %s) ~ %s + %s", outcome, exposure, adj)
  fit <- tryCatch(svycoxph(as.formula(fml_str), design = design),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  # cox.zph requires a coxph object — svycoxph extends coxph so cox.zph works
  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (is.null(zph)) return(NULL)
  zt <- as.data.frame(zph$table)
  zt$term <- rownames(zt)
  zt$exposure <- exposure
  zt$outcome  <- outcome
  zt$model    <- model_label
  zt
}

for (oc in outcomes) {
  for (exp in c("sum_dehp_mol_z", "URXMIB_z", "URXMHP_z")) {  # 主 exposures
    if (!exp %in% phth_exposures) next
    for (m in c("M1", "M2")) {
      r <- run_schoenfeld(design_mort, exp, oc, m, specs[[m]])
      if (!is.null(r)) schoen_rows[[length(schoen_rows)+1]] <- r
    }
  }
}

if (length(schoen_rows) > 0) {
  schoen_df <- do.call(rbind, schoen_rows)
  # Order columns: exposure / outcome / model / term / chisq / df / p
  cols_keep <- c("exposure","outcome","model","term","chisq","df","p")
  cols_present <- intersect(cols_keep, names(schoen_df))
  schoen_df <- schoen_df[, c(cols_present, setdiff(names(schoen_df), cols_present))]
  write.csv(schoen_df, "output/tables/cox_schoenfeld.csv", row.names = FALSE)
  cat(sprintf("→ output/tables/cox_schoenfeld.csv (rows=%d)\n", nrow(schoen_df)))
  cat("--- Schoenfeld preview ---\n")
  print(head(schoen_df, 20), row.names = FALSE)
  # Flag any global P < 0.05 violations
  global_rows <- schoen_df[schoen_df$term == "GLOBAL", ]
  if (nrow(global_rows) > 0) {
    cat("\n--- GLOBAL PH test (potential violations P < 0.05) ---\n")
    print(global_rows[, c("exposure","outcome","model","chisq","p")], row.names = FALSE)
  }
} else {
  cat("[WARN] No Schoenfeld outputs collected\n")
}

cat("\n========================================\n")
cat("Cox done.\n")
cat("========================================\n")
