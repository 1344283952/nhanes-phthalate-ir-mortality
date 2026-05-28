# ============================================
# 009 / scripts/14_evalue.R
# E-value (VanderWeele & Ding 2017) for residual confounding sensitivity
#
# 对各 Phthalate (8 metabolites + Σ-DEHP) × 2 outcome (IR binary logistic +
# all-cause mortality Cox) 计算 point estimate 的 E-value。
#
# 解读 (VanderWeele 2017 Ann Intern Med):
#   - E-value > 2 = 需要一个未测量混杂同时与暴露 & 结局关联强度 ≥ 2 才能 explain away
#   - E-value > 1.5 = "moderate-to-strong" robustness
#   - E-value ≈ 1 = 几乎任何未测量混杂都能 explain away
#
# 输入: data/processed/nhanes_design.RData
# 输出: output/tables/evalue_phth_ir.csv
#        output/tables/evalue_phth_mort.csv
# ============================================

suppressPackageStartupMessages({
  library(survey); library(survival); library(dplyr); library(EValue); library(broom)
})

cat("========================================\n")
cat("009 / 14_evalue.R: E-value (point estimate)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

# 暴露 panel: 8 metabolites + Σ-DEHP
exposures <- c(
  MEP    = "URXMEP_z",
  MnBP   = "URXMBP_z",
  MiBP   = "URXMIB_z",
  MBzP   = "URXMZP_z",
  MEHP   = "URXMHP_z",
  MEHHP  = "URXMHH_z",
  MEOHP  = "URXMOH_z",
  MECPP  = "URXECP_z",
  `Sum-DEHP` = "sum_dehp_mol_z"
)
exposures <- exposures[exposures %in% names(nhanes_final)]
cat(sprintf("可用暴露 (%d):\n", length(exposures)))
print(exposures)

# W16 R-Causal CI-1 fix:
# Pair E-value with the headline causal estimand (CMAverse Rte) — must use the
# Pearl-backdoor cov_pre set, NOT the CDE-style M2 cov_base that includes
# bmi/waist/hypertension (which are mediator-confounders downstream of X).
# cov_pre = pre-exposure baseline confounders ONLY.
# 协变量 (Pearl-backdoor cov_pre)
cov_pre <- c("age","sex_male","race","education","pir","smoke","cotinine_log")
# Keep cov_base as documentation backup of CDE-style M2 (NOT used in E-value)
cov_base <- cov_pre  # alias for backward compatibility with rest of script
# (drink 全 NA, 跳过; bmi/waist/hypertension dropped per Pearl-backdoor)

# ------------------------------------------------------------------
# Outcome 1: IR binary (logistic, OR)
# ------------------------------------------------------------------
cat("\n--- Outcome 1: IR binary (HOMA-IR >= 2.5) ---\n")

run_ir_logistic <- function(exp_lab, exp_var) {
  f <- as.formula(paste("ir_binary ~", exp_var, "+",
                        paste(cov_base, collapse = " + ")))
  m <- tryCatch(
    svyglm(f, design = design_main, family = quasibinomial()),
    error = function(e) { cat("  [err]", exp_lab, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  tt <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term == exp_var)
  if (nrow(tt) == 0) return(NULL)

  # VanderWeele 2017: for logistic with common outcome (>15%) use sqrt(OR) approx
  # IR prevalence ~47% -> rare = FALSE
  ev <- EValue::evalues.OR(est = tt$estimate,
                           lo  = tt$conf.low,
                           hi  = tt$conf.high,
                           rare = FALSE)
  data.frame(
    exposure  = exp_lab,
    outcome   = "IR_binary",
    OR        = tt$estimate,
    CI_low    = tt$conf.low,
    CI_high   = tt$conf.high,
    p_value   = tt$p.value,
    # ev row 1 = RR approx, row 2 = E-values
    RR_approx = ev["RR","point"],
    Evalue_point = ev["E-values","point"],
    Evalue_CIbound = ev["E-values","lower"],
    stringsAsFactors = FALSE
  )
}

ir_rows <- purrr::map_dfr(names(exposures), ~ run_ir_logistic(.x, exposures[[.x]]))
print(ir_rows, row.names = FALSE)

# ------------------------------------------------------------------
# Outcome 2: All-cause mortality (Cox HR via svycoxph)
# ------------------------------------------------------------------
cat("\n--- Outcome 2: All-cause mortality (svycoxph) ---\n")

# Build new design for mortality cohort (must have permth + mort_allcause non-NA)
df_mort_use <- df_mort %>%
  dplyr::filter(!is.na(permth), !is.na(mort_allcause), permth > 0)
cat(sprintf("Mortality analytic N = %d (events = %d)\n",
            nrow(df_mort_use), sum(df_mort_use$mort_allcause == 1, na.rm=TRUE)))

design_mort2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                          weights = ~wt_pooled, data = df_mort_use, nest = TRUE)

# Rare outcome (events ~8.6%) -> rare = TRUE
run_mort_cox <- function(exp_lab, exp_var) {
  f <- as.formula(paste("Surv(permth, mort_allcause) ~", exp_var, "+",
                        paste(cov_base, collapse = " + ")))
  m <- tryCatch(
    svycoxph(f, design = design_mort2),
    error = function(e) { cat("  [err]", exp_lab, ":", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  s <- summary(m)
  coef_row <- which(rownames(s$coefficients) == exp_var)
  if (length(coef_row) == 0) return(NULL)
  hr   <- exp(s$coefficients[coef_row, "coef"])
  se   <- s$coefficients[coef_row, "se(coef)"]
  ci_lo <- exp(s$coefficients[coef_row, "coef"] - 1.96 * se)
  ci_hi <- exp(s$coefficients[coef_row, "coef"] + 1.96 * se)
  p    <- s$coefficients[coef_row, "Pr(>|z|)"]
  ev   <- EValue::evalues.HR(est = hr, lo = ci_lo, hi = ci_hi, rare = TRUE)
  data.frame(
    exposure  = exp_lab,
    outcome   = "All_cause_mortality",
    HR        = hr,
    CI_low    = ci_lo,
    CI_high   = ci_hi,
    p_value   = p,
    RR_approx = ev["RR","point"],
    Evalue_point = ev["E-values","point"],
    Evalue_CIbound = ev["E-values","lower"],
    stringsAsFactors = FALSE
  )
}

mort_rows <- purrr::map_dfr(names(exposures), ~ run_mort_cox(.x, exposures[[.x]]))
print(mort_rows, row.names = FALSE)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(ir_rows,   "output/tables/evalue_phth_ir.csv",   row.names = FALSE)
write.csv(mort_rows, "output/tables/evalue_phth_mort.csv", row.names = FALSE)

cat("\n========================================\n")
cat(sprintf("已保存 %d IR + %d mortality E-value rows\n",
            nrow(ir_rows), nrow(mort_rows)))
cat("  output/tables/evalue_phth_ir.csv\n")
cat("  output/tables/evalue_phth_mort.csv\n")
cat("\nE-value > 2 = robust to moderate unmeasured confounding (VanderWeele 2017).\n")
cat("========================================\n")
