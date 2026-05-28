# ============================================
# 009 / scripts/18_evalue_ci_bound.R
# E-value (CI lower bound, VanderWeele 2017 Supplementary)
#
# 同 14_evalue 但额外汇出"CI lower bound E-value":
#   - reviewer 标准武器 (point estimate E-value 单看不够,
#     因为 point 估高 + CI 跨 1 时 E-value 反而虚高)
#   - CI-bound E-value = OR/HR 的 95% CI 下界对应的 E-value
#   - 若 CI 跨 null (lower < 1) -> E-value(CI lower) = 1.00 (无 robustness)
#   - 若 CI lower > 1 -> 报 E-value of that bound
#
# 输出: output/tables/evalue_ci_bound.csv (Point + CIbound 双套)
# ============================================

suppressPackageStartupMessages({
  library(survey); library(survival); library(dplyr); library(EValue); library(broom); library(purrr)
})

cat("========================================\n")
cat("009 / 18_evalue_ci_bound.R: E-value (Point + CI lower bound)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

exposures <- c(
  MEP        = "URXMEP_z",
  MnBP       = "URXMBP_z",
  MiBP       = "URXMIB_z",
  MBzP       = "URXMZP_z",
  MEHP       = "URXMHP_z",
  MEHHP      = "URXMHH_z",
  MEOHP      = "URXMOH_z",
  MECPP      = "URXECP_z",
  `Sum-DEHP` = "sum_dehp_mol_z"
)
exposures <- exposures[exposures %in% names(nhanes_final)]

# W16 R-Causal CI-1 fix: cov_pre Pearl-backdoor set (drop bmi/waist/hypertension)
cov_pre <- c("age","sex_male","race","education","pir","smoke","cotinine_log")
cov_base <- cov_pre  # alias

# ------------------------------------------------------------------
# Helper: pull "point" and "lower" E-values from evalues output
# ------------------------------------------------------------------
ev_safe <- function(ev_obj, which_col) {
  v <- tryCatch(ev_obj["E-values", which_col], error=function(e) NA_real_)
  if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)
}

# ------------------------------------------------------------------
# Outcome 1: IR binary (logistic, OR, rare=FALSE since ~47% prevalence)
# ------------------------------------------------------------------
cat("\n--- Outcome 1: IR binary (logistic OR) ---\n")
run_ir_logistic <- function(exp_lab, exp_var) {
  f <- as.formula(paste("ir_binary ~", exp_var, "+", paste(cov_base, collapse=" + ")))
  m <- tryCatch(svyglm(f, design=design_main, family=quasibinomial()),
                error=function(e) NULL)
  if (is.null(m)) return(NULL)
  tt <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term == exp_var)
  if (nrow(tt) == 0) return(NULL)
  ev <- EValue::evalues.OR(est = tt$estimate, lo = tt$conf.low, hi = tt$conf.high,
                           rare = FALSE)
  data.frame(
    exposure = exp_lab, outcome = "IR_binary", measure = "OR",
    est_point = tt$estimate, CI_low = tt$conf.low, CI_high = tt$conf.high,
    p_value = tt$p.value,
    RR_approx_point = ev_safe(ev, "point"),
    RR_approx_CIlower = ev_safe(ev, "lower"),
    Evalue_point = ev_safe(ev, "point"),
    Evalue_CIbound = ev_safe(ev, "lower"),
    crosses_null = ifelse(tt$conf.low <= 1 & tt$conf.high >= 1, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}
ir_rows <- map_dfr(names(exposures), ~ run_ir_logistic(.x, exposures[[.x]]))

# Fix: evalues.OR returns RR-approx in row "RR" cols, E-value in row "E-values" cols
# We want RR_point = ev["RR","point"], E_point = ev["E-values","point"], etc.
fix_ev_cols <- function(rows_df, fitter) {
  out <- map_dfr(seq_len(nrow(rows_df)), function(i) {
    r <- rows_df[i, ]
    ev <- if (r$measure == "OR") {
      EValue::evalues.OR(est = r$est_point, lo = r$CI_low, hi = r$CI_high, rare = FALSE)
    } else if (r$measure == "HR") {
      EValue::evalues.HR(est = r$est_point, lo = r$CI_low, hi = r$CI_high, rare = TRUE)
    } else NA
    r$RR_approx_point   <- ev["RR","point"]
    r$RR_approx_CIlower <- ev["RR","lower"]
    r$Evalue_point      <- ev["E-values","point"]
    r$Evalue_CIbound    <- ev["E-values","lower"]
    r
  })
  out
}
ir_rows <- fix_ev_cols(ir_rows)
print(ir_rows, row.names=FALSE)

# ------------------------------------------------------------------
# Outcome 2: All-cause mortality (Cox, HR, rare=TRUE since ~8.6%)
# ------------------------------------------------------------------
cat("\n--- Outcome 2: All-cause mortality (svycoxph HR) ---\n")
df_mort_use <- df_mort %>% dplyr::filter(!is.na(permth), !is.na(mort_allcause), permth > 0)
cat(sprintf("Mortality analytic N = %d (events = %d)\n",
            nrow(df_mort_use), sum(df_mort_use$mort_allcause==1, na.rm=TRUE)))
design_mort2 <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_pooled,
                          data=df_mort_use, nest=TRUE)

run_mort_cox <- function(exp_lab, exp_var) {
  f <- as.formula(paste("Surv(permth, mort_allcause) ~", exp_var, "+",
                        paste(cov_base, collapse=" + ")))
  m <- tryCatch(svycoxph(f, design=design_mort2), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)
  ci <- which(rownames(s$coefficients) == exp_var)
  if (length(ci) == 0) return(NULL)
  hr <- exp(s$coefficients[ci,"coef"]); se <- s$coefficients[ci,"se(coef)"]
  ci_lo <- exp(s$coefficients[ci,"coef"] - 1.96*se)
  ci_hi <- exp(s$coefficients[ci,"coef"] + 1.96*se)
  p_val <- s$coefficients[ci,"Pr(>|z|)"]
  data.frame(
    exposure = exp_lab, outcome = "All_cause_mortality", measure = "HR",
    est_point = hr, CI_low = ci_lo, CI_high = ci_hi, p_value = p_val,
    RR_approx_point = NA, RR_approx_CIlower = NA,
    Evalue_point = NA, Evalue_CIbound = NA,
    crosses_null = ifelse(ci_lo <= 1 & ci_hi >= 1, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}
mort_rows <- map_dfr(names(exposures), ~ run_mort_cox(.x, exposures[[.x]]))
mort_rows <- fix_ev_cols(mort_rows)
print(mort_rows, row.names=FALSE)

# ------------------------------------------------------------------
# Combine + save
# ------------------------------------------------------------------
all_rows <- bind_rows(ir_rows, mort_rows)

# Round
all_rows$est_point         <- round(all_rows$est_point, 4)
all_rows$CI_low            <- round(all_rows$CI_low, 4)
all_rows$CI_high           <- round(all_rows$CI_high, 4)
all_rows$p_value           <- signif(all_rows$p_value, 4)
all_rows$RR_approx_point   <- round(all_rows$RR_approx_point, 4)
all_rows$RR_approx_CIlower <- round(all_rows$RR_approx_CIlower, 4)
all_rows$Evalue_point      <- round(all_rows$Evalue_point, 4)
all_rows$Evalue_CIbound    <- round(all_rows$Evalue_CIbound, 4)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(all_rows, "output/tables/evalue_ci_bound.csv", row.names = FALSE)

cat("\n========================================\n")
cat(sprintf("Saved %d rows: output/tables/evalue_ci_bound.csv\n", nrow(all_rows)))
cat("\nInterpretation (VanderWeele 2017):\n")
cat("  Evalue_point   = E-value for point estimate (OR/HR).\n")
cat("  Evalue_CIbound = E-value for CI lower bound (more conservative,\n")
cat("                   reviewer-standard武器 -- 若 CI 跨 null = 1.00 / NA -> 无 robustness).\n")
cat("  Both > 2 = robust to moderate unmeasured confounding.\n")
cat("========================================\n")
