# ============================================
# 009 / scripts/12_cmaverse_4way.R
# CMAverse 4-way decomposition (VanderWeele 2014 Epidemiology
# DOI 10.1097/EDE.0000000000000121 + Shi 2021 Stat Methods Med Res
# DOI 10.1177/09622802211009243)
#
# 输入: data/processed/nhanes_final.RData
# 输出:
#   output/tables/cmaverse_phth_ir_4way.csv     (Y = IR binary; X = sum_dehp_mol_z;
#                                                 M = adiposity composite z-score)
#   output/tables/cmaverse_phth_mort_4way.csv   (Y = all-cause mortality;
#                                                 X = sum_dehp_mol_z; M = HOMA-IR log)
#   output/tables/cmaverse_phth_4way.RData      (full cm_fit objects)
#
# 4-way decomposition components (VanderWeele 2014):
#   CDE      — Controlled Direct Effect (M held at reference)
#   INT_ref  — Reference Interaction (exposure-only effect at M=ref)
#   INT_med  — Mediated Interaction (interaction × M change)
#   PIE      — Pure Indirect Effect
#   Total    = CDE + INT_ref + INT_med + PIE
#
# Mediator choice per task:
#   IR outcome  (logistic): M1=HSCRP / M2=adiposity / M3=HOMA-IR
#     ⇒ HSCRP only 368/2239 complete (16%) → using BMI/waist composite as M
#       (adiposity proxy retains N=2,200+ analytic sample)
#     ⇒ M3=HOMA-IR is the outcome proxy (collinear with IR binary) → skipped
#       per task instruction "此处 IR 是 outcome, skip M3"
#   Mortality outcome (Cox): M = HOMA-IR (log) is canonical IR mediator
#
# Method = "rb" (regression-based) — postc empty (no exposure-induced
#   mediator-outcome confounders explicitly modeled; pre-exposure C only).
#   Rationale: with M = adiposity proxy or HOMA-IR proxy, the BMI/waist/HOMA-IR
#   ARE the mediators; the L1 mediator-confounder structure that forced gformula
#   in 005 (BMI affecting M=FIB-4) does not apply here. R-Causal-Methods
#   limitation: rb is consistent under Pearl 2009 Theorem 4.5.1 if no exposure-
#   induced M-Y confounders. Document in Limitations section.
#
# Bootstrap: nboot = 500
# CMAverse has no native svy support → main 4-way unweighted;
#   weighted total-effect svyglm comparison saved as sensitivity row.
# Seed: 20260524
# ============================================

set.seed(20260524)

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survey)
  library(CMAverse)
})

cat("========================================\n")
cat("009 / 12_cmaverse_4way — Phthalate × {IR, mortality} 4-way decomposition\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final n=%d ; df_mort n=%d (all-cause=%d / CM=%d)\n",
            nrow(nhanes_final), nrow(df_mort),
            sum(df_mort$mort_allcause, na.rm = TRUE),
            sum(df_mort$mort_cm, na.rm = TRUE)))

# ---------------------------------------------------------------
# Step 0: Common covariates + exposure construct
# ---------------------------------------------------------------
# X (main exposure): sum_dehp_mol_z (z-score of mole-weighted Σ-DEHP)
# Pre-exposure baseline C (basec):
#   age, sex_male, race, education, pir, smoke_ever, htn_med, cotinine_log
# Note: 'drink' completely missing in 03_clean → not used
nhanes_final$smoke_ever <- ifelse(!is.na(nhanes_final$smoke) &
                                    nhanes_final$smoke == "Ever", 1L, 0L)
df_mort$smoke_ever      <- ifelse(!is.na(df_mort$smoke) &
                                    df_mort$smoke == "Ever", 1L, 0L)

# Median impute cotinine_log
med_cot_n <- median(nhanes_final$cotinine_log, na.rm = TRUE)
nhanes_final$cotinine_log[is.na(nhanes_final$cotinine_log)] <- med_cot_n
df_mort$cotinine_log[is.na(df_mort$cotinine_log)]           <- med_cot_n

# race factor → numeric levels not needed; CMAverse handles factors
# Coerce race to factor for regression compatibility
for (df_name in c("nhanes_final", "df_mort")) {
  d <- get(df_name)
  if (!is.factor(d$race)) d$race <- factor(d$race)
  if (!is.factor(d$education)) d$education <- factor(d$education)
  assign(df_name, d)
}

# W16 R-Causal C2 / R-Stats C7 fix:
# basec = pre-exposure baseline confounders ONLY (Pearl-backdoor cov_pre).
# Drop hypertension — it is L1 exposure-induced mediator-confounder per
# Methods §2.4 + Supp §S8 DAG declaration (downstream of X). Including it
# in basec would partially block X → adiposity → IR pathway through HTN.
basec_vars <- c("age", "sex_male", "race", "education", "pir",
                "smoke_ever", "cotinine_log")

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ===============================================================
# PART 1: IR outcome (logistic) — X → adiposity → IR binary
# ===============================================================
cat("\n========== PART 1: X=Σ-DEHP-z, M=adiposity, Y=IR binary ==========\n")

# Build adiposity composite (BMI + waist standardized average)
nhanes_final$bmi_z   <- as.numeric(scale(nhanes_final$bmi))
if ("waist" %in% names(nhanes_final)) {
  nhanes_final$waist_z <- as.numeric(scale(nhanes_final$waist))
  nhanes_final$adiposity_z <- rowMeans(cbind(nhanes_final$bmi_z,
                                              nhanes_final$waist_z),
                                        na.rm = TRUE)
} else {
  nhanes_final$adiposity_z <- nhanes_final$bmi_z
}
nhanes_final$adiposity_z[is.nan(nhanes_final$adiposity_z)] <- NA

ir_vars <- c("sum_dehp_mol_z", "adiposity_z", "ir_binary",
             basec_vars, "wt_pooled", "SDMVPSU", "SDMVSTRA")

df_ir <- nhanes_final %>%
  select(any_of(ir_vars)) %>%
  filter(if_all(c("sum_dehp_mol_z", "adiposity_z", "ir_binary",
                  "age", "sex_male", "race", "education", "pir",
                  "smoke_ever"), ~ !is.na(.)))
cat(sprintf("CMAverse IR analytic n=%d (cases=%d)\n",
            nrow(df_ir), sum(df_ir$ir_binary == 1)))

# W16 R-Stats C8 fix: Q4-vs-Q1 contrast (not IQR P25-vs-P75)
# Manuscript §2.5 (v) + Figure 5 + §3.5 all claim "Q4 vs Q1 contrast".
# Compute Q4 median (top quartile centroid) vs Q1 median (bottom quartile centroid).
qs_dehp <- quantile(df_ir$sum_dehp_mol_z, probs = c(0, .25, .5, .75, 1.0), na.rm = TRUE)
q_idx_ir <- cut(df_ir$sum_dehp_mol_z, breaks = qs_dehp, include.lowest = TRUE, labels = FALSE)
q1_med <- median(df_ir$sum_dehp_mol_z[q_idx_ir == 1], na.rm = TRUE)
q4_med <- median(df_ir$sum_dehp_mol_z[q_idx_ir == 4], na.rm = TRUE)
cat(sprintf("[R-Stats C8 fix] Q1 median = %.4f, Q4 median = %.4f (Q4-vs-Q1 contrast spans %.4f)\n",
            q1_med, q4_med, q4_med - q1_med))

cat("\n[1/2] cmest (rb method, 4-way) → IR binary [B=500] ...\n")
t0 <- Sys.time()
cm_ir <- tryCatch({
  CMAverse::cmest(
    data       = df_ir,
    model      = "rb",
    outcome    = "ir_binary",
    exposure   = "sum_dehp_mol_z",
    mediator   = "adiposity_z",
    basec      = basec_vars,
    EMint      = TRUE,
    mreg       = list("linear"),
    yreg       = "logistic",
    # W16 R-Stats C8 fix: Q4-vs-Q1 contrast (Q1 median → Q4 median)
    astar      = q1_med,
    a          = q4_med,
    mval       = list(median(df_ir$adiposity_z, na.rm = TRUE)),
    estimation = "imputation",
    inference  = "bootstrap",
    nboot      = 500,
    boot.ci.type = "per",
    multimp    = FALSE
  )
}, error = function(e) {
  cat("ERROR cmest IR:", conditionMessage(e), "\n")
  NULL
})
cat(sprintf("  fit time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Extract 4-way components
extract_cmest_table <- function(cm_fit, outcome_lbl) {
  if (is.null(cm_fit)) return(data.frame())
  ss <- tryCatch(summary(cm_fit), error = function(e) NULL)
  if (is.null(ss)) return(data.frame())
  if (is.null(ss$effect.pe)) return(data.frame())
  ci_lo <- if (!is.null(ss$effect.ci.low))  as.numeric(ss$effect.ci.low)  else rep(NA, length(ss$effect.pe))
  ci_hi <- if (!is.null(ss$effect.ci.high)) as.numeric(ss$effect.ci.high) else rep(NA, length(ss$effect.pe))
  pvs   <- if (!is.null(ss$effect.pval))    as.numeric(ss$effect.pval)    else rep(NA, length(ss$effect.pe))
  data.frame(
    outcome    = outcome_lbl,
    component  = names(ss$effect.pe),
    estimate   = as.numeric(ss$effect.pe),
    ci_lcl     = ci_lo,
    ci_ucl     = ci_hi,
    p          = pvs,
    stringsAsFactors = FALSE
  )
}

ir_4way <- extract_cmest_table(cm_ir, "IR binary (HOMA>=2.5)")
if (nrow(ir_4way) > 0) {
  # Format effect_str
  ir_4way$effect_str <- sprintf("%.3f (%.3f to %.3f), p=%.4g",
                                 ir_4way$estimate, ir_4way$ci_lcl,
                                 ir_4way$ci_ucl, ir_4way$p)
  write.csv(ir_4way, "output/tables/cmaverse_phth_ir_4way.csv", row.names = FALSE)
  cat("\n→ output/tables/cmaverse_phth_ir_4way.csv\n")
  print(ir_4way[, c("component","estimate","ci_lcl","ci_ucl","p")])
} else {
  write.csv(data.frame(note = "cmest IR fit failed or summary unavailable"),
            "output/tables/cmaverse_phth_ir_4way.csv", row.names = FALSE)
}

# ===============================================================
# PART 2: Mortality outcome (Cox) — X → HOMA-IR (log) → mortality
# ===============================================================
cat("\n========== PART 2: X=Σ-DEHP-z, M=HOMA-IR(log), Y=all-cause mort (Cox) ==========\n")

df_mort$homa_ir_log <- if ("homa_ir_log" %in% names(df_mort)) {
  df_mort$homa_ir_log
} else {
  log(pmax(df_mort$homa_ir, 0.01))
}
df_mort$sum_dehp_mol_z <- if ("sum_dehp_mol_z" %in% names(df_mort)) {
  df_mort$sum_dehp_mol_z
} else {
  nhanes_final$sum_dehp_mol_z[match(df_mort$SEQN, nhanes_final$SEQN)]
}

mort_vars <- c("sum_dehp_mol_z", "homa_ir_log", "permth", "mort_allcause",
               basec_vars, "wt_pooled", "SDMVPSU", "SDMVSTRA")

df_m <- df_mort %>%
  select(any_of(mort_vars)) %>%
  filter(if_all(c("sum_dehp_mol_z", "homa_ir_log", "permth", "mort_allcause",
                  "age", "sex_male", "race", "education", "pir",
                  "smoke_ever"), ~ !is.na(.)))
cat(sprintf("CMAverse mort analytic n=%d (events=%d)\n",
            nrow(df_m), sum(df_m$mort_allcause)))

# W16 R-Stats C8 fix (mortality): Q4-vs-Q1 contrast
qs_dehp_m <- quantile(df_m$sum_dehp_mol_z, probs = c(0, .25, .5, .75, 1.0), na.rm = TRUE)
q_idx_m <- cut(df_m$sum_dehp_mol_z, breaks = qs_dehp_m, include.lowest = TRUE, labels = FALSE)
q1_med_m <- median(df_m$sum_dehp_mol_z[q_idx_m == 1], na.rm = TRUE)
q4_med_m <- median(df_m$sum_dehp_mol_z[q_idx_m == 4], na.rm = TRUE)
cat(sprintf("[R-Stats C8 fix] Q1 median = %.4f, Q4 median = %.4f (Q4-vs-Q1 contrast spans %.4f)\n",
            q1_med_m, q4_med_m, q4_med_m - q1_med_m))

cat("\n[2/2] cmest (rb method, 4-way) → all-cause mortality (Cox) [B=500] ...\n")
t0 <- Sys.time()
cm_mort <- tryCatch({
  CMAverse::cmest(
    data       = df_m,
    model      = "rb",
    outcome    = "permth",
    event      = "mort_allcause",
    exposure   = "sum_dehp_mol_z",
    mediator   = "homa_ir_log",
    basec      = basec_vars,
    EMint      = TRUE,
    mreg       = list("linear"),
    yreg       = "coxph",
    # W16 R-Stats C8 fix: Q4-vs-Q1 contrast (Q1 median → Q4 median)
    astar      = q1_med_m,
    a          = q4_med_m,
    mval       = list(median(df_m$homa_ir_log, na.rm = TRUE)),
    estimation = "imputation",
    inference  = "bootstrap",
    nboot      = 500,
    boot.ci.type = "per",
    multimp    = FALSE
  )
}, error = function(e) {
  cat("ERROR cmest mort:", conditionMessage(e), "\n")
  NULL
})
cat(sprintf("  fit time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

mort_4way <- extract_cmest_table(cm_mort, "All-cause mortality")
if (nrow(mort_4way) > 0) {
  mort_4way$effect_str <- sprintf("%.3f (%.3f to %.3f), p=%.4g",
                                   mort_4way$estimate, mort_4way$ci_lcl,
                                   mort_4way$ci_ucl, mort_4way$p)
  write.csv(mort_4way, "output/tables/cmaverse_phth_mort_4way.csv", row.names = FALSE)
  cat("\n→ output/tables/cmaverse_phth_mort_4way.csv\n")
  print(mort_4way[, c("component","estimate","ci_lcl","ci_ucl","p")])
} else {
  write.csv(data.frame(note = "cmest mort fit failed or summary unavailable"),
            "output/tables/cmaverse_phth_mort_4way.csv", row.names = FALSE)
}

# ===============================================================
# Save full + weighted comparison
# ===============================================================
# Weighted Cox sensitivity (total effect, unweighted vs svy)
# W16 R-Stats C7 fix: regressors match basec (Pearl-backdoor; drop hypertension)
cat("\n[3/3] Weighted (svycoxph) vs unweighted Cox total effect — sensitivity check ...\n")
comparison <- tryCatch({
  des_w <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                     weights = ~wt_pooled, data = df_m, nest = TRUE)
  un <- coxph(Surv(permth, mort_allcause) ~ sum_dehp_mol_z + age + sex_male + race +
                education + pir + smoke_ever, data = df_m)
  wt <- survey::svycoxph(Surv(permth, mort_allcause) ~ sum_dehp_mol_z + age + sex_male + race +
                          education + pir + smoke_ever, design = des_w)
  data.frame(
    model = c("Unweighted Cox", "Weighted svycoxph"),
    hr = c(exp(coef(un)["sum_dehp_mol_z"]), exp(coef(wt)["sum_dehp_mol_z"])),
    hr_lcl = c(exp(confint(un)["sum_dehp_mol_z", 1]),
               exp(confint(wt)["sum_dehp_mol_z", 1])),
    hr_ucl = c(exp(confint(un)["sum_dehp_mol_z", 2]),
               exp(confint(wt)["sum_dehp_mol_z", 2])),
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  cat("WARN sensitivity Cox failed:", conditionMessage(e), "\n")
  data.frame()
})
if (nrow(comparison) > 0) {
  write.csv(comparison, "output/tables/cmaverse_phth_wt_unwt_compare.csv",
            row.names = FALSE)
  cat("\nWeighted vs unweighted total effect:\n"); print(comparison)
}

# W16 R-Stats M6 fix: Weighted CMAverse IR 4-way sensitivity
# Use IPTW-style pseudo-weighting via subsample, given CMAverse has no native svy.
# Approach: weighted logistic + linear regressions in svyglm to recover Rte/Rpnie
# components manually as comparison row (no full 4-way decomposition, but Rte +
# direction check). Saves to cmaverse_phth_ir_weighted_sensitivity.csv.
cat("\n[Extra] Weighted IR Rte vs unweighted Rte comparison (R-Stats M6) ...\n")
ir_wt_comp <- tryCatch({
  des_ir <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                      weights = ~wt_pooled, data = df_ir, nest = TRUE)
  # Total effect (Rte) — exposure on outcome, adjusted for basec only
  fml_un <- as.formula(paste("ir_binary ~ sum_dehp_mol_z +",
                             paste(basec_vars, collapse = " + ")))
  un_ir <- glm(fml_un, data = df_ir, family = quasibinomial())
  wt_ir <- survey::svyglm(fml_un, design = des_ir, family = quasibinomial())
  # Q4 vs Q1 contrast on linear predictor
  beta_un  <- coef(un_ir)["sum_dehp_mol_z"]
  beta_wt  <- coef(wt_ir)["sum_dehp_mol_z"]
  se_un    <- summary(un_ir)$coefficients["sum_dehp_mol_z", "Std. Error"]
  se_wt    <- summary(wt_ir)$coefficients["sum_dehp_mol_z", "Std. Error"]
  contrast <- q4_med - q1_med
  data.frame(
    model = c("Unweighted total effect (IR)", "Weighted svyglm total effect (IR)"),
    contrast = c(sprintf("Q4-Q1 (%.3f units)", contrast)),
    or_per_contrast = c(exp(beta_un * contrast), exp(beta_wt * contrast)),
    lcl = c(exp((beta_un - 1.96*se_un) * contrast),
            exp((beta_wt - 1.96*se_wt) * contrast)),
    ucl = c(exp((beta_un + 1.96*se_un) * contrast),
            exp((beta_wt + 1.96*se_wt) * contrast)),
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  cat("WARN weighted IR Rte sensitivity failed:", conditionMessage(e), "\n")
  data.frame()
})
if (nrow(ir_wt_comp) > 0) {
  write.csv(ir_wt_comp, "output/tables/cmaverse_phth_ir_weighted_sensitivity.csv",
            row.names = FALSE)
  cat("\nIR weighted-vs-unweighted Rte comparison:\n"); print(ir_wt_comp)
}

save(cm_ir, cm_mort, ir_4way, mort_4way, comparison,
     file = "output/tables/cmaverse_phth_4way.RData")

cat("\n保存:\n")
cat("  output/tables/cmaverse_phth_ir_4way.csv\n")
cat("  output/tables/cmaverse_phth_mort_4way.csv\n")
cat("  output/tables/cmaverse_phth_wt_unwt_compare.csv\n")
cat("  output/tables/cmaverse_phth_4way.RData\n")
cat("\nDONE 12_cmaverse_4way.R\n")
