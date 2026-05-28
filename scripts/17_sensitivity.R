# ============================================
# 009 / scripts/17_sensitivity.R
# 8 套敏感性分析 (Plan M13-α 详尽版)
#
#   S1: 排除前 2 年内死亡 (immortal time bias, mortality only)
#   S2: IR cutoff HOMA-IR >= 2.5 vs >= 3.6 (Asian) vs >= 4.65 (RTH lower)
#   S3: Creatinine-adjusted Phthalate (per 100 mg/dL creatinine)
#   S4: MICE 多重插补 m=20 (mice 包) — IR binary
#   S5: 排除前 5 年慢病诊断 (CHD/CHF/stroke history, MCQ160 系列)
#   S6: 单 PHTHTE metabolite 排除 (leave-one-out, driver detection)
#   S7: 仅 fasting >= 10h (more strict)
#   S8: 仅 NHANES 2011-2018 (post-EPA awareness)
#
# 主暴露: sum_dehp_mol_z (Σ-DEHP) for S1-S5/S7/S8;  S6 = leave-one-out 8 metabolites
# 主结局: IR binary (logistic) + all-cause mortality (Cox, for S1)
# 主结果对比: 主分析 M2 svyglm OR for Σ-DEHP z-score
#
# 输出: output/tables/sensitivity_8sets.csv
# ============================================

suppressPackageStartupMessages({
  library(survey); library(survival); library(dplyr); library(broom); library(purrr)
  library(mice)
})

cat("========================================\n")
cat("009 / 17_sensitivity.R: 8 sensitivity analyses\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

PRIMARY_EXP <- "sum_dehp_mol_z"
cov_base <- c("age","sex_male","race","education","pir","bmi","waist","smoke","hypertension")

# ------------------------------------------------------------------
# 主分析 baseline OR (Σ-DEHP -> IR binary, M2)
# ------------------------------------------------------------------
get_baseline_ir <- function() {
  f <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
  m <- svyglm(f, design = design_main, family = quasibinomial())
  tt <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term == PRIMARY_EXP)
  data.frame(scenario = "Main analysis (baseline)", n = nrow(nhanes_final),
             events = sum(nhanes_final$ir_binary, na.rm=TRUE),
             outcome = "IR_binary", measure = "OR",
             est = tt$estimate, CI_low = tt$conf.low, CI_high = tt$conf.high,
             p_value = tt$p.value, stringsAsFactors = FALSE)
}
baseline <- get_baseline_ir()
cat("--- Baseline (IR binary) ---\n"); print(baseline, row.names=FALSE)

# Generic fit helpers
fit_logistic <- function(data_df, formula, weights_col = "wt_pooled") {
  d <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=as.formula(paste0("~",weights_col)),
                 data=data_df, nest=TRUE)
  svyglm(formula, design=d, family=quasibinomial())
}

extract_or <- function(model, term_name) {
  tt <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term == term_name)
  if (nrow(tt) == 0) return(list(est=NA, lo=NA, hi=NA, p=NA))
  list(est=tt$estimate, lo=tt$conf.low, hi=tt$conf.high, p=tt$p.value)
}

results <- list(baseline)

# ------------------------------------------------------------------
# S1: 排除前 2 年内死亡 (Cox HR)
# ------------------------------------------------------------------
cat("\n--- S1: Exclude deaths within 2 years (immortal time bias) ---\n")
df_s1 <- df_mort %>% dplyr::filter(!is.na(permth), !is.na(mort_allcause))
n_before <- nrow(df_s1)
df_s1 <- df_s1 %>% dplyr::filter(!(mort_allcause == 1 & permth < 24))  # permth in months
n_after <- nrow(df_s1)
cat(sprintf("  Removed %d early deaths\n", n_before - n_after))

res_s1 <- tryCatch({
  d_s1 <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_pooled,
                    data=df_s1, nest=TRUE)
  f <- as.formula(paste("Surv(permth, mort_allcause) ~", PRIMARY_EXP, "+",
                        paste(cov_base, collapse=" + ")))
  m <- svycoxph(f, design=d_s1)
  s <- summary(m)
  ci <- which(rownames(s$coefficients) == PRIMARY_EXP)
  hr <- exp(s$coefficients[ci,"coef"]); se <- s$coefficients[ci,"se(coef)"]
  data.frame(scenario = "S1: Exclude deaths within 2 yrs", n = nrow(df_s1),
             events = sum(df_s1$mort_allcause, na.rm=TRUE),
             outcome = "All_cause_mortality", measure = "HR",
             est = hr, CI_low = exp(s$coefficients[ci,"coef"] - 1.96*se),
             CI_high = exp(s$coefficients[ci,"coef"] + 1.96*se),
             p_value = s$coefficients[ci,"Pr(>|z|)"], stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S1: Exclude deaths within 2 yrs", n=NA, events=NA,
             outcome="All_cause_mortality", measure="HR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
results[[length(results)+1]] <- res_s1
print(res_s1, row.names=FALSE)

# ------------------------------------------------------------------
# S2: IR cutoff sensitivity (>=2.5 / >=3.6 / >=4.65)
# ------------------------------------------------------------------
cat("\n--- S2: IR cutoff sensitivity ---\n")
for (cutoff in c(2.5, 3.6, 4.65)) {
  d <- nhanes_final
  d$ir_var <- as.integer(d$homa_ir >= cutoff)
  ev_n <- sum(d$ir_var, na.rm=TRUE)
  f <- as.formula(paste("ir_var ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
  m <- tryCatch(fit_logistic(d, f), error=function(e) NULL)
  o <- if (!is.null(m)) extract_or(m, PRIMARY_EXP) else list(est=NA,lo=NA,hi=NA,p=NA)
  results[[length(results)+1]] <- data.frame(
    scenario = sprintf("S2: HOMA-IR >= %.2f", cutoff), n = nrow(d), events = ev_n,
    outcome = "IR_binary", measure = "OR",
    est = o$est, CI_low = o$lo, CI_high = o$hi, p_value = o$p,
    stringsAsFactors = FALSE)
  cat(sprintf("  HOMA-IR >= %.2f: OR=%.3f [%.3f, %.3f]  (events=%d)\n",
              cutoff, o$est %||% NA, o$lo %||% NA, o$hi %||% NA, ev_n))
}

# ------------------------------------------------------------------
# S3: Creatinine-adjusted Phthalate
# ------------------------------------------------------------------
cat("\n--- S3: Creatinine-adjusted Sum-DEHP ---\n")
res_s3 <- tryCatch({
  d <- nhanes_final %>% dplyr::filter(!is.na(URXUCR), URXUCR > 0)
  # Build Sum-DEHP creatinine-adjusted (per 100 mg/dL creatinine)
  dehp_imp <- c("URXMHP_imp","URXMHH_imp","URXMOH_imp","URXECP_imp")
  dehp_imp_avail <- intersect(dehp_imp, names(d))
  if (length(dehp_imp_avail) < 3) stop("DEHP imputed cols not enough")
  # mole-weighted sum, then divided by creatinine
  d$sum_dehp_mol_cr <- with(d,
    (URXMHP_imp/278 + URXMHH_imp/294 + URXMOH_imp/292 + URXECP_imp/308) / URXUCR * 100)
  d$sum_dehp_mol_cr_log2 <- log2(pmax(d$sum_dehp_mol_cr, 1e-6))
  d$sum_dehp_mol_cr_z <- as.numeric(scale(d$sum_dehp_mol_cr_log2))
  f <- as.formula(paste("ir_binary ~ sum_dehp_mol_cr_z +", paste(cov_base, collapse=" + ")))
  m <- fit_logistic(d, f)
  o <- extract_or(m, "sum_dehp_mol_cr_z")
  data.frame(scenario="S3: Creatinine-adjusted Sum-DEHP", n=nrow(d),
             events=sum(d$ir_binary, na.rm=TRUE), outcome="IR_binary", measure="OR",
             est=o$est, CI_low=o$lo, CI_high=o$hi, p_value=o$p, stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S3: Creatinine-adjusted Sum-DEHP", n=NA, events=NA,
             outcome="IR_binary", measure="OR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
`%||%` <- function(a,b) if(!is.null(a) && !is.na(a)) a else b
results[[length(results)+1]] <- res_s3
print(res_s3, row.names=FALSE)

# ------------------------------------------------------------------
# S4: MICE m=20 (IR binary)
# ------------------------------------------------------------------
cat("\n--- S4: MICE m=20 multiple imputation ---\n")
res_s4 <- tryCatch({
  # Choose minimal subset to impute (limit dimensionality)
  vars_keep <- c("SEQN","SDMVPSU","SDMVSTRA","wt_pooled","ir_binary",
                 PRIMARY_EXP, cov_base)
  vars_keep <- intersect(vars_keep, names(nhanes_final))
  d4 <- nhanes_final[, vars_keep]
  # Check missing rate
  miss_rate <- sapply(d4, function(x) mean(is.na(x)))
  cat("  Missing rate by var:\n"); print(round(miss_rate, 3))

  set.seed(20260523)
  m_imp <- mice(d4, m = 20, printFlag = FALSE, seed = 20260523)
  # Manual loop over m=20 imputations
  ests <- numeric(20); ses <- numeric(20); ns <- integer(20)
  for (i in 1:20) {
    di <- complete(m_imp, i)
    d_des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_pooled,
                       data=di, nest=TRUE)
    fi <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
    mi <- svyglm(fi, design=d_des, family=quasibinomial())
    tt <- broom::tidy(mi) %>% dplyr::filter(term == PRIMARY_EXP)
    if (nrow(tt) == 0) stop(sprintf("imp %d: term %s not found", i, PRIMARY_EXP))
    ests[i] <- tt$estimate
    ses[i]  <- tt$std.error
    ns[i]   <- nrow(di)
  }
  # Rubin's rules
  q_bar <- mean(ests)
  u_bar <- mean(ses^2)
  b     <- var(ests)
  T_var <- u_bar + (1 + 1/20) * b
  se_pooled <- sqrt(T_var)
  or_pooled <- exp(q_bar)
  ci_lo <- exp(q_bar - 1.96 * se_pooled)
  ci_hi <- exp(q_bar + 1.96 * se_pooled)
  p_val <- 2 * pnorm(-abs(q_bar / se_pooled))
  data.frame(scenario="S4: MICE m=20 multiple imputation", n=round(mean(ns)),
             events=sum(nhanes_final$ir_binary, na.rm=TRUE),
             outcome="IR_binary", measure="OR",
             est=or_pooled, CI_low=ci_lo, CI_high=ci_hi, p_value=p_val,
             stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S4: MICE m=20 multiple imputation", n=NA, events=NA,
             outcome="IR_binary", measure="OR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
results[[length(results)+1]] <- res_s4
print(res_s4, row.names=FALSE)

# ------------------------------------------------------------------
# S5: 排除既往慢病 (CHD/CHF/stroke, MCQ160 系列)
# ------------------------------------------------------------------
cat("\n--- S5: Exclude prior CHD/CHF/stroke (MCQ160 series) ---\n")
res_s5 <- tryCatch({
  d <- nhanes_final
  # MCQ160B (CHF), MCQ160C (CHD), MCQ160E (heart attack), MCQ160F (stroke)
  cvd_cols <- intersect(c("MCQ160B","MCQ160C","MCQ160E","MCQ160F"), names(d))
  if (length(cvd_cols) == 0) stop("No MCQ160 cols available — skip")
  cat(sprintf("  Using cols: %s\n", paste(cvd_cols, collapse=", ")))
  has_cvd <- rep(FALSE, nrow(d))
  for (c in cvd_cols) {
    v <- d[[c]]
    has_cvd <- has_cvd | (!is.na(v) & v == 1)
  }
  n_excl <- sum(has_cvd)
  cat(sprintf("  Excluded %d with prior CVD\n", n_excl))
  d <- d[!has_cvd, ]
  f <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
  m <- fit_logistic(d, f)
  o <- extract_or(m, PRIMARY_EXP)
  data.frame(scenario="S5: Exclude prior CHD/CHF/stroke", n=nrow(d),
             events=sum(d$ir_binary, na.rm=TRUE), outcome="IR_binary", measure="OR",
             est=o$est, CI_low=o$lo, CI_high=o$hi, p_value=o$p, stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S5: Exclude prior CHD/CHF/stroke", n=NA, events=NA,
             outcome="IR_binary", measure="OR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
results[[length(results)+1]] <- res_s5
print(res_s5, row.names=FALSE)

# ------------------------------------------------------------------
# S6: TRUE leave-one-out (drop each metabolite, refit mixture WQS on remaining 7)
# W16 R-DataChain P0 fix: previously this section ran single-metabolite-alone
# (NOT leave-one-out). Manuscript §3.7 claims "removing MiBP nullified the mixture
# signal" — that requires actual leave-one-out. Reimplement with WQS-style
# unweighted mixture index (mean of z-scores of remaining 7 metabolites).
# ------------------------------------------------------------------
cat("\n--- S6: TRUE leave-one-out mixture (drop one metabolite, refit on remaining 7) ---\n")
single_phth <- c(MEP="URXMEP_z", MnBP="URXMBP_z", MiBP="URXMIB_z", MBzP="URXMZP_z",
                 MEHP="URXMHP_z", MEHHP="URXMHH_z", MEOHP="URXMOH_z", MECPP="URXECP_z")
single_phth <- single_phth[single_phth %in% names(nhanes_final)]

# Build mixture from all 8 (full reference)
all_phth_z <- as.character(single_phth)
d_full <- nhanes_final
# Subset to complete cases on all 8
keep_idx <- complete.cases(d_full[, all_phth_z, drop = FALSE])
d_full <- d_full[keep_idx, ]
d_full$mix_full_z <- rowMeans(d_full[, all_phth_z, drop = FALSE], na.rm = TRUE)
f_full <- as.formula(paste("ir_binary ~ mix_full_z +", paste(cov_base, collapse=" + ")))
m_full <- tryCatch(fit_logistic(d_full, f_full), error=function(e) NULL)
o_full <- if (!is.null(m_full)) extract_or(m_full, "mix_full_z") else list(est=NA,lo=NA,hi=NA,p=NA)
results[[length(results)+1]] <- data.frame(
  scenario = "S6: Mixture full (all 8 metab)",
  n = nrow(d_full), events = sum(d_full$ir_binary, na.rm=TRUE),
  outcome = "IR_binary", measure = "OR",
  est = o_full$est, CI_low = o_full$lo, CI_high = o_full$hi, p_value = o_full$p,
  stringsAsFactors = FALSE)
cat(sprintf("  [full mix all 8]: OR=%.3f [%.3f, %.3f] p=%.3g\n",
            o_full$est %||% NA, o_full$lo %||% NA, o_full$hi %||% NA, o_full$p %||% NA))

# Drop each metabolite in turn — refit mixture on remaining 7
for (lab in names(single_phth)) {
  drop_var <- single_phth[[lab]]
  remain_z <- setdiff(all_phth_z, drop_var)
  d <- d_full   # complete-case from above
  d$mix_loo_z <- rowMeans(d[, remain_z, drop = FALSE], na.rm = TRUE)
  f <- as.formula(paste("ir_binary ~ mix_loo_z +", paste(cov_base, collapse=" + ")))
  m <- tryCatch(fit_logistic(d, f), error=function(e) NULL)
  o <- if (!is.null(m)) extract_or(m, "mix_loo_z") else list(est=NA,lo=NA,hi=NA,p=NA)
  results[[length(results)+1]] <- data.frame(
    scenario = sprintf("S6: Mixture LOO (drop %s, 7 remain)", lab),
    n = nrow(d), events = sum(d$ir_binary, na.rm=TRUE),
    outcome = "IR_binary", measure = "OR",
    est = o$est, CI_low = o$lo, CI_high = o$hi, p_value = o$p,
    stringsAsFactors = FALSE)
  cat(sprintf("  drop %-6s: OR=%.3f [%.3f, %.3f] p=%.3g\n",
              lab,
              ifelse(is.null(o$est)||is.na(o$est), NA, o$est),
              ifelse(is.null(o$lo)||is.na(o$lo), NA, o$lo),
              ifelse(is.null(o$hi)||is.na(o$hi), NA, o$hi),
              ifelse(is.null(o$p)||is.na(o$p), NA, o$p)))
}

# ------------------------------------------------------------------
# S7: 仅 fasting >= 10h
# ------------------------------------------------------------------
cat("\n--- S7: Fasting >= 10 h (strict) ---\n")
res_s7 <- tryCatch({
  d <- nhanes_final %>% dplyr::filter(!is.na(fasting_hours), fasting_hours >= 10)
  cat(sprintf("  N after >=10h filter: %d\n", nrow(d)))
  if (nrow(d) < 200) stop("subset too small")
  f <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
  m <- fit_logistic(d, f)
  o <- extract_or(m, PRIMARY_EXP)
  data.frame(scenario="S7: Fasting >= 10h strict", n=nrow(d),
             events=sum(d$ir_binary, na.rm=TRUE), outcome="IR_binary", measure="OR",
             est=o$est, CI_low=o$lo, CI_high=o$hi, p_value=o$p, stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S7: Fasting >= 10h strict", n=NA, events=NA,
             outcome="IR_binary", measure="OR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
results[[length(results)+1]] <- res_s7
print(res_s7, row.names=FALSE)

# ------------------------------------------------------------------
# S8: 仅 NHANES 2011-2018 (post-EPA awareness)
# ------------------------------------------------------------------
cat("\n--- S8: NHANES 2011-2018 only ---\n")
res_s8 <- tryCatch({
  late_cycles <- c("NHANES_2011_2012","NHANES_2013_2014",
                   "NHANES_2015_2016","NHANES_2017_2018")
  if (!"cycle_tag" %in% names(nhanes_final)) stop("cycle_tag not in nhanes_final")
  d <- nhanes_final %>% dplyr::filter(cycle_tag %in% late_cycles)
  cat(sprintf("  N for 2011-2018: %d\n", nrow(d)))
  if (nrow(d) < 200) stop("subset too small")
  f <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "+", paste(cov_base, collapse=" + ")))
  m <- fit_logistic(d, f)
  o <- extract_or(m, PRIMARY_EXP)
  data.frame(scenario="S8: 2011-2018 cycles only", n=nrow(d),
             events=sum(d$ir_binary, na.rm=TRUE), outcome="IR_binary", measure="OR",
             est=o$est, CI_low=o$lo, CI_high=o$hi, p_value=o$p, stringsAsFactors=FALSE)
}, error = function(e) {
  cat("  [err]", conditionMessage(e), "\n")
  data.frame(scenario="S8: 2011-2018 cycles only", n=NA, events=NA,
             outcome="IR_binary", measure="OR",
             est=NA, CI_low=NA, CI_high=NA, p_value=NA, stringsAsFactors=FALSE)
})
results[[length(results)+1]] <- res_s8
print(res_s8, row.names=FALSE)

# ------------------------------------------------------------------
# Combine + save
# ------------------------------------------------------------------
all_sens <- bind_rows(results)
all_sens$est     <- round(all_sens$est, 4)
all_sens$CI_low  <- round(all_sens$CI_low, 4)
all_sens$CI_high <- round(all_sens$CI_high, 4)
all_sens$p_value <- signif(all_sens$p_value, 4)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(all_sens, "output/tables/sensitivity_8sets.csv", row.names = FALSE)

cat("\n========================================\n")
cat("--- Final sensitivity table (8 sets) ---\n")
print(all_sens, row.names=FALSE)
cat(sprintf("\nSaved %d sensitivity rows: output/tables/sensitivity_8sets.csv\n", nrow(all_sens)))
cat("========================================\n")
