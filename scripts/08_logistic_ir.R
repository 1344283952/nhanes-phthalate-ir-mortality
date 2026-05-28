# ============================================
# 009 / scripts/08_logistic_ir.R
# Logistic IR binary (HOMA≥2.5) primary outcome + Linear HOMA-IR continuous (log)
# Stack 1 主分析 (N=2,239)
# 输入: data/processed/nhanes_final.RData + nhanes_design.RData
# 输出: output/tables/logistic_ir_binary.csv
#       output/tables/logistic_ir_binary_asian.csv   (sensitivity HOMA≥3.6)
#       output/tables/linear_homa_ir.csv             (secondary, log-HOMA continuous)
#
# 3 model 递进调整: Crude / M1 (age+sex+race) / M2 (M1 + edu+pir+bmi+smoke+hypertension)
# 各 Phthalate 单独 (8 metabolites + Σ-DEHP/HMW/LMW)
# OR (95% CI) per-SD + 四分位 OR + P-trend
#
# 用 survey::svyglm (binomial quasi-binomial for IR binary, gaussian for log-HOMA)
# ============================================

suppressPackageStartupMessages({
  library(survey); library(dplyr); library(broom); library(purrr)
})

cat("========================================\n")
cat("009 / Logistic IR (Stack 1, N=2,239)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

# Re-derive drink fallback (same as 07)
if (all(is.na(nhanes_final$drink)) && "ALQ101" %in% names(nhanes_final)) {
  nhanes_final$drink <- factor(case_when(
    nhanes_final$ALQ101 == 1 ~ "Yes",
    nhanes_final$ALQ101 == 2 ~ "No"
  ), levels = c("No","Yes"))
  cat(sprintf("[fallback] drink ← ALQ101: non-NA = %d\n",
              sum(!is.na(nhanes_final$drink))))
}
use_drink <- !all(is.na(nhanes_final$drink)) && length(unique(na.omit(nhanes_final$drink))) > 1

# Cast factors
for (v in c("race","education","smoke","drink","sex_male","hypertension")) {
  if (v %in% names(nhanes_final) && !is.factor(nhanes_final[[v]])) {
    nhanes_final[[v]] <- factor(nhanes_final[[v]])
  }
}

# Rebuild design with mutations
design_main <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

# ----------------------------------------------------------
# Exposure list
# ----------------------------------------------------------
phth_exposures <- c(
  "URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
  "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z",
  "sum_dehp_mol_z","sum_hmw_z","sum_lmw_z"
)
phth_exposures <- intersect(phth_exposures, names(nhanes_final))

cat(sprintf("Stack 1 main: N=%d, IR binary cases=%d (%.1f%%)\n",
            nrow(nhanes_final),
            sum(nhanes_final$ir_binary == 1, na.rm = TRUE),
            100 * mean(nhanes_final$ir_binary == 1, na.rm = TRUE)))
cat(sprintf("IR binary Asian cutoff (>=3.6) cases=%d (%.1f%%)\n",
            sum(nhanes_final$ir_binary_asian == 1, na.rm = TRUE),
            100 * mean(nhanes_final$ir_binary_asian == 1, na.rm = TRUE)))

# Adjustment sets
# W16 R-Stats C7 fix — Path A (text fix to match code):
#   M2 includes BMI + hypertension as "mediator-confounders" (sensitivity over-adjustment
#   block per Westreich-Greenland 2013); this is the EXACT controlled-direct-effect (CDE)
#   block. The Pearl-backdoor MSAS (Methods §2.4 primary) drops BMI/waist/hypertension
#   and is recovered through the CMAverse 4-way decomposition (12_cmaverse_4way.R Rte).
#   The headline single-spec OR 1.065 (P=0.262 / 0.334) is from this CDE-style M2,
#   intentionally — the 47.7% mediation finding (CMAverse) depends on this contrast.
#   kcal/fish/drink omission: kcal/fish in §S3.1 of Pearl-backdoor are reserved for
#   the cov_pre Pearl-backdoor sensitivity block; drink is empirically all-NA across
#   the D-J cohort (sub-agent A 2026-05-23 实查 ALQ111 全 NA 教训).
cov_M1 <- "age + sex_male + race"
cov_M2_base <- "age + sex_male + race + education + pir + bmi + smoke + hypertension"
cov_M2 <- if (use_drink) paste(cov_M2_base, "+ drink") else cov_M2_base
if (!use_drink) cat("[note] drink 排除自 M2 (全 NA / no variation)\n")

specs <- list(
  Crude = "1",
  M1    = cov_M1,
  M2    = cov_M2
)

# ----------------------------------------------------------
# Helper: 跑一个 phth × model × outcome
# fam: "binomial" → svyglm quasibinomial, OR per-SD
#      "gaussian" → svyglm gaussian, beta per-SD
# ----------------------------------------------------------
run_logit <- function(design, exposure, outcome, model_label, adj, fam = "binomial") {
  fml_str <- if (adj == "1") sprintf("%s ~ %s", outcome, exposure) else
                             sprintf("%s ~ %s + %s", outcome, exposure, adj)
  fam_obj <- if (fam == "binomial") quasibinomial() else gaussian()
  fit <- tryCatch(svyglm(as.formula(fml_str), design = design, family = fam_obj),
                  error = function(e) {
                    cat(sprintf("[err] %s/%s/%s: %s\n",
                                exposure, outcome, model_label, conditionMessage(e)))
                    NULL
                  })
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  if (!exposure %in% rownames(s$coefficients)) return(NULL)
  beta <- s$coefficients[exposure, "Estimate"]
  se   <- s$coefficients[exposure, "Std. Error"]
  p    <- s$coefficients[exposure, "Pr(>|t|)"]
  if (fam == "binomial") {
    data.frame(
      exposure = exposure, outcome = outcome, model = model_label,
      OR_per_SD = exp(beta),
      lo        = exp(beta - 1.96 * se),
      hi        = exp(beta + 1.96 * se),
      p_per_SD  = p,
      n_obs     = sum(weights(fit, type = "prior") > 0),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      exposure = exposure, outcome = outcome, model = model_label,
      beta_per_SD = beta,
      lo          = beta - 1.96 * se,
      hi          = beta + 1.96 * se,
      p_per_SD    = p,
      n_obs       = sum(weights(fit, type = "prior") > 0),
      stringsAsFactors = FALSE
    )
  }
}

run_logit_qtrend <- function(design, exposure_raw, outcome, model_label, adj, fam = "binomial") {
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

  fam_obj <- if (fam == "binomial") quasibinomial() else gaussian()

  # P-trend
  fml_t <- if (adj == "1")
    sprintf("%s ~ %s", outcome, paste0(exposure_raw,"_qmed")) else
    sprintf("%s ~ %s + %s", outcome, paste0(exposure_raw,"_qmed"), adj)
  fit_t <- tryCatch(svyglm(as.formula(fml_t), design = design, family = fam_obj),
                    error = function(e) NULL)
  qmed_var <- paste0(exposure_raw,"_qmed")
  p_trend <- if (!is.null(fit_t) && qmed_var %in% rownames(summary(fit_t)$coefficients))
    summary(fit_t)$coefficients[qmed_var, "Pr(>|t|)"] else NA_real_

  # Q-level OR (Q2/3/4 vs Q1)
  fml_q <- if (adj == "1")
    sprintf("%s ~ %s", outcome, paste0(exposure_raw,"_q")) else
    sprintf("%s ~ %s + %s", outcome, paste0(exposure_raw,"_q"), adj)
  fit_q <- tryCatch(svyglm(as.formula(fml_q), design = design, family = fam_obj),
                    error = function(e) NULL)
  q_rows <- NULL
  if (!is.null(fit_q)) {
    sq <- summary(fit_q)$coefficients
    for (qlab in c("Q2","Q3","Q4")) {
      term <- paste0(exposure_raw, "_q", qlab)
      if (term %in% rownames(sq)) {
        beta <- sq[term, "Estimate"]; se <- sq[term, "Std. Error"]; p <- sq[term, "Pr(>|t|)"]
        if (fam == "binomial") {
          q_rows <- rbind(q_rows, data.frame(
            exposure = exposure_raw, outcome = outcome, model = model_label,
            term = paste0("q_",qlab,"_vs_Q1"),
            OR = exp(beta), lo = exp(beta - 1.96*se), hi = exp(beta + 1.96*se), p = p,
            stringsAsFactors = FALSE))
        } else {
          q_rows <- rbind(q_rows, data.frame(
            exposure = exposure_raw, outcome = outcome, model = model_label,
            term = paste0("q_",qlab,"_vs_Q1"),
            beta = beta, lo = beta - 1.96*se, hi = beta + 1.96*se, p = p,
            stringsAsFactors = FALSE))
        }
      }
    }
  }
  if (!is.null(q_rows)) q_rows$p_trend <- p_trend
  q_rows
}

# ----------------------------------------------------------
# Run: IR binary (HOMA>=2.5) + IR binary Asian (HOMA>=3.6) + log-HOMA continuous
# ----------------------------------------------------------
outcomes_binom <- c("ir_binary","ir_binary_asian")
outcomes_cont  <- c("homa_ir_log")

# Binom per-SD
res_b <- list()
for (oc in outcomes_binom) for (exp in phth_exposures) for (m in names(specs)) {
  r <- run_logit(design_main, exp, oc, m, specs[[m]], "binomial")
  if (!is.null(r)) res_b[[length(res_b)+1]] <- r
}
res_b_df <- do.call(rbind, res_b)

# Binom quartile
res_bq <- list()
for (oc in outcomes_binom) for (exp in phth_exposures) for (m in names(specs)) {
  r <- run_logit_qtrend(design_main, exp, oc, m, specs[[m]], "binomial")
  if (!is.null(r)) res_bq[[length(res_bq)+1]] <- r
}
res_bq_df <- do.call(rbind, res_bq)

# Linear log-HOMA per-SD
res_l <- list()
for (oc in outcomes_cont) for (exp in phth_exposures) for (m in names(specs)) {
  r <- run_logit(design_main, exp, oc, m, specs[[m]], "gaussian")
  if (!is.null(r)) res_l[[length(res_l)+1]] <- r
}
res_l_df <- do.call(rbind, res_l)

# Linear log-HOMA quartile
res_lq <- list()
for (oc in outcomes_cont) for (exp in phth_exposures) for (m in names(specs)) {
  r <- run_logit_qtrend(design_main, exp, oc, m, specs[[m]], "gaussian")
  if (!is.null(r)) res_lq[[length(res_lq)+1]] <- r
}
res_lq_df <- do.call(rbind, res_lq)

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------
out_ir         <- res_b_df  %>% filter(outcome == "ir_binary")
out_ir_asian   <- res_b_df  %>% filter(outcome == "ir_binary_asian")
out_ir_q       <- res_bq_df %>% filter(outcome == "ir_binary")
out_ir_asian_q <- res_bq_df %>% filter(outcome == "ir_binary_asian")
out_homa       <- res_l_df
out_homa_q     <- res_lq_df

cat("\n--- IR binary (HOMA>=2.5) per-SD preview ---\n"); print(out_ir, row.names = FALSE)
cat("\n--- IR binary (HOMA>=2.5) quartile preview ---\n"); print(out_ir_q, row.names = FALSE)
cat("\n--- log-HOMA-IR per-SD preview ---\n"); print(out_homa, row.names = FALSE)
cat("\n--- IR binary Asian (HOMA>=3.6) per-SD preview ---\n"); print(out_ir_asian, row.names = FALSE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(out_ir,         "output/tables/logistic_ir_binary.csv",          row.names = FALSE)
write.csv(out_ir_q,       "output/tables/logistic_ir_binary_quartile.csv", row.names = FALSE)
write.csv(out_ir_asian,   "output/tables/logistic_ir_binary_asian.csv",    row.names = FALSE)
write.csv(out_ir_asian_q, "output/tables/logistic_ir_binary_asian_quartile.csv", row.names = FALSE)
write.csv(out_homa,       "output/tables/linear_homa_ir.csv",              row.names = FALSE)
write.csv(out_homa_q,     "output/tables/linear_homa_ir_quartile.csv",     row.names = FALSE)

cat("\n[OK] output/tables/logistic_ir_binary.csv\n")
cat("[OK] output/tables/logistic_ir_binary_quartile.csv\n")
cat("[OK] output/tables/logistic_ir_binary_asian.csv\n")
cat("[OK] output/tables/logistic_ir_binary_asian_quartile.csv\n")
cat("[OK] output/tables/linear_homa_ir.csv\n")
cat("[OK] output/tables/linear_homa_ir_quartile.csv\n")

cat("\n========================================\n")
cat("Logistic done.\n")
cat("========================================\n")
