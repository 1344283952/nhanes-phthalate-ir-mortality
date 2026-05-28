# ============================================
# 009 / scripts/10_qgcomp.R
# Quantile g-computation (Keil 2020 BES, DOI 10.1289/EHP5838)
# for Phthalate mixture × IR outcome (HOMA-IR continuous + IR binary)
#
# Pipeline (Keil 2020 §3 recommendation):
#   - qgcomp.glm.noboot → component weights (pos/neg split)
#   - qgcomp.glm.boot   → total mixture effect ψ + bootstrap CI (B=500)
#
# 输入: data/processed/nhanes_final.RData
# 输出: output/tables/qgcomp_phth_homa.csv (HOMA-IR continuous)
#       output/tables/qgcomp_phth_ir_binary.csv (IR binary)
#       output/tables/qgcomp_phth_results.RData (raw fit objects)
#       output/tables/qgcomp_phth_weights.csv (pos+neg direction split from noboot)
#       output/tables/qgcomp_phth_summary.csv (ψ + CI compact)
#
# Mixture (8 Phthalate metabolite z-scores):
#   URXMEP_z (MEP, DEP - cosmetics)
#   URXMBP_z (MnBP, DBP - PVC)
#   URXMIB_z (MiBP - PVC)
#   URXMZP_z (MBzP, BBzP - flooring)
#   URXMHP_z (MEHP, DEHP M1)
#   URXMHH_z (MEHHP, DEHP M2)
#   URXMOH_z (MEOHP, DEHP M3)
#   URXECP_z (MECPP, DEHP M4)
#
# Covariates: age, sex_male, race, education, pir, bmi, waist,
#             smoke_ever (binary), cotinine_log
# Note: 'drink' completely missing (0/2239) per 03_clean — excluded.
#       cotinine_log 1095/2239 → median imputed to keep N stable.
#
# Weights: wt_pooled (NHANES pooled fasting-subsample weight per Series 2 No. 190)
# Bootstrap: B = 500 (Keil 2020 §3.2 recommends ≥ 500 for percentile CI)
# Set seed 20260524 (009 启动 + 1d)
# ============================================

set.seed(20260524)

suppressPackageStartupMessages({
  library(dplyr)
  library(qgcomp)
})

cat("========================================\n")
cat("009 / 10_qgcomp — Phthalate mixture × IR (continuous + binary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final loaded: n=%d ; IR cases=%d (%.1f%%)\n",
            nrow(nhanes_final),
            sum(nhanes_final$ir_binary == 1, na.rm = TRUE),
            100 * mean(nhanes_final$ir_binary == 1, na.rm = TRUE)))

# ---------------------------------------------------------------
# Step 1: Build mixture & covariate set
# ---------------------------------------------------------------
mixture <- c("URXMEP_z", "URXMBP_z", "URXMIB_z", "URXMZP_z",
             "URXMHP_z", "URXMHH_z", "URXMOH_z", "URXECP_z")
cov_set <- c("age", "sex_male", "race", "education", "pir",
             "bmi", "waist", "smoke_ever", "cotinine_log")

# Construct smoke_ever (binary 0/1) from existing smoke factor
if (!"smoke_ever" %in% names(nhanes_final)) {
  nhanes_final$smoke_ever <- ifelse(!is.na(nhanes_final$smoke) &
                                      nhanes_final$smoke == "Ever", 1L, 0L)
}

# Median-impute cotinine_log (keep N=2,239)
if ("cotinine_log" %in% names(nhanes_final)) {
  median_cot <- median(nhanes_final$cotinine_log, na.rm = TRUE)
  nhanes_final$cotinine_log <- ifelse(is.na(nhanes_final$cotinine_log),
                                      median_cot, nhanes_final$cotinine_log)
}

core_keep <- c("SEQN", mixture, cov_set, "homa_ir_log", "ir_binary",
               "wt_pooled", "SDMVPSU", "SDMVSTRA")
df_q <- nhanes_final %>%
  select(any_of(core_keep)) %>%
  filter(if_all(all_of(c(mixture, "age", "sex_male", "race", "education",
                          "pir", "bmi", "homa_ir_log", "ir_binary", "wt_pooled")),
                ~ !is.na(.)))
for (cv in c("waist")) {
  df_q[[cv]][is.na(df_q[[cv]])] <- median(df_q[[cv]], na.rm = TRUE)
}
df_q$smoke_ever[is.na(df_q$smoke_ever)] <- 0L

cat(sprintf("Analytic n (post complete-case): %d ; IR cases = %d (%.1f%%)\n",
            nrow(df_q),
            sum(df_q$ir_binary == 1),
            100 * mean(df_q$ir_binary == 1)))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

# ---------------------------------------------------------------
# Helper: extract ψ + CI + p from qgcomp.*.boot output
# qgcomp boot output structure:
#   $psi      length 1 — mixture effect on link scale
#   $ci       length 2 — c(lcl, ucl)
#   $pval     length 2 — c(intercept_p, psi_p)
#   $var.psi  length 1
#   $coef     length 2 — c(intercept, psi)
#   $B, $n    NULL in boot variant — get from qg$bootsamps + nrow(qg$qx)
# ---------------------------------------------------------------
.safe_num <- function(x, i = 1L) {
  if (is.null(x) || length(x) < i) NA_real_ else as.numeric(x)[i]
}

summarize_qg <- function(qg, outcome_lbl, scale_lbl, n_obs) {
  if (is.null(qg)) {
    return(data.frame(outcome = outcome_lbl, psi = NA_real_, psi_se = NA_real_,
                      psi_lcl = NA_real_, psi_ucl = NA_real_, p = NA_real_,
                      scale = scale_lbl, boot_B = NA_integer_, n_obs = n_obs,
                      stringsAsFactors = FALSE))
  }
  psi <- .safe_num(qg$psi, 1)
  ci  <- if (!is.null(qg$ci)) as.numeric(qg$ci) else c(NA_real_, NA_real_)
  pv  <- .safe_num(qg$pval, 2)  # qgcomp boot returns c(intercept_p, psi_p)
  se  <- if (!is.null(qg$var.psi)) sqrt(.safe_num(qg$var.psi, 1)) else NA_real_
  data.frame(
    outcome = outcome_lbl,
    psi = psi, psi_se = se,
    psi_lcl = if (length(ci) >= 2) ci[1] else NA_real_,
    psi_ucl = if (length(ci) >= 2) ci[2] else NA_real_,
    p = pv,
    scale = scale_lbl,
    boot_B = if (!is.null(qg$bootsamps)) nrow(qg$bootsamps) else NA_integer_,
    n_obs = n_obs,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------
# Helper: extract pos/neg weights from qgcomp.*.noboot
# noboot variant has $pos.weights / $neg.weights as named vectors
# ---------------------------------------------------------------
extract_weights <- function(qg_nb, outcome_lbl) {
  if (is.null(qg_nb)) return(data.frame())
  pw <- if (!is.null(qg_nb$pos.weights)) qg_nb$pos.weights else numeric(0)
  nw <- if (!is.null(qg_nb$neg.weights)) qg_nb$neg.weights else numeric(0)
  parts <- list()
  if (length(pw) > 0) {
    parts[[length(parts) + 1]] <- data.frame(
      outcome = outcome_lbl, direction = "positive",
      metabolite = names(pw), weight = as.numeric(pw),
      stringsAsFactors = FALSE)
  }
  if (length(nw) > 0) {
    parts[[length(parts) + 1]] <- data.frame(
      outcome = outcome_lbl, direction = "negative",
      metabolite = names(nw), weight = as.numeric(nw),
      stringsAsFactors = FALSE)
  }
  if (length(parts) == 0) return(data.frame())
  do.call(rbind, parts)
}

write_outcome_csv <- function(qg_boot, qg_nb, outcome_lbl, scale_lbl, csv_path, n_obs) {
  psi_row <- summarize_qg(qg_boot, outcome_lbl, scale_lbl, n_obs)
  if (scale_lbl == "OR") {
    psi_row$psi_exp <- exp(psi_row$psi)
    psi_row$psi_exp_lcl <- exp(psi_row$psi_lcl)
    psi_row$psi_exp_ucl <- exp(psi_row$psi_ucl)
    psi_row$effect_str <- sprintf("OR=%.3f (%.3f-%.3f), p=%.4g",
                                  psi_row$psi_exp, psi_row$psi_exp_lcl,
                                  psi_row$psi_exp_ucl, psi_row$p)
  } else {
    psi_row$effect_str <- sprintf("psi=%.3f (%.3f to %.3f), p=%.4g",
                                  psi_row$psi, psi_row$psi_lcl,
                                  psi_row$psi_ucl, psi_row$p)
  }
  # Append top 3 pos / neg from noboot for clarity
  if (!is.null(qg_nb)) {
    pw <- qg_nb$pos.weights; nw <- qg_nb$neg.weights
    pw_str <- if (length(pw) > 0)
      paste(head(sprintf("%s=%.2f", names(sort(pw, decreasing = TRUE)),
                          sort(pw, decreasing = TRUE)), 3), collapse = "; ")
      else ""
    nw_str <- if (length(nw) > 0)
      paste(head(sprintf("%s=%.2f", names(sort(nw, decreasing = TRUE)),
                          sort(nw, decreasing = TRUE)), 3), collapse = "; ")
      else ""
    psi_row$top_pos_weights <- pw_str
    psi_row$top_neg_weights <- nw_str
  }
  write.csv(psi_row, csv_path, row.names = FALSE)
  cat(sprintf("  → %s\n", csv_path))
  cat(sprintf("  %s\n", psi_row$effect_str))
  if (!is.null(qg_nb)) {
    cat(sprintf("  top pos weights: %s\n", psi_row$top_pos_weights))
    cat(sprintf("  top neg weights: %s\n", psi_row$top_neg_weights))
  }
}

# ---------------------------------------------------------------
# Step 2: HOMA-IR continuous (gaussian)
# ---------------------------------------------------------------
cat("\n[1a/2] qgcomp.glm.noboot (weights) → HOMA-IR log ...\n")
f_homa <- as.formula(paste("homa_ir_log ~",
                            paste(c(mixture, cov_set), collapse = " + ")))
qg_homa_nb <- tryCatch({
  qgcomp::qgcomp.glm.noboot(
    f = f_homa, expnms = mixture, data = df_q,
    family = gaussian(), q = 4,
    weights = df_q$wt_pooled
  )
}, error = function(e) {
  cat("ERROR HOMA noboot:", conditionMessage(e), "\n"); NULL
})

cat("[1b/2] qgcomp.glm.boot (psi + CI) → HOMA-IR log (B=500) ...\n")
t0 <- Sys.time()
qg_homa <- tryCatch({
  qgcomp::qgcomp.glm.boot(
    f = f_homa, expnms = mixture, data = df_q,
    family = gaussian(), q = 4, B = 500, seed = 20260524,
    weights = df_q$wt_pooled
  )
}, error = function(e) {
  cat("ERROR HOMA boot:", conditionMessage(e), "\n"); NULL
})
cat(sprintf("  HOMA boot fit time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------
# Step 3: IR binary (logistic, OR scale)
# ---------------------------------------------------------------
cat("\n[2a/2] qgcomp.glm.noboot (weights) → IR binary ...\n")
f_ir <- as.formula(paste("ir_binary ~",
                          paste(c(mixture, cov_set), collapse = " + ")))
qg_ir_nb <- tryCatch({
  qgcomp::qgcomp.glm.noboot(
    f = f_ir, expnms = mixture, data = df_q,
    family = binomial(link = "logit"), q = 4,
    weights = df_q$wt_pooled
  )
}, error = function(e) {
  cat("ERROR IR-binary noboot:", conditionMessage(e), "\n"); NULL
})

cat("[2b/2] qgcomp.glm.boot (OR + CI) → IR binary (B=500) ...\n")
t0 <- Sys.time()
qg_ir <- tryCatch({
  qgcomp::qgcomp.glm.boot(
    f = f_ir, expnms = mixture, data = df_q,
    family = binomial(link = "logit"), q = 4, B = 500, seed = 20260524,
    weights = df_q$wt_pooled
  )
}, error = function(e) {
  cat("ERROR IR-binary boot:", conditionMessage(e), "\n"); NULL
})
cat(sprintf("  IR boot fit time: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------
# Step 4: Write per-outcome CSVs
# ---------------------------------------------------------------
write_outcome_csv(qg_homa, qg_homa_nb, "HOMA-IR (log, continuous)",
                  "beta (log-scale)",
                  "output/tables/qgcomp_phth_homa.csv", nrow(df_q))
write_outcome_csv(qg_ir, qg_ir_nb, "IR binary (HOMA >= 2.5)", "OR",
                  "output/tables/qgcomp_phth_ir_binary.csv", nrow(df_q))

# ---------------------------------------------------------------
# Step 5: Weights long-format CSV (from noboot)
# ---------------------------------------------------------------
weights_df <- rbind(
  extract_weights(qg_homa_nb, "HOMA-IR (log)"),
  extract_weights(qg_ir_nb,   "IR binary (logit)")
)
if (nrow(weights_df) > 0) {
  write.csv(weights_df, "output/tables/qgcomp_phth_weights.csv", row.names = FALSE)
  cat("\n→ output/tables/qgcomp_phth_weights.csv\n")
  print(weights_df)
}

# ---------------------------------------------------------------
# Step 6: Combined compact summary
# ---------------------------------------------------------------
summary_df <- rbind(
  summarize_qg(qg_homa, "HOMA-IR (log, gaussian)", "beta (log-scale)", nrow(df_q)),
  summarize_qg(qg_ir,   "IR binary (HOMA >= 2.5)", "log-OR", nrow(df_q))
)
write.csv(summary_df, "output/tables/qgcomp_phth_summary.csv", row.names = FALSE)
cat("\n→ output/tables/qgcomp_phth_summary.csv\n")
print(summary_df)

save(qg_homa, qg_homa_nb, qg_ir, qg_ir_nb,
     weights_df, summary_df,
     file = "output/tables/qgcomp_phth_results.RData")

cat("\n保存:\n")
cat("  output/tables/qgcomp_phth_homa.csv\n")
cat("  output/tables/qgcomp_phth_ir_binary.csv\n")
cat("  output/tables/qgcomp_phth_weights.csv\n")
cat("  output/tables/qgcomp_phth_summary.csv\n")
cat("  output/tables/qgcomp_phth_results.RData\n")
cat("\nDONE 10_qgcomp.R\n")
