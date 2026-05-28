# ============================================
# 009 / scripts/32_retrodesign.R
# Retrospective design analysis: Type S / Type M (Gelman-Carlin 2014)
#
# W16 Round 1 reset — R-Stats M5 + R-Causal MR-5 fix:
# Manuscript §3.8 null mortality finding (CM events n=64-68) is in the
# underpowered zone. Per Gelman-Carlin 2014 *Perspect Psychol Sci*, an under-
# powered null should report Type S (sign error probability) and Type M
# (magnitude exaggeration ratio).
#
# Prior plausible HR per Trasande 2022 (NHANES 2001-2014, all-cause mortality)
# = 1.10-1.15 per SD-z Σ-DEHP. Use HR = 1.10 (log HR ≈ 0.0953) as prior plausible
# effect; SE = se_obs from svycoxph CM model output.
#
# W16 Round 4 SA-B4 extension (R-Bias Ground 18):
# All-cause mortality Σ-DEHP M2 HR 1.69 (n_event=141) is *positive* but
# substantially larger than literature priors (Zhang 2025 HR 1.35 per SD,
# Trasande 2022 HR 1.14 per SD). Extend Type-M analysis to quantify the
# magnitude-exaggeration risk for the all-cause finding.
#  - Zhang 2025 prior: D = log(1.35)
#  - Trasande 2022 prior: D = log(1.14)
#  - s = SE(log HR) derived from cox_mortality_allcause.csv (sum_dehp_mol_z M2)
#  - WARN if Type M ≥ 1.5x for either prior
#
# Output: output/tables/retrodesign_results.csv (CM)
#       + output/tables/retrodesign_all_cause.csv (all-cause, SA-B4)
# ============================================

suppressPackageStartupMessages({
  library(dplyr)
})

cat("========================================\n")
cat("009 / 32_retrodesign.R — Type S / Type M analysis (R-Stats M5 fix)\n")
cat("========================================\n\n")

# Reimplement retro_design (retrodesign package may not be installed)
# Closed-form per Gelman-Carlin 2014 Box 1
retro_design <- function(A, s, alpha = 0.05, df = Inf) {
  # A = true effect size (assumed); s = SE of estimate; alpha; df
  z <- qt(1 - alpha/2, df)
  # Probability that |estimate| > z*s AND has wrong sign
  # = P(estimate < -z*s | true = A) for A > 0
  # Power
  power <- pt(z - A/s, df, lower.tail = FALSE) +
           pt(-z - A/s, df, lower.tail = TRUE)
  # Type S: P(estimate < 0 and significant | true = A)
  type_S <- pt(-z - A/s, df, lower.tail = TRUE) / power
  # Type M: E(|estimate| | |estimate| > z*s) / |A|
  # Approximated via numerical: simulate from t with non-centrality A/s
  set.seed(20260523)
  n_sim <- 1e5
  # For large df, use normal; here z ~ N(A/s, 1) -> estimate ~ N(A, s^2)
  est <- rnorm(n_sim, mean = A, sd = s)
  signif_idx <- abs(est) > z * s
  type_M <- if (sum(signif_idx) > 0) mean(abs(est[signif_idx])) / abs(A) else NA
  list(power = power, type_s = type_S, type_m = type_M)
}

# ----------------------------------------------------------
# Pull observed SE from cox_mortality_cm.csv for Σ-DEHP M2
# ----------------------------------------------------------
cm_df <- tryCatch(read.csv("output/tables/cox_mortality_cm.csv",
                           stringsAsFactors = FALSE),
                  error = function(e) NULL)

dehp_row <- if (!is.null(cm_df)) {
  cm_df[cm_df$exposure == "sum_dehp_mol_z" & cm_df$model == "M2", ]
} else NULL

if (!is.null(dehp_row) && nrow(dehp_row) > 0) {
  hr  <- dehp_row$HR_per_SD[1]
  lo  <- dehp_row$lo[1]
  hi  <- dehp_row$hi[1]
  # SE of log HR
  se_log_hr <- (log(hi) - log(lo)) / (2 * 1.96)
  n_event <- dehp_row$n_event[1]
  cat(sprintf("Observed Σ-DEHP-z CM M2: HR=%.3f (%.3f-%.3f), n_event=%d, SE(log HR)=%.4f\n",
              hr, lo, hi, n_event, se_log_hr))
} else {
  cat("[WARN] CM Σ-DEHP M2 row not found in cox_mortality_cm.csv — using fallback SE 0.27\n")
  se_log_hr <- 0.27
  n_event <- NA
}

# ----------------------------------------------------------
# Prior plausible effects per literature
# ----------------------------------------------------------
prior_hrs <- c(
  "Trasande 2022 NHANES all-cause" = 1.10,
  "Trasande 2022 high estimate"    = 1.15,
  "Zhang 2024 phth CVD est"        = 1.12,
  "Conservative null-ish"          = 1.05
)

retro_rows <- list()
for (lab in names(prior_hrs)) {
  A_log <- log(prior_hrs[[lab]])  # true effect on log HR scale
  rd <- retro_design(A = A_log, s = se_log_hr, alpha = 0.05, df = Inf)
  retro_rows[[length(retro_rows)+1]] <- data.frame(
    prior_label = lab,
    prior_HR    = prior_hrs[[lab]],
    log_HR      = round(A_log, 4),
    SE_obs      = round(se_log_hr, 4),
    power_80pct = round(rd$power, 3),
    type_S      = round(rd$type_s, 3),
    type_M      = round(rd$type_m, 2),
    stringsAsFactors = FALSE
  )
}

retro_df <- do.call(rbind, retro_rows)
retro_df$interpretation <- with(retro_df,
  ifelse(power_80pct >= 0.80,
         "Adequate power",
         ifelse(type_M > 2,
                sprintf("Underpowered — Type M ~%.1fx exaggeration", type_M),
                "Marginal — interpret cautiously")))

cat("\n--- Retro-design results (Gelman-Carlin 2014) ---\n")
print(retro_df, row.names = FALSE)

# Add metadata row
metadata <- data.frame(
  prior_label = "[metadata]",
  prior_HR = NA, log_HR = NA, SE_obs = se_log_hr,
  power_80pct = NA, type_S = NA, type_M = NA,
  interpretation = sprintf("CM mortality n_event=%s; below Peduzzi EPV-10 (=68/9 = 7.6); R-Stats M5 fix",
                           ifelse(is.na(n_event), "?", as.character(n_event))),
  stringsAsFactors = FALSE
)
retro_df_out <- rbind(retro_df, metadata)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(retro_df_out, "output/tables/retrodesign_results.csv", row.names = FALSE)

cat(sprintf("\n[OK] output/tables/retrodesign_results.csv (rows=%d)\n", nrow(retro_df_out)))
cat("\nInterpretation (Gelman-Carlin 2014):\n")
cat("  Power <0.30 + Type M >2  → 'inconclusive', NOT 'no-effect'\n")
cat("  Type S ~0.10              → 10%% probability of sign error if significant\n")

# ============================================================
# W16 Round 4 SA-B4 — All-cause mortality Type-M analysis
# ============================================================
# Pull observed SE from cox_mortality_allcause.csv for Σ-DEHP M2
# (HR 1.69 is the published primary all-cause finding for sum_dehp_mol_z M2)
cat("\n--------------------------------------------------------\n")
cat("W16 Round 4 SA-B4 — All-cause HR 1.69 Type-M analysis\n")
cat("R-Bias Ground 18 extension (Zhang 2025 + Trasande 2022 priors)\n")
cat("--------------------------------------------------------\n\n")

allcause_df <- tryCatch(read.csv("output/tables/cox_mortality_allcause.csv",
                                 stringsAsFactors = FALSE),
                        error = function(e) NULL)

dehp_ac_row <- if (!is.null(allcause_df)) {
  allcause_df[allcause_df$exposure == "sum_dehp_mol_z" &
              allcause_df$model == "M2", ]
} else NULL

if (!is.null(dehp_ac_row) && nrow(dehp_ac_row) > 0) {
  hr_ac  <- dehp_ac_row$HR_per_SD[1]
  lo_ac  <- dehp_ac_row$lo[1]
  hi_ac  <- dehp_ac_row$hi[1]
  se_log_hr_ac <- (log(hi_ac) - log(lo_ac)) / (2 * 1.96)
  n_event_ac   <- dehp_ac_row$n_event[1]
  cat(sprintf("Observed Σ-DEHP-z all-cause M2: HR=%.3f (%.3f-%.3f), n_event=%d, SE(log HR)=%.4f\n",
              hr_ac, lo_ac, hi_ac, n_event_ac, se_log_hr_ac))
} else {
  cat("[WARN] All-cause Σ-DEHP M2 row not found in cox_mortality_allcause.csv — using fallback SE 0.12\n")
  se_log_hr_ac <- 0.12
  n_event_ac   <- NA
  hr_ac        <- NA
}

# Two literature priors per R-Bias Ground 18
prior_hrs_ac <- c(
  "Zhang 2025 LancetPH (HR 1.35, N=8378)"    = 1.35,
  "Trasande 2022 EnvironPollut (HR 1.14, N=5303)" = 1.14
)

retro_ac_rows <- list()
for (lab in names(prior_hrs_ac)) {
  A_log <- log(prior_hrs_ac[[lab]])
  rd <- retro_design(A = A_log, s = se_log_hr_ac, alpha = 0.05, df = Inf)
  # sign-error-rate is the raw P(wrong sign | true = A) without conditioning on
  # significance; equals pt(-A/s, df). Manuscript-friendly companion to Type S
  # (which is conditional on significance).
  sign_err_rate <- pt(-A_log / se_log_hr_ac, df = Inf, lower.tail = TRUE)
  retro_ac_rows[[length(retro_ac_rows)+1]] <- data.frame(
    prior_label       = lab,
    prior_D_log_HR    = round(A_log, 4),
    prior_HR          = prior_hrs_ac[[lab]],
    SE_obs            = round(se_log_hr_ac, 4),
    power             = round(rd$power, 3),
    type_S            = round(rd$type_s, 4),
    type_M            = round(rd$type_m, 2),
    sign_error_rate   = round(sign_err_rate, 4),
    exaggeration_ratio = round(rd$type_m, 2),
    stringsAsFactors  = FALSE
  )
}
retro_ac_df <- do.call(rbind, retro_ac_rows)

# Trigger warning if Type M >= 1.5x for either prior
any_exag <- any(retro_ac_df$type_M >= 1.5, na.rm = TRUE)
retro_ac_df$warning <- ifelse(
  retro_ac_df$type_M >= 1.5,
  sprintf("All-cause HR 1.69 Type-M exaggeration concern (Type M ~%.1fx)",
          retro_ac_df$type_M),
  "OK (Type M < 1.5x)"
)

cat("\n--- All-cause retro-design results (Gelman-Carlin 2014) ---\n")
print(retro_ac_df, row.names = FALSE)

# Add metadata footer row
metadata_ac <- data.frame(
  prior_label       = "[metadata]",
  prior_D_log_HR    = NA,
  prior_HR          = NA,
  SE_obs            = se_log_hr_ac,
  power             = NA,
  type_S            = NA,
  type_M            = NA,
  sign_error_rate   = NA,
  exaggeration_ratio = NA,
  warning = sprintf("All-cause n_event=%s; observed HR=%s; W16 R4 SA-B4 R-Bias Ground 18 extension",
                    ifelse(is.na(n_event_ac), "?", as.character(n_event_ac)),
                    ifelse(is.na(hr_ac), "?", sprintf("%.3f", hr_ac))),
  stringsAsFactors  = FALSE
)
retro_ac_df_out <- rbind(retro_ac_df, metadata_ac)

write.csv(retro_ac_df_out, "output/tables/retrodesign_all_cause.csv",
          row.names = FALSE)
cat(sprintf("\n[OK] output/tables/retrodesign_all_cause.csv (rows=%d)\n",
            nrow(retro_ac_df_out)))

if (any_exag) {
  cat("\n[WARN] All-cause HR 1.69 Type-M exaggeration concern\n")
  cat(sprintf("       Zhang 2025 prior (HR 1.35): Type M = %.2fx\n",
              retro_ac_df$type_M[1]))
  cat(sprintf("       Trasande 2022 prior (HR 1.14): Type M = %.2fx\n",
              retro_ac_df$type_M[2]))
  cat("       (Type M >= 1.5x indicates magnitude exaggeration risk;\n")
  cat("        observed HR 1.69 may overstate the underlying effect.)\n")
} else {
  cat("\n[OK] Both priors yield Type M < 1.5x — no exaggeration concern.\n")
}

cat("\n========================================\n")
cat("Retrodesign done.\n")
cat("========================================\n")
