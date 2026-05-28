# ============================================
# 009 / scripts/27_lag_analysis.R
# Reverse-causation evaluation: Lag analysis + Negative control outcome
# Phthalate semi-lives are short (12-48 h), but chronic exposure is the proxy.
#
# Two complementary tests:
#
# (1) Lag-style sensitivity (proxy for short-term reverse causation):
#     Exclude participants with self-reported recent (<= 1 year) changes
#     in alcohol intake (ALQ151 / ALQ160) or dietary intake
#     (DR1TKCAL extreme < 5th or > 95th centile, proxy for unstable diet).
#     Re-estimate Sigma-DEHP -> IR association in the restricted set.
#
# (2) Negative control outcome (Lipsitch 2010):
#     Phthalate exposure should be unrelated to:
#       LBXHGB    (hemoglobin)
#       LBXRBCSI  (red blood cell count)
#       LBXMCVSI  (mean corpuscular volume)
#       LBXPLTSI  (platelet count)
#     None of these is biologically downstream of phthalate -> IR.
#     A non-null association raises concern about residual confounding.
#
# Output:
#   output/tables/lag_negative_control.csv
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(survey); library(broom); library(tidyr); library(ggplot2)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 27 Lag + Negative Control Outcome\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

df <- nhanes_final

# Re-derive Q4 vs Q1-Q3 indicator
q_breaks <- quantile(df$sum_dehp_mol, probs = seq(0, 1, 0.25), na.rm = TRUE)
df$dehp_q <- cut(df$sum_dehp_mol, breaks = q_breaks, include.lowest = TRUE,
                 labels = c("Q1","Q2","Q3","Q4"))
df$high_dehp <- as.integer(df$dehp_q == "Q4")

# ============================================
# (1) Lag-restricted sample (exclude unstable diet/alcohol cases)
# ============================================
cat("\n[1] Lag-restricted re-estimation (stable diet + alcohol proxy)...\n")

# Build "unstable" flag.  We don't have direct repeat measures, but use:
#   - extreme kcal_day (< P5 or > P95 within stack) -> dietary instability proxy
#   - extreme BMI (z-score abs > 3) -> chronic vs acute weight change proxy
# (alcohol variable is all-NA in this cohort, so we cannot use it.)
if ("kcal_day" %in% names(df)) {
  kcal_lo <- quantile(df$kcal_day, 0.05, na.rm = TRUE)
  kcal_hi <- quantile(df$kcal_day, 0.95, na.rm = TRUE)
  df$kcal_unstable <- (!is.na(df$kcal_day)) & (df$kcal_day < kcal_lo | df$kcal_day > kcal_hi)
} else {
  df$kcal_unstable <- FALSE
}

# bmi-based proxy: extreme BMI z (top/bottom 2%) often indicates short-term change
bmi_lo <- quantile(df$bmi, 0.02, na.rm = TRUE)
bmi_hi <- quantile(df$bmi, 0.98, na.rm = TRUE)
df$bmi_extreme <- (!is.na(df$bmi)) & (df$bmi < bmi_lo | df$bmi > bmi_hi)

# Subset 1: lag-restricted (stable diet AND non-extreme BMI)
df_lag <- df %>% filter(!kcal_unstable, !bmi_extreme)

cat(sprintf("Original analytic n     : %d\n", nrow(df)))
cat(sprintf("Lag-restricted n        : %d\n", nrow(df_lag)))

# Re-estimate main effect
re_estimate_main <- function(dat, label) {
  needed <- c("ir_binary","high_dehp","age","RIAGENDR","race","education",
              "pir","bmi","waist","smoke","SDMVPSU","SDMVSTRA","wt_pooled")
  dat_cc <- dat[complete.cases(dat[, needed]), ]
  # Drop any unused factor levels to avoid contrasts error
  for (v in c("race","education","smoke")) {
    if (is.factor(dat_cc[[v]])) dat_cc[[v]] <- droplevels(dat_cc[[v]])
  }
  dsn <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                   weights = ~wt_pooled, data = dat_cc, nest = TRUE)
  fit <- svyglm(
    ir_binary ~ high_dehp + age + RIAGENDR + race + education +
                pir + bmi + waist + smoke,
    design = dsn, family = quasibinomial()
  )
  tmp <- tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "high_dehp")
  data.frame(
    set       = label,
    n         = nrow(dat_cc),
    estimate  = tmp$estimate,
    conf.low  = tmp$conf.low,
    conf.high = tmp$conf.high,
    p.value   = tmp$p.value,
    stringsAsFactors = FALSE
  )
}

main_full <- re_estimate_main(df,      "Full sample")
main_lag  <- re_estimate_main(df_lag,  "Lag-restricted (stable diet)")

cat(sprintf("  Full      : OR = %.3f (%.3f-%.3f), p=%.3g, n=%d\n",
            main_full$estimate, main_full$conf.low, main_full$conf.high,
            main_full$p.value, main_full$n))
cat(sprintf("  Lag-rest  : OR = %.3f (%.3f-%.3f), p=%.3g, n=%d\n",
            main_lag$estimate,  main_lag$conf.low,  main_lag$conf.high,
            main_lag$p.value,  main_lag$n))

lag_pass <- TRUE
ratio <- main_lag$estimate / main_full$estimate
if (abs(log(ratio)) > log(1.5)) {  # >50% shift
  cat("  WARNING: Lag-restricted OR shifts >50% relative to full sample\n")
  lag_pass <- FALSE
} else {
  cat("  -> Lag sensitivity PASS (effect estimate stable within 50%)\n")
}

# ============================================
# (2) Negative control outcomes
# ============================================
cat("\n[2] Negative control outcomes (Lipsitch 2010)...\n")

# Candidate NCO vars: biologically distant from phthalate -> IR axis
nco_candidates <- c(
  "LBXHGB"   = "Hemoglobin (g/dL)",
  "LBXRBCSI" = "Red blood cell count (10^6/uL)",
  "LBXMCVSI" = "Mean corpuscular volume (fL)",
  "LBXPLTSI" = "Platelet count (10^3/uL)"
)

nco_avail <- intersect(names(nco_candidates), names(df))
cat(sprintf("NCO candidates available: %s\n", paste(nco_avail, collapse = ", ")))

# Fit linear model: NCO ~ high_dehp + same covariates as main
nco_rows <- list()
for (v in nco_avail) {
  # Skip if too few obs
  if (sum(!is.na(df[[v]])) < 200) next

  needed_nco <- c(v, "high_dehp","age","RIAGENDR","race","education",
                  "pir","bmi","waist","smoke","SDMVPSU","SDMVSTRA","wt_pooled")
  dat_v <- df[complete.cases(df[, needed_nco]), ]
  for (fv in c("race","education","smoke")) {
    if (is.factor(dat_v[[fv]])) dat_v[[fv]] <- droplevels(dat_v[[fv]])
  }
  if (nrow(dat_v) < 200) next

  dsn_v <- tryCatch(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                              weights = ~wt_pooled, data = dat_v, nest = TRUE),
                    error = function(e) NULL)
  if (is.null(dsn_v)) next

  fmla <- as.formula(paste0(v, " ~ high_dehp + age + RIAGENDR + race + education +
                              pir + bmi + waist + smoke"))
  fit <- tryCatch(
    svyglm(fmla, design = dsn_v, family = gaussian()),
    error = function(e) NULL)
  if (is.null(fit)) next

  tt <- tryCatch(tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) next
  row <- tt %>% filter(term == "high_dehp")
  if (nrow(row) == 0) next

  pass_flag <- (row$p.value > 0.05)  # PASS = no significant association
  nco_rows[[length(nco_rows) + 1]] <- data.frame(
    NCO         = v,
    description = nco_candidates[v],
    beta        = round(row$estimate, 4),
    conf.low    = round(row$conf.low, 4),
    conf.high   = round(row$conf.high, 4),
    p.value     = signif(row$p.value, 3),
    n_used      = nrow(dat_v),
    NCE_pass    = pass_flag,
    stringsAsFactors = FALSE
  )
}

nco_df <- do.call(rbind, nco_rows)
print(nco_df)

n_pass <- sum(nco_df$NCE_pass)
n_total <- nrow(nco_df)
cat(sprintf("\nNegative control overall: %d / %d pass (p>0.05)\n", n_pass, n_total))

# ============================================
# Output
# ============================================
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

combined <- bind_rows(
  main_full %>% transmute(
    test = "Main full sample",
    n = n,
    estimate = round(estimate, 4),
    conf.low = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    p.value = signif(p.value, 3),
    pass = NA,
    note = "reference"),
  main_lag %>% transmute(
    test = "Lag-restricted (stable diet)",
    n = n,
    estimate = round(estimate, 4),
    conf.low = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    p.value = signif(p.value, 3),
    pass = lag_pass,
    note = sprintf("ratio_lag/full=%.2f", ratio)),
  nco_df %>% transmute(
    test = paste0("NCO ", NCO, " ", description),
    n = n_used,
    estimate = beta,
    conf.low = conf.low,
    conf.high = conf.high,
    p.value = p.value,
    pass = NCE_pass,
    note = "NCE: PASS = p>0.05")
)

write.csv(combined, "output/tables/lag_negative_control.csv", row.names = FALSE)
cat("Saved: output/tables/lag_negative_control.csv\n")

# Optional: small NCO forest plot
nco_plot <- nco_df %>%
  ggplot(aes(x = beta, y = description, color = NCE_pass)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d7263d"),
                     labels = c("TRUE" = "PASS (p>0.05)", "FALSE" = "FAIL (p<0.05)"),
                     name = "NCE result") +
  labs(title = "Negative control outcomes: high Sigma-DEHP",
       subtitle = sprintf("%d / %d NCOs pass (expected no association)",
                          n_pass, n_total),
       x = "Beta (NCO units)", y = NULL) +
  theme_minimal(base_size = 11)

ggsave("output/figures/lag_negative_control.png", nco_plot,
       width = 9, height = 4, dpi = 150)
cat("Saved: output/figures/lag_negative_control.png\n")

cat("\n========================================\n")
cat("27 Lag + Negative Control complete\n")
cat(sprintf("  -> Lag sensitivity: %s (full OR %.3f vs lag OR %.3f)\n",
            ifelse(lag_pass, "PASS", "FAIL"),
            main_full$estimate, main_lag$estimate))
cat(sprintf("  -> NCO: %d / %d pass\n", n_pass, n_total))
cat("========================================\n")
