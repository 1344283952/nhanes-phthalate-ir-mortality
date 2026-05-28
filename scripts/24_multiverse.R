# ============================================
# 009 / scripts/24_multiverse.R
# Multiverse + Specification curve analysis (Simonsohn 2020 Nat Hum Behav; Steegen 2016)
# Phthalate × IR robustness
#
# Design (target ~144 specs):
#   exposure (3): sum_dehp_mol / sum_hmw / sum_lmw
#   cutoff (3): Q2/Q3/Q4 vs Q1 reference
#   outcome (2): ir_binary / homa_ir_log (continuous)
#   model spec (2): M1 minimal / M2 full
#   subset (4): Crude / fasting>=10h / age<60 / male only
#   = 3 * 3 * 2 * 2 * 4 = 144 specifications
#
# Output:
#   output/tables/multiverse_results.csv
#   output/figures/multiverse_curve.png
#   output/figures/multiverse_pvalue_dist.png
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(survey); library(broom)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 24 Multiverse + Specification curve (Phthalate x IR)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

# ---- Spec dimensions ----
exposures <- c("sum_dehp_mol", "sum_hmw", "sum_lmw")
exposure_labels <- c("sum_dehp_mol" = "Sigma-DEHP",
                     "sum_hmw"      = "Sigma-HMW",
                     "sum_lmw"      = "Sigma-LMW")
cutoff_levels <- c("Q2vsQ1", "Q3vsQ1", "Q4vsQ1")
outcomes  <- c("ir_binary", "homa_ir_log")
model_specs <- c("M1", "M2")
subsets   <- c("Crude", "Fasting>=10h", "Age<60", "Male")

# ---- Build quartile factors for each exposure ----
df <- nhanes_final
for (ex in exposures) {
  q_breaks <- quantile(df[[ex]], probs = seq(0, 1, 0.25), na.rm = TRUE)
  q_breaks <- unique(q_breaks)
  if (length(q_breaks) >= 5) {
    df[[paste0(ex, "_q")]] <- cut(df[[ex]], breaks = q_breaks,
                                   include.lowest = TRUE,
                                   labels = c("Q1","Q2","Q3","Q4"))
  } else {
    # fallback: rank-based
    df[[paste0(ex, "_q")]] <- factor(ntile(df[[ex]], 4),
                                      levels = 1:4,
                                      labels = c("Q1","Q2","Q3","Q4"))
  }
}

# ---- Covariate sets ----
# M1 minimal: age + sex + race
# M2 full: M1 + edu + pir + bmi + waist + smoke + drink + cycle
cov_m1 <- c("age", "RIAGENDR", "race")
cov_m2 <- c("age", "RIAGENDR", "race", "education", "pir",
            "bmi", "waist", "smoke", "cycle_tag")

# Filter subsets in advance to keep clean
df_all <- df
df_fasting10 <- df %>% filter(fasting_hours >= 10)
df_age60     <- df %>% filter(age < 60)
df_male      <- df %>% filter(RIAGENDR == 1)

subset_data <- list(
  "Crude"        = df_all,
  "Fasting>=10h" = df_fasting10,
  "Age<60"       = df_age60,
  "Male"         = df_male
)

cat("Subset N:\n")
for (s in names(subset_data)) cat(sprintf("  %-13s n=%d\n", s, nrow(subset_data[[s]])))

# ---- Run one specification ----
run_spec <- function(exposure, cutoff, outcome, model_spec, subset_name) {

  d <- subset_data[[subset_name]]
  if (nrow(d) < 100) return(NULL)

  qcol <- paste0(exposure, "_q")
  d$expq <- factor(d[[qcol]], levels = c("Q1","Q2","Q3","Q4"))

  if (sum(table(d$expq) > 5) < 2) return(NULL)

  cov_set <- if (model_spec == "M1") cov_m1 else cov_m2

  # Build formula
  rhs <- paste(c("expq", cov_set), collapse = " + ")
  fmla <- as.formula(paste0(outcome, " ~ ", rhs))

  # Survey design
  dsn <- tryCatch(
    svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
              weights = ~wt_pooled, data = d, nest = TRUE),
    error = function(e) NULL)
  if (is.null(dsn)) return(NULL)

  # Fit model
  fit <- tryCatch({
    if (outcome == "ir_binary") {
      svyglm(fmla, design = dsn, family = quasibinomial())
    } else {
      svyglm(fmla, design = dsn, family = gaussian())
    }
  }, error = function(e) NULL)

  if (is.null(fit)) return(NULL)

  # Extract effect for target cutoff (Q2/Q3/Q4 vs Q1)
  target_term <- paste0("expq", substr(cutoff, 1, 2))
  est <- tryCatch(tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(est)) return(NULL)
  row <- est %>% filter(term == target_term)
  if (nrow(row) == 0) return(NULL)

  exp_flag <- (outcome == "ir_binary")
  est_val <- if (exp_flag) exp(row$estimate) else row$estimate
  lo_val  <- if (exp_flag) exp(row$conf.low)  else row$conf.low
  hi_val  <- if (exp_flag) exp(row$conf.high) else row$conf.high

  data.frame(
    exposure   = exposure,
    cutoff     = cutoff,
    outcome    = outcome,
    model_spec = model_spec,
    subset     = subset_name,
    estimate   = est_val,
    conf.low   = lo_val,
    conf.high  = hi_val,
    p.value    = row$p.value,
    n          = nrow(d),
    stringsAsFactors = FALSE
  )
}

# ---- Run all specs ----
all_specs <- expand.grid(
  exposure   = exposures,
  cutoff     = cutoff_levels,
  outcome    = outcomes,
  model_spec = model_specs,
  subset     = subsets,
  stringsAsFactors = FALSE
)

cat(sprintf("\nTotal specifications to run: %d\n", nrow(all_specs)))

t0 <- Sys.time()
results <- pmap_dfr(
  all_specs,
  function(exposure, cutoff, outcome, model_spec, subset) {
    res <- run_spec(exposure, cutoff, outcome, model_spec, subset)
    if (is.null(res)) {
      data.frame(
        exposure   = exposure, cutoff = cutoff, outcome = outcome,
        model_spec = model_spec, subset = subset,
        estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        p.value  = NA_real_, n = NA_integer_,
        stringsAsFactors = FALSE)
    } else res
  })
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Completed %d specs in %.1fs\n", nrow(results), elapsed))

# ---- Summary stats ----
results_ok <- results %>% filter(!is.na(estimate))
n_total <- nrow(results)
n_ok    <- nrow(results_ok)

# For binary outcomes: OR>1 = positive (higher IR)
# For continuous outcomes: estimate>0 = positive (higher log HOMA-IR)
results_ok <- results_ok %>%
  mutate(direction = case_when(
    outcome == "ir_binary"   & estimate > 1 ~ "positive",
    outcome == "ir_binary"   & estimate < 1 ~ "negative",
    outcome == "homa_ir_log" & estimate > 0 ~ "positive",
    outcome == "homa_ir_log" & estimate < 0 ~ "negative",
    TRUE ~ "null"
  ),
  sig = (!is.na(p.value) & p.value < 0.05),
  effect_sign = sign(case_when(
    outcome == "ir_binary"   ~ estimate - 1,
    outcome == "homa_ir_log" ~ estimate,
    TRUE ~ 0
  ))) %>%
  arrange(effect_sign * estimate)

pct_pos <- mean(results_ok$direction == "positive") * 100
pct_neg <- mean(results_ok$direction == "negative") * 100
pct_sig_pos <- mean(results_ok$direction == "positive" & results_ok$sig) * 100
pct_sig_neg <- mean(results_ok$direction == "negative" & results_ok$sig) * 100

cat("\n=== Multiverse summary ===\n")
cat(sprintf("Total specs run     : %d\n", n_total))
cat(sprintf("Valid specs         : %d\n", n_ok))
cat(sprintf("Positive direction  : %.1f%%\n", pct_pos))
cat(sprintf("Negative direction  : %.1f%%\n", pct_neg))
cat(sprintf("Sig (p<0.05) pos    : %.1f%%\n", pct_sig_pos))
cat(sprintf("Sig (p<0.05) neg    : %.1f%%\n", pct_sig_neg))

# ---- Save results table ----
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# Add ordering column for plot
results_ok$spec_id <- seq_len(nrow(results_ok))

# Save table
results_ok %>%
  mutate(estimate_pretty = sprintf("%.3f (%.3f-%.3f)", estimate, conf.low, conf.high),
         p_pretty        = sprintf("%.3g", p.value)) %>%
  write.csv("output/tables/multiverse_results.csv", row.names = FALSE)

cat("\nSaved: output/tables/multiverse_results.csv\n")

# ---- Specification curve plot ----
# Two-panel layout: top = estimates with CIs ordered by effect size
#                   bottom = factor matrix indicating which factor levels are in each spec

# Top panel data
top_df <- results_ok %>%
  mutate(
    sort_metric = case_when(
      outcome == "ir_binary"   ~ log(pmax(estimate, 1e-6)),
      outcome == "homa_ir_log" ~ estimate,
      TRUE ~ 0
    )) %>%
  arrange(outcome, sort_metric) %>%
  group_by(outcome) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# For visualization, split by outcome
p_top_binary <- top_df %>% filter(outcome == "ir_binary") %>%
  ggplot(aes(x = rank, y = estimate, color = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, alpha = 0.4) +
  geom_point(size = 1) +
  scale_color_manual(values = c("TRUE" = "#d7263d", "FALSE" = "#1b4965"),
                     labels = c("TRUE" = "p<0.05", "FALSE" = "p>=0.05"),
                     name = "Significance") +
  scale_y_continuous(trans = "log10") +
  labs(title = "Specification curve: IR binary (OR scale, log)",
       x = "Specification (ordered by effect)",
       y = "Odds Ratio") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")

p_top_cont <- top_df %>% filter(outcome == "homa_ir_log") %>%
  ggplot(aes(x = rank, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, alpha = 0.4) +
  geom_point(size = 1) +
  scale_color_manual(values = c("TRUE" = "#d7263d", "FALSE" = "#1b4965"),
                     labels = c("TRUE" = "p<0.05", "FALSE" = "p>=0.05"),
                     name = "Significance") +
  labs(title = "Specification curve: log(HOMA-IR) continuous",
       x = "Specification (ordered by effect)",
       y = "Beta (log HOMA-IR scale)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")

# Save spec curve combined
png("output/figures/multiverse_curve.png", width = 1200, height = 900, res = 150)
gridExtra::grid.arrange(p_top_binary, p_top_cont, nrow = 2)
dev.off()

cat("Saved: output/figures/multiverse_curve.png\n")

# ---- P-value distribution histogram ----
p_pdist <- results_ok %>%
  ggplot(aes(x = p.value, fill = outcome)) +
  geom_histogram(binwidth = 0.025, boundary = 0, alpha = 0.7,
                 position = "identity", color = "white") +
  facet_wrap(~ outcome, scales = "free_y") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("ir_binary" = "#1b9e77", "homa_ir_log" = "#7570b3")) +
  labs(title = "P-value distribution across 144 specifications",
       subtitle = sprintf("Pct positive direction = %.1f%% / sig pos = %.1f%%",
                          pct_pos, pct_sig_pos),
       x = "p-value", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave("output/figures/multiverse_pvalue_dist.png", p_pdist,
       width = 9, height = 4, dpi = 150)
cat("Saved: output/figures/multiverse_pvalue_dist.png\n")

# ---- Save summary stats ----
summary_df <- data.frame(
  metric = c("total_specs", "valid_specs", "pct_positive_direction",
             "pct_negative_direction", "pct_sig_positive", "pct_sig_negative",
             "median_estimate_binary", "median_estimate_continuous"),
  value  = c(n_total, n_ok, round(pct_pos, 2),
             round(pct_neg, 2), round(pct_sig_pos, 2), round(pct_sig_neg, 2),
             round(median(results_ok$estimate[results_ok$outcome == "ir_binary"], na.rm = TRUE), 3),
             round(median(results_ok$estimate[results_ok$outcome == "homa_ir_log"], na.rm = TRUE), 3))
)
write.csv(summary_df, "output/tables/multiverse_summary.csv", row.names = FALSE)

cat("\n========================================\n")
cat("24 Multiverse + Spec curve 完成\n")
cat(sprintf("  -> %d specs / %d valid / %.1f%% positive\n",
            n_total, n_ok, pct_pos))
cat("========================================\n")
