# ============================================
# 009 / scripts/26_episensr.R
# Probabilistic Bias Analysis (Lash 2009; episensr R package 2.1.0+)
# Phthalate (high Sigma-DEHP Q4 vs Q1-Q3) -> IR binary
#
# 3 sensitivity targets:
#   (a) Phthalate measurement error  (CV 30-50%, non-differential)
#   (b) Outcome misclassification    (HOMA-IR sens 0.80 / spec 0.95)
#   (c) Unmeasured confounder (E-value-equivalent assumed RR 1.5-2.5)
#
# Monte Carlo iter: 100,000 (100k)
#
# Output:
#   output/tables/episensr_phth_ir.csv
#   output/figures/episensr_bias_adjusted.png
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(episensr); library(tidyr)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 26 Probabilistic Bias Analysis (episensr 100k iter)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

# ---- Build 2x2 table: high_dehp vs ir_binary ----
df <- nhanes_final
q_breaks <- quantile(df$sum_dehp_mol, probs = seq(0, 1, 0.25), na.rm = TRUE)
df$dehp_q <- cut(df$sum_dehp_mol, breaks = q_breaks, include.lowest = TRUE,
                 labels = c("Q1","Q2","Q3","Q4"))
df$high_dehp <- as.integer(df$dehp_q == "Q4")

# 2x2 table: rows = exposure (1=high, 0=low), cols = outcome (1=IR, 0=no IR)
tab <- with(df, table(factor(high_dehp, levels = c(1, 0)),
                      factor(ir_binary, levels = c(1, 0))))
print(tab)

a <- tab[1, 1]; b <- tab[1, 2]; c <- tab[2, 1]; d <- tab[2, 2]
obs_or <- (a * d) / (b * c)
obs_rr <- (a / (a + b)) / (c / (c + d))
cat(sprintf("\nObserved OR (high vs not-high DEHP -> IR): %.3f\n", obs_or))
cat(sprintf("Observed RR (high vs not-high DEHP -> IR): %.3f\n", obs_rr))

# Helper to pull OR -- total error row from episensr 2.1.0+ output
# adj_measures rows (probsens):
#   1: Relative Risk -- systematic error
#   2: Relative Risk -- total error
#   3: Odds Ratio   -- systematic error
#   4: Odds Ratio   -- total error
get_or_total <- function(res) {
  m <- res$adj_measures
  if (is.null(m) || nrow(m) < 4) return(c(NA, NA, NA))
  as.numeric(m[4, ])  # row 4 = OR total error
}
# For probsens_conf, OR rows differ; pick last OR row (OR SMR -- systematic + random)
get_or_total_conf <- function(res) {
  m <- res$adj_measures
  if (is.null(m) || nrow(m) < 4) return(c(NA, NA, NA))
  as.numeric(m[4, ])  # row 4 = OR SMR -- systematic + random
}

N_ITER <- 100000L

# ----------------------------------------------------------
# (a) Phthalate exposure measurement error (non-differential)
# ----------------------------------------------------------
cat("\n[a] Exposure measurement error (non-differential)...\n")
res_a <- tryCatch(probsens(
  case          = tab,
  type          = "exposure",
  reps          = N_ITER,
  seca          = list("trapezoidal", c(0.70, 0.80, 0.90, 0.95)),
  spca          = list("trapezoidal", c(0.80, 0.85, 0.95, 0.99))
), error = function(e) { cat("  episensr error:", conditionMessage(e), "\n"); NULL })

if (!is.null(res_a)) {
  adj_a <- get_or_total(res_a)
  cat(sprintf("  Adj OR (total error, median): %.3f (95%%SI %.3f, %.3f)\n",
              adj_a[1], adj_a[2], adj_a[3]))
}

# ----------------------------------------------------------
# (b) Outcome misclassification (HOMA-IR cutoff: sens 0.80 / spec 0.95)
# ----------------------------------------------------------
cat("\n[b] Outcome misclassification (HOMA-IR cutoff)...\n")
res_b <- tryCatch(probsens(
  case          = tab,
  type          = "outcome",
  reps          = N_ITER,
  seca          = list("trapezoidal", c(0.75, 0.80, 0.85, 0.90)),
  spca          = list("trapezoidal", c(0.90, 0.93, 0.97, 0.99))
), error = function(e) { cat("  episensr error:", conditionMessage(e), "\n"); NULL })

if (!is.null(res_b)) {
  adj_b <- get_or_total(res_b)
  cat(sprintf("  Adj OR (total error, median): %.3f (95%%SI %.3f, %.3f)\n",
              adj_b[1], adj_b[2], adj_b[3]))
}

# ----------------------------------------------------------
# (c) Unmeasured confounder
# ----------------------------------------------------------
cat("\n[c] Unmeasured confounder...\n")
res_c <- tryCatch(probsens_conf(
  case          = tab,
  reps          = N_ITER,
  prev_exp      = list("trapezoidal", c(0.35, 0.40, 0.45, 0.50)),
  prev_nexp     = list("trapezoidal", c(0.15, 0.20, 0.25, 0.30)),
  risk          = list("trapezoidal", c(1.3, 1.5, 2.0, 2.5))
), error = function(e) { cat("  episensr error:", conditionMessage(e), "\n"); NULL })

if (!is.null(res_c)) {
  adj_c <- get_or_total_conf(res_c)
  cat(sprintf("  Adj OR (total error, median): %.3f (95%%SI %.3f, %.3f)\n",
              adj_c[1], adj_c[2], adj_c[3]))
}

# ---- Save summary table ----
summary_rows <- list(
  c("Crude observed OR",     sprintf("%.3f", obs_or), NA, NA, "raw 2x2"),
  c("Crude observed RR",     sprintf("%.3f", obs_rr), NA, NA, "raw 2x2")
)

if (!is.null(res_a)) {
  v <- get_or_total(res_a)
  summary_rows[[length(summary_rows) + 1]] <- c(
    "Adj OR / exposure misclass (a)",
    sprintf("%.3f", v[1]),
    sprintf("%.3f", v[2]),
    sprintf("%.3f", v[3]),
    "Se [0.70-0.95] / Sp [0.80-0.99]"
  )
}

if (!is.null(res_b)) {
  v <- get_or_total(res_b)
  summary_rows[[length(summary_rows) + 1]] <- c(
    "Adj OR / outcome misclass (b)",
    sprintf("%.3f", v[1]),
    sprintf("%.3f", v[2]),
    sprintf("%.3f", v[3]),
    "Se [0.75-0.90] / Sp [0.90-0.99]"
  )
}

if (!is.null(res_c)) {
  v <- get_or_total_conf(res_c)
  summary_rows[[length(summary_rows) + 1]] <- c(
    "Adj OR / unmeasured conf (c)",
    sprintf("%.3f", v[1]),
    sprintf("%.3f", v[2]),
    sprintf("%.3f", v[3]),
    "RR_CD [1.3-2.5], prev_exp [0.35-0.50]"
  )
}

summary_df <- do.call(rbind, summary_rows) %>% as.data.frame(stringsAsFactors = FALSE)
names(summary_df) <- c("estimand", "median_OR", "ci_lo", "ci_hi", "params")

if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

write.csv(summary_df, "output/tables/episensr_phth_ir.csv", row.names = FALSE)
cat("\nSaved: output/tables/episensr_phth_ir.csv\n")

# ---- Plot bias-adjusted OR forest ----
plot_df <- summary_df %>%
  filter(!is.na(ci_lo)) %>%
  mutate(median_OR = as.numeric(median_OR),
         ci_lo     = as.numeric(ci_lo),
         ci_hi     = as.numeric(ci_hi))

p_forest <- plot_df %>%
  ggplot(aes(x = median_OR, y = estimand)) +
  geom_vline(xintercept = obs_or, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey80") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2,
                color = "#1b4965", orientation = "y") +
  geom_point(size = 3, color = "#d7263d") +
  scale_x_log10() +
  labs(title = "Probabilistic Bias Analysis: Sigma-DEHP Q4 vs Q1-Q3 -> IR",
       subtitle = sprintf("Crude OR = %.3f; 100k Monte Carlo iter", obs_or),
       x = "Bias-adjusted OR (log10)", y = NULL) +
  theme_minimal(base_size = 11)

ggsave("output/figures/episensr_bias_adjusted.png", p_forest,
       width = 9, height = 4, dpi = 150)
cat("Saved: output/figures/episensr_bias_adjusted.png\n")

# ---- Save underlying RData for inspection ----
save(res_a, res_b, res_c, obs_or, obs_rr, tab,
     file = "data/processed/episensr_009.RData")

cat("\n========================================\n")
cat("26 episensr Probabilistic Bias Analysis complete\n")
cat(sprintf("  -> Crude OR = %.3f\n", obs_or))
if (!is.null(res_a)) {
  v <- get_or_total(res_a); cat(sprintf("  -> Adj OR (a exposure ME)  = %.3f (95%%SI %.3f-%.3f)\n", v[1], v[2], v[3]))
}
if (!is.null(res_b)) {
  v <- get_or_total(res_b); cat(sprintf("  -> Adj OR (b outcome ME)   = %.3f (95%%SI %.3f-%.3f)\n", v[1], v[2], v[3]))
}
if (!is.null(res_c)) {
  v <- get_or_total_conf(res_c); cat(sprintf("  -> Adj OR (c uncontrolled) = %.3f (95%%SI %.3f-%.3f)\n", v[1], v[2], v[3]))
}
cat("========================================\n")
