# ============================================
# 009 / scripts/21_dca.R
# Decision Curve Analysis (Vickers 2006)
# 输入: data/processed/tripod_ir_models.RData
# 输出: output/figures/dca_phth_ir.png + dca_phth_ir.pdf
#       output/tables/dca_net_benefit.csv
#
# Strategies compared:
#   - Treat all
#   - Treat none
#   - Covariates only model
#   - Covariates + Phthalate (Σ-DEHP / Σ-HMW / Σ-LMW)
# Threshold range: 5% – 50%
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(dcurves)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 21 DCA: covariates vs covariates + Phthalate\n")
cat("========================================\n\n")

load("data/processed/tripod_ir_models.RData")
cat(sprintf("Cohort N = %d, IR=1: %d\n", nrow(df), sum(y==1)))

# Re-fit baseline + augmented (on full cohort for prediction)
f_base <- as.formula(paste("ir_binary ~", paste(covars, collapse = " + ")))
aug_vars <- c(covars, "sum_dehp_mol_z", "sum_hmw_z", "sum_lmw_z")
f_aug  <- as.formula(paste("ir_binary ~", paste(aug_vars, collapse = " + ")))

fit_base <- suppressWarnings(glm(f_base, data = df, family = binomial()))
fit_aug  <- suppressWarnings(glm(f_aug,  data = df, family = binomial()))

df$p_base <- predict(fit_base, type = "response")
df$p_aug  <- predict(fit_aug,  type = "response")

# ============================================
# dcurves::dca
# ============================================
dca_res <- dcurves::dca(
  formula = ir_binary ~ p_base + p_aug,
  data = df,
  thresholds = seq(0.05, 0.50, by = 0.01),
  label = list(p_base = "Covariates only",
               p_aug  = "Covariates + Phthalate")
)

# Net benefit table
nb_df <- as.data.frame(dca_res$dca)
cat("\n--- DCA 表头 (Net benefit) ---\n")
print(head(nb_df, 10))

# ============================================
# Plot
# ============================================
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

p_dca <- plot(dca_res, smooth = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Decision Curve Analysis: Insulin Resistance Prediction",
    subtitle = sprintf("NHANES 2005-2018  N = %d  IR cases = %d (%.1f%%)",
                       nrow(df), sum(df$ir_binary == 1),
                       100*mean(df$ir_binary == 1)),
    x = "Threshold probability",
    y = "Net benefit"
  ) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("output/figures/dca_phth_ir.png", p_dca,
       width = 7.5, height = 5.5, dpi = 300)
ggsave("output/figures/dca_phth_ir.pdf", p_dca,
       width = 7.5, height = 5.5)

# Save table
write.csv(nb_df, "output/tables/dca_net_benefit.csv", row.names = FALSE)

# Summary: max NB and incremental NB at 3 key thresholds
key_thr <- c(0.10, 0.20, 0.30, 0.40)
summary_rows <- list()
for (th in key_thr) {
  s_base <- nb_df[abs(nb_df$threshold - th) < 1e-6 & nb_df$label == "Covariates only", ]
  s_aug  <- nb_df[abs(nb_df$threshold - th) < 1e-6 & nb_df$label == "Covariates + Phthalate", ]
  if (nrow(s_base) > 0 && nrow(s_aug) > 0) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      threshold = th,
      nb_baseline = s_base$net_benefit[1],
      nb_augmented = s_aug$net_benefit[1],
      delta_nb = s_aug$net_benefit[1] - s_base$net_benefit[1]
    )
  }
}
summary_df <- do.call(rbind, summary_rows)
cat("\n--- ΔNB at key thresholds ---\n")
print(round(summary_df, 5))
write.csv(summary_df, "output/tables/dca_summary_thresholds.csv", row.names = FALSE)

cat("\n输出: output/figures/dca_phth_ir.png\n")
cat("输出: output/figures/dca_phth_ir.pdf\n")
cat("输出: output/tables/dca_net_benefit.csv\n")
cat("输出: output/tables/dca_summary_thresholds.csv\n")
cat("========================================\n21 完成\n========================================\n")
