# ============================================
# 009 / scripts/22_calibration.R
# Calibration belt (Nattino 2014 givitiR) + Hosmer-Lemeshow
# 输入: data/processed/tripod_ir_models.RData
# 输出: output/figures/calibration_belt.png + calibration_belt.pdf
#       output/tables/calibration_metrics.csv
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(givitiR); library(ResourceSelection)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 22 Calibration belt + HL test\n")
cat("========================================\n\n")

load("data/processed/tripod_ir_models.RData")

# Re-fit baseline + augmented + grab OOF probabilities from tripod_ir_models.RData (oof)
f_base <- as.formula(paste("ir_binary ~", paste(covars, collapse = " + ")))
aug_vars <- c(covars, "sum_dehp_mol_z", "sum_hmw_z", "sum_lmw_z")
f_aug  <- as.formula(paste("ir_binary ~", paste(aug_vars, collapse = " + ")))

fit_base <- suppressWarnings(glm(f_base, data = df, family = binomial()))
fit_aug  <- suppressWarnings(glm(f_aug,  data = df, family = binomial()))

df$p_base <- predict(fit_base, type = "response")
df$p_aug  <- predict(fit_aug,  type = "response")

# ============================================
# givitiR::givitiCalibrationBelt
# ============================================
make_belt <- function(p, y, name) {
  cb <- givitiR::givitiCalibrationBelt(o = y, e = p, devel = "internal")
  list(cb = cb, name = name)
}

cb_base <- make_belt(df$p_base, df$ir_binary, "Baseline (covariates only)")
cb_aug  <- make_belt(df$p_aug,  df$ir_binary, "Augmented (+ Phthalate)")

# Save calibration belt plot (combined panel)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# ============================================
# Hosmer-Lemeshow (g = 10) via ResourceSelection — compute FIRST so we can
# put hl_p in calibration belt title (W16 R-Figure P0 #2 fix)
# ============================================
hl_base <- ResourceSelection::hoslem.test(df$ir_binary, df$p_base, g = 10)
hl_aug  <- ResourceSelection::hoslem.test(df$ir_binary, df$p_aug,  g = 10)

# Two panel plot — title now shows Hosmer-Lemeshow P (matches manuscript)
# (Nattino 2014 belt_p is also reported but not as figure title — was misleading at 0.001)
png("output/figures/calibration_belt.png",
    width = 12, height = 5.5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(cb_base$cb, main = sprintf("%s\nHosmer-Lemeshow P = %.3f",
                                 cb_base$name, hl_base$p.value),
     xlab = "Predicted IR risk", ylab = "Observed IR")
plot(cb_aug$cb, main = sprintf("%s\nHosmer-Lemeshow P = %.3f",
                                cb_aug$name, hl_aug$p.value),
     xlab = "Predicted IR risk", ylab = "Observed IR")
dev.off()

pdf("output/figures/calibration_belt.pdf", width = 12, height = 5.5)
par(mfrow = c(1, 2))
plot(cb_base$cb, main = sprintf("%s\nHosmer-Lemeshow P = %.3f",
                                 cb_base$name, hl_base$p.value),
     xlab = "Predicted IR risk", ylab = "Observed IR")
plot(cb_aug$cb, main = sprintf("%s\nHosmer-Lemeshow P = %.3f",
                                cb_aug$name, hl_aug$p.value),
     xlab = "Predicted IR risk", ylab = "Observed IR")
dev.off()

# ============================================
# Calibration in the large + slope
# ============================================
cal_metrics <- function(p, y) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  fit <- suppressWarnings(glm(y ~ lp, family = binomial()))
  data.frame(
    cal_intercept = coef(fit)[1],
    cal_slope     = coef(fit)[2],
    brier         = mean((p - y)^2)
  )
}
cm_base <- cal_metrics(df$p_base, df$ir_binary)
cm_aug  <- cal_metrics(df$p_aug,  df$ir_binary)

metrics_df <- data.frame(
  model = c("Baseline (covariates only)", "Augmented (+ Phthalate)"),
  cal_intercept = c(cm_base$cal_intercept, cm_aug$cal_intercept),
  cal_slope     = c(cm_base$cal_slope,     cm_aug$cal_slope),
  brier         = c(cm_base$brier,         cm_aug$brier),
  hl_chisq      = c(as.numeric(hl_base$statistic), as.numeric(hl_aug$statistic)),
  hl_df         = c(hl_base$parameter, hl_aug$parameter),
  hl_p          = c(hl_base$p.value, hl_aug$p.value),
  belt_p        = c(cb_base$cb$p.value, cb_aug$cb$p.value)
)

cat("\n--- Calibration metrics ---\n")
print(round(metrics_df[, -1], 4))

write.csv(metrics_df, "output/tables/calibration_metrics.csv", row.names = FALSE)

cat("\n输出: output/figures/calibration_belt.png\n")
cat("输出: output/figures/calibration_belt.pdf\n")
cat("输出: output/tables/calibration_metrics.csv\n")
cat("========================================\n22 完成\n========================================\n")
