# ============================================
# 009 / scripts/20_nri_idi.R
# Net Reclassification Index + Integrated Discrimination Improvement
# 输入: data/processed/tripod_ir_models.RData
# 输出: output/tables/nri_idi_phth.csv
#
# Baseline:   covariates only (a)
# Augmented:  covariates + Σ-DEHP / Σ-HMW / Σ-LMW
#
# Metrics:
#   - Continuous NRI (Pencina 2008)
#   - Categorical NRI (3-stratum: <0.30 / 0.30-0.60 / >0.60)
#   - IDI
#   - 1000 bootstrap CI
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(nricens); library(PredictABEL)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 20 NRI + IDI: covariates vs covariates + Phthalate Σ scores\n")
cat("========================================\n\n")

load("data/processed/tripod_ir_models.RData")
cat(sprintf("Cohort N = %d, IR=1: %d\n", nrow(df), sum(y==1)))

# Baseline model: covariates only
f_base <- as.formula(paste("ir_binary ~", paste(covars, collapse = " + ")))
# Augmented model: covariates + Σ-DEHP/HMW/LMW
aug_vars <- c(covars, "sum_dehp_mol_z", "sum_hmw_z", "sum_lmw_z")
f_aug  <- as.formula(paste("ir_binary ~", paste(aug_vars, collapse = " + ")))

fit_base <- suppressWarnings(glm(f_base, data = df, family = binomial()))
fit_aug  <- suppressWarnings(glm(f_aug,  data = df, family = binomial()))

p_base <- predict(fit_base, type = "response")
p_aug  <- predict(fit_aug,  type = "response")

cat(sprintf("\nBaseline 系数 n=%d, Aug 系数 n=%d\n",
            length(coef(fit_base)), length(coef(fit_aug))))

# ============================================
# Continuous NRI + IDI via nricens
# ============================================
# nricens::nribin needs event = 0/1 and the two probability vectors
nri_cont_obj <- tryCatch({
  nricens::nribin(
    event = df$ir_binary,
    p.std = p_base,
    p.new = p_aug,
    updown = "diff",         # continuous (any change)
    cut = 0,                 # continuous: any movement
    niter = 1000,
    msg = FALSE
  )
}, error = function(e) { cat("nribin diff err:", e$message, "\n"); NULL })

# Categorical NRI: cut at 0.30 / 0.60 (3 risk strata)
nri_cat_obj <- tryCatch({
  nricens::nribin(
    event = df$ir_binary,
    p.std = p_base,
    p.new = p_aug,
    updown = "category",
    cut = c(0.30, 0.60),
    niter = 1000,
    msg = FALSE
  )
}, error = function(e) { cat("nribin cat err:", e$message, "\n"); NULL })

# ============================================
# IDI via PredictABEL (also gives discrimination slope)
# ============================================
idi_obj <- tryCatch({
  PredictABEL::reclassification(
    data = data.frame(ir_binary = df$ir_binary),
    cOutcome = 1,
    predrisk1 = p_base,
    predrisk2 = p_aug,
    cutoff = c(0, 0.30, 0.60, 1.0)
  )
}, error = function(e) { cat("PredictABEL err:", e$message, "\n"); NULL })

# ============================================
# 手动 bootstrap CI for IDI (PredictABEL 不带 CI)
# IDI = (mean p_aug | y=1) - (mean p_aug | y=0)
#     - (mean p_base | y=1) - (mean p_base | y=0))
# ============================================
idi_stat <- function(pb, pa, yy) {
  (mean(pa[yy == 1]) - mean(pa[yy == 0])) -
    (mean(pb[yy == 1]) - mean(pb[yy == 0]))
}
idi_obs <- idi_stat(p_base, p_aug, df$ir_binary)

B <- 1000
boot_idi <- numeric(B)
boot_nri_cont <- numeric(B)
boot_nri_cat  <- numeric(B)
n <- nrow(df)
for (b in 1:B) {
  idx <- sample.int(n, replace = TRUE)
  d_b <- df[idx, ]
  fit_b1 <- suppressWarnings(glm(f_base, data = d_b, family = binomial()))
  fit_b2 <- suppressWarnings(glm(f_aug,  data = d_b, family = binomial()))
  pb_b <- predict(fit_b1, newdata = d_b, type = "response")
  pa_b <- predict(fit_b2, newdata = d_b, type = "response")
  yb <- d_b$ir_binary
  boot_idi[b] <- idi_stat(pb_b, pa_b, yb)

  # NRI continuous
  up_evt   <- mean((pa_b > pb_b)[yb == 1]) - mean((pa_b < pb_b)[yb == 1])
  up_noevt <- mean((pa_b < pb_b)[yb == 0]) - mean((pa_b > pb_b)[yb == 0])
  boot_nri_cont[b] <- up_evt + up_noevt

  # NRI categorical (3 strata)
  cat_b <- cut(pb_b, breaks = c(-Inf, 0.30, 0.60, Inf), labels = FALSE)
  cat_a <- cut(pa_b, breaks = c(-Inf, 0.30, 0.60, Inf), labels = FALSE)
  up_e <- mean((cat_a > cat_b)[yb == 1]) - mean((cat_a < cat_b)[yb == 1])
  up_n <- mean((cat_a < cat_b)[yb == 0]) - mean((cat_a > cat_b)[yb == 0])
  boot_nri_cat[b] <- up_e + up_n
}

# Point estimates (in-sample)
cat_obs1 <- cut(p_base, breaks = c(-Inf, 0.30, 0.60, Inf), labels = FALSE)
cat_obs2 <- cut(p_aug,  breaks = c(-Inf, 0.30, 0.60, Inf), labels = FALSE)
nri_cont_obs <- with(df, {
  ev <- mean((p_aug > p_base)[ir_binary == 1]) -
        mean((p_aug < p_base)[ir_binary == 1])
  ne <- mean((p_aug < p_base)[ir_binary == 0]) -
        mean((p_aug > p_base)[ir_binary == 0])
  ev + ne
})
nri_cat_obs <- with(df, {
  ev <- mean((cat_obs2 > cat_obs1)[ir_binary == 1]) -
        mean((cat_obs2 < cat_obs1)[ir_binary == 1])
  ne <- mean((cat_obs2 < cat_obs1)[ir_binary == 0]) -
        mean((cat_obs2 > cat_obs1)[ir_binary == 0])
  ev + ne
})

ci_q <- function(v) c(quantile(v, 0.025, na.rm=TRUE),
                      quantile(v, 0.975, na.rm=TRUE))

res <- data.frame(
  metric = c("Continuous NRI (any diff)",
             "Categorical NRI (cut 0.30 / 0.60)",
             "IDI"),
  estimate = c(nri_cont_obs, nri_cat_obs, idi_obs),
  ci_lo = c(ci_q(boot_nri_cont)[1], ci_q(boot_nri_cat)[1], ci_q(boot_idi)[1]),
  ci_hi = c(ci_q(boot_nri_cont)[2], ci_q(boot_nri_cat)[2], ci_q(boot_idi)[2]),
  bootstrap_iter = c(B, B, B)
)
# Two-sided p value: 2 * min(P(boot < 0), P(boot > 0))
two_sided_p <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(NA_real_)
  p_lo <- mean(v <= 0)
  p_hi <- mean(v >= 0)
  min(1, 2 * min(p_lo, p_hi))
}
res$p_value <- c(
  two_sided_p(boot_nri_cont),
  two_sided_p(boot_nri_cat),
  two_sided_p(boot_idi)
)

cat("\n--- NRI + IDI 结果 (manual bootstrap) ---\n")
print(round(res[, -1], 4))

# Save with extra context column
res$source <- "manual bootstrap (1000 iter, 2-sided P)"

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "output/tables/nri_idi_phth.csv", row.names = FALSE)

cat("\n输出: output/tables/nri_idi_phth.csv\n")
cat("========================================\n20 完成\n========================================\n")
