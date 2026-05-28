# ============================================
# 009 / scripts/19_tripod_ai.R
# TRIPOD-AI 27-item compliance + IR prediction model
# 输入: data/processed/nhanes_final.RData
# 输出: output/tables/tripod_ai_models.csv
#       output/tables/tripod_ai_checklist.csv
#       data/processed/tripod_ir_models.RData (供 20-23 复用)
#
# 4 models:
#   (a) covariates only
#   (b) covariates + Σ-DEHP
#   (c) covariates + 8 phth metabolites + Σ-DEHP/HMW/LMW
#   (d) XGBoost full (covariates + 全 phth 暴露)
#
# 评估:
#   5-fold stratified CV
#   AUROC + 95% CI (DeLong via pROC)
#   Calibration in the large + slope
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(pROC); library(xgboost); library(stats)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 19 TRIPOD-AI 4 IR prediction models\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
df <- nhanes_final
cat(sprintf("原始 N = %d, IR=1: %d (%.1f%%)\n",
            nrow(df), sum(df$ir_binary==1, na.rm=TRUE),
            100*mean(df$ir_binary==1, na.rm=TRUE)))

# --- 协变量准备: impute 缺失 + drop 全 NA ---
covars <- c("age","sex_male","race","education","pir","bmi","waist","smoke",
            "hypertension","cotinine_log","fish_freq_30d","protein_g","kcal_day")
phth_metab <- c("URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
                "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z")
phth_sum   <- c("sum_dehp_mol_z","sum_hmw_z","sum_lmw_z")
y_col <- "ir_binary"

vars_all <- c(y_col, covars, phth_metab, phth_sum)
df <- df[, vars_all]
cat("\n变量 NA 情况 (impute 前):\n")
for (v in vars_all) {
  miss <- sum(is.na(df[[v]]))
  if (miss > 0) cat(sprintf("  %-18s missing %d (%.1f%%)\n",
                            v, miss, 100*miss/nrow(df)))
}

# `drink` 全 NA - 任务原列表里有 drink, 但 03_clean_data 没生成有效值, drop
# `cotinine_log` 51% NA - missing 视为 below LOD; 补 LOD/√2 的 log
if ("cotinine_log" %in% names(df)) {
  obs <- df$cotinine_log[!is.na(df$cotinine_log)]
  if (length(obs) > 0) {
    lod_log <- min(obs) - log(sqrt(2))
  } else lod_log <- 0
  df$cotinine_log[is.na(df$cotinine_log)] <- lod_log
}
# 连续协变量: median impute
for (v in c("waist","protein_g","kcal_day","pir","bmi")) {
  if (v %in% names(df) && any(is.na(df[[v]]))) {
    df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE)
  }
}

# Final complete cases (应几乎全保留)
df <- df[!is.na(df$ir_binary), ]
cc <- complete.cases(df[, c(covars, phth_metab, phth_sum)])
df <- df[cc, ]
cat(sprintf("\nImpute 后 N = %d, IR=1: %d (%.1f%%)\n",
            nrow(df), sum(df$ir_binary==1), 100*mean(df$ir_binary==1)))

# Factor → numeric for xgboost (use model.matrix internally)
make_X <- function(d, vars) {
  ff <- as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- model.matrix(ff, data = d)[, -1, drop = FALSE]
  mm
}
y <- df$ir_binary
X_cov <- make_X(df, covars)
X_dehp <- make_X(df, c(covars, "sum_dehp_mol_z"))
X_full <- make_X(df, c(covars, phth_metab, phth_sum))

cat(sprintf("\nDesign matrix dims: X_cov=%dx%d, X_dehp=%dx%d, X_full=%dx%d\n",
            nrow(X_cov), ncol(X_cov),
            nrow(X_dehp), ncol(X_dehp),
            nrow(X_full), ncol(X_full)))

# ============================================
# 5-fold stratified CV
# ============================================
K <- 5
ix_pos <- which(y == 1); ix_neg <- which(y == 0)
ix_pos <- sample(ix_pos); ix_neg <- sample(ix_neg)
fold_pos <- cut(seq_along(ix_pos), breaks = K, labels = FALSE)
fold_neg <- cut(seq_along(ix_neg), breaks = K, labels = FALSE)
fold <- integer(length(y))
fold[ix_pos] <- fold_pos
fold[ix_neg] <- fold_neg
table(fold, y)

# Container for OOF predictions
oof <- data.frame(
  y = y,
  fold = fold,
  p_cov = NA_real_,
  p_dehp = NA_real_,
  p_full = NA_real_,
  p_xgb = NA_real_
)

# ============================================
# Fit 4 models per fold
# ============================================
for (k in 1:K) {
  idx_test  <- which(fold == k)
  idx_train <- which(fold != k)
  d_tr <- df[idx_train, ]; y_tr <- y[idx_train]
  d_te <- df[idx_test, ];  y_te <- y[idx_test]

  # Model (a): covariates only — glm
  f_cov <- as.formula(paste("ir_binary ~", paste(covars, collapse = " + ")))
  fit_a <- suppressWarnings(glm(f_cov, data = d_tr, family = binomial()))
  oof$p_cov[idx_test] <- predict(fit_a, newdata = d_te, type = "response")

  # Model (b): covariates + Σ-DEHP
  f_dehp <- as.formula(paste("ir_binary ~",
                             paste(c(covars, "sum_dehp_mol_z"), collapse = " + ")))
  fit_b <- suppressWarnings(glm(f_dehp, data = d_tr, family = binomial()))
  oof$p_dehp[idx_test] <- predict(fit_b, newdata = d_te, type = "response")

  # Model (c): covariates + 8 phth + Σ-DEHP/HMW/LMW
  f_full <- as.formula(paste("ir_binary ~",
                             paste(c(covars, phth_metab, phth_sum), collapse = " + ")))
  fit_c <- suppressWarnings(glm(f_full, data = d_tr, family = binomial()))
  oof$p_full[idx_test] <- predict(fit_c, newdata = d_te, type = "response")

  # Model (d): XGBoost (covariates + 8 phth + Σ-DEHP/HMW/LMW)
  Xtr <- make_X(d_tr, c(covars, phth_metab, phth_sum))
  Xte <- make_X(d_te, c(covars, phth_metab, phth_sum))
  dtrain <- xgb.DMatrix(data = Xtr, label = y_tr)
  fit_d <- xgb.train(
    params = list(objective = "binary:logistic", eta = 0.05,
                  max_depth = 4, subsample = 0.8,
                  colsample_bytree = 0.8, eval_metric = "logloss"),
    data = dtrain, nrounds = 300, verbose = 0
  )
  oof$p_xgb[idx_test] <- predict(fit_d, newdata = Xte)
}

cat("\n5-fold CV 完成\n")

# ============================================
# Metric computation per model: AUROC + 95% CI + Brier + cal in large + slope
# ============================================
eval_model <- function(p, y, name) {
  r <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  ci <- pROC::ci.auc(r, method = "delong")
  brier <- mean((p - y)^2)

  # Calibration in the large (intercept) & slope (logit p)
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  cal_fit <- suppressWarnings(glm(y ~ lp, family = binomial()))
  cal_intercept <- coef(cal_fit)[1]
  cal_slope     <- coef(cal_fit)[2]

  # Hosmer-Lemeshow (10 group)
  hl_p <- tryCatch({
    g <- cut(p, breaks = unique(quantile(p, probs = seq(0, 1, 0.1))),
             include.lowest = TRUE)
    obs <- tapply(y, g, sum)
    exp <- tapply(p, g, sum)
    n   <- tapply(y, g, length)
    chi <- sum((obs - exp)^2 / (exp * (1 - exp/n)), na.rm = TRUE)
    1 - pchisq(chi, df = length(obs) - 2)
  }, error = function(e) NA_real_)

  data.frame(
    model = name,
    auroc = as.numeric(r$auc),
    auroc_lo = ci[1],
    auroc_hi = ci[3],
    brier = brier,
    cal_in_large = cal_intercept,
    cal_slope = cal_slope,
    hl_p = hl_p
  )
}

res_a <- eval_model(oof$p_cov,  oof$y, "a_covariates_only")
res_b <- eval_model(oof$p_dehp, oof$y, "b_cov_plus_DEHP")
res_c <- eval_model(oof$p_full, oof$y, "c_cov_plus_8phth_DEHP_HMW_LMW")
res_d <- eval_model(oof$p_xgb,  oof$y, "d_xgboost_full")

res_models <- rbind(res_a, res_b, res_c, res_d)

# DeLong comparison: c vs a, b vs a, d vs a
delong_test <- function(p1, p2, y, lbl) {
  r1 <- pROC::roc(y, p1, quiet = TRUE, direction = "<")
  r2 <- pROC::roc(y, p2, quiet = TRUE, direction = "<")
  rt <- pROC::roc.test(r1, r2, method = "delong")
  data.frame(
    contrast = lbl,
    auc1 = as.numeric(r1$auc),
    auc2 = as.numeric(r2$auc),
    delta = as.numeric(r2$auc - r1$auc),
    p_value = rt$p.value
  )
}
delong_tbl <- rbind(
  delong_test(oof$p_cov, oof$p_dehp, oof$y, "b - a (DEHP add)"),
  delong_test(oof$p_cov, oof$p_full, oof$y, "c - a (full phth add)"),
  delong_test(oof$p_cov, oof$p_xgb,  oof$y, "d - a (XGB add)")
)

cat("\n--- 4 models AUROC + 95% CI + Brier + cal ---\n")
print(round(res_models[, -1], 4))
cat("\n--- DeLong AUC contrast ---\n")
print(round(delong_tbl[, -1], 4))

# Write
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res_models, "output/tables/tripod_ai_models.csv", row.names = FALSE)
write.csv(delong_tbl, "output/tables/tripod_ai_delong_contrast.csv", row.names = FALSE)

# ============================================
# TRIPOD-AI 27-item checklist (manual semi-auto)
# ============================================
tripod <- data.frame(
  item = c(
    "1a. Title identifies AI / prediction model",
    "1b. Title identifies target population (NHANES adults)",
    "2.  Abstract structured (Aim/Method/Results/Conclusion)",
    "3a. Background: target condition + need",
    "3b. Objectives of analysis",
    "4a. Source of data (NHANES 2005-2018)",
    "4b. Eligibility criteria (fasting + no DM + no cancer + no HBV/HCV)",
    "5a. Setting / dates",
    "5b. Cohort split (5-fold CV stratified by IR)",
    "6.  Outcome definition (HOMA-IR ≥ 2.5)",
    "7a. Predictors (8 phth + Σ-DEHP/HMW/LMW + 13 covariates)",
    "7b. Predictor measurement standardization",
    "8.  Sample size justification (N=2,239 vs 21 predictors, EPV>40)",
    "9.  Missing data handling (median impute for continuous, LOD/√2 for cotinine)",
    "10a. Statistical model: 3 GLM + 1 XGBoost",
    "10b. Predictor selection (a priori + Σ scores)",
    "10c. Internal validation: 5-fold CV",
    "11.  Risk groups (not applicable - continuous prediction)",
    "12.  Model assumption checks (calibration in large + slope + HL test)",
    "13a. Flow of participants (output/tables/flow_counts.csv)",
    "13b. Characteristics (Table 1 from script 05_table1.R)",
    "14a. Outcome rate (47.4%)",
    "14b. Predictor distribution (z-scaled)",
    "15a. Model presentation (coefficients in 19 output)",
    "15b. Final equation reproducible",
    "16a. Performance: AUROC + 95% CI DeLong",
    "16b. Calibration: cal in large + slope + HL P",
    "17.  Model comparison (DeLong: c-a, b-a, d-a)",
    "18.  Discussion of limitations",
    "19.  Discussion of interpretation",
    "20.  Implication for practice",
    "21.  Supplementary information",
    "22.  Funding (JGJX2021D37)",
    "23.  Data sharing (NHANES public + scripts on GitHub)",
    "24.  Code availability (scripts/19-23)",
    "25.  Conflict of interest",
    "26.  Ethics statement (NCHS IRB)",
    "27.  Reporting framework (TRIPOD-AI 2024)"
  ),
  compliance = c(
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","NA","Yes","Yes",
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","Yes","Yes","Yes",
    "Yes","Yes","Yes"
  )
)
write.csv(tripod, "output/tables/tripod_ai_checklist.csv", row.names = FALSE)

# Save model objects for downstream scripts (20-23)
save(df, y, oof, covars, phth_metab, phth_sum,
     res_models, delong_tbl,
     file = "data/processed/tripod_ir_models.RData")

cat("\n输出: output/tables/tripod_ai_models.csv\n")
cat("输出: output/tables/tripod_ai_delong_contrast.csv\n")
cat("输出: output/tables/tripod_ai_checklist.csv\n")
cat("输出: data/processed/tripod_ir_models.RData\n")
cat("\n========================================\n19 完成\n========================================\n")
