# ============================================
# 009 / scripts/09_rcs.R
# RCS (Hmisc::rcspline.eval 5-knot Harrell) Phthalate × HOMA-IR 剂量-反应
# 输入: data/processed/nhanes_final.RData + nhanes_design.RData
# 输出: output/tables/rcs_phth_homa.csv     (linear + non-linear P-value per metabolite × outcome)
#       output/figures/rcs_phth_homa.png   (8-panel grid)
#
# 8 metabolites × 2 outcomes (HOMA-IR continuous log-transformed, IR binary)
# Knots at 5/27.5/50/72.5/95 percentile (5-knot Harrell positions, Harrell 2015 RMS §2.4)
# W16 R-Stats C6 fix: switched from 4-knot rms-default (5/35/65/95) to
# 5-knot Harrell positions — manuscript §2.6 + Figure 3 caption declare 5-knot.
# 调整 M2: age + sex_male + race + education + pir + bmi + smoke + hypertension
# 用 svyglm with rcs() basis (复杂抽样兼容): 通过 design$variables 抽 rcs 列, 然后 svyglm 拟合
# ============================================

suppressPackageStartupMessages({
  library(survey); library(dplyr); library(rms); library(Hmisc); library(ggplot2); library(purrr)
})

cat("========================================\n")
cat("009 / RCS: Phthalate × HOMA-IR + IR binary\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

# Cast factors
for (v in c("race","education","smoke","sex_male","hypertension")) {
  if (v %in% names(nhanes_final) && !is.factor(nhanes_final[[v]])) {
    nhanes_final[[v]] <- factor(nhanes_final[[v]])
  }
}

design_main <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

# ----------------------------------------------------------
# Targets
# ----------------------------------------------------------
metabs <- c("URXMEP","URXMBP","URXMIB","URXMZP",
            "URXMHP","URXMHH","URXMOH","URXECP")
metabs <- intersect(metabs, names(nhanes_final))

# Use log2-transformed exposure (continuous, normal-shape, recommended for skewed pollutant)
phth_log_cols <- paste0(metabs, "_log2")
phth_log_cols <- intersect(phth_log_cols, names(nhanes_final))

cat(sprintf("RCS will run %d metabolites × 2 outcomes = %d models\n",
            length(phth_log_cols), length(phth_log_cols) * 2))

# Adjustment
adj <- "age + sex_male + race + education + pir + bmi + smoke + hypertension"

# W16 R-Stats C6 fix:
# 5-knot Harrell percentiles (5/27.5/50/72.5/95) — Harrell 2015 RMS §2.4
# Manuscript §2.6 + Figure 3 caption declare "5-knot Harrell positions".
# With n=2,239 and IR events 1,062, K=5 is well-supported.
knot_pct <- c(.05, .275, .50, .725, .95)

# ----------------------------------------------------------
# Helper: 跑一个 metabolite × outcome RCS
# 返回: linear P + non-linear P, plus prediction grid
# ----------------------------------------------------------
run_rcs <- function(design, log_var, outcome, fam = "binomial") {
  x <- design$variables[[log_var]]
  if (sum(!is.na(x)) < 100) return(NULL)
  knots <- quantile(x, knot_pct, na.rm = TRUE)
  # 5 knots: require ≥4 unique values for non-degenerate basis
  if (length(unique(knots)) < 4) return(NULL)

  # Build RCS basis manually using Hmisc::rcspline.eval (vector → matrix of K-1 cols)
  rcs_basis <- Hmisc::rcspline.eval(x, knots = knots, inclx = TRUE)
  # rcs_basis 第 1 列是原 x (inclx=TRUE), 后续是非线性 basis
  colnames(rcs_basis) <- c(paste0(log_var, "_rcs1"),
                          paste0(log_var, "_rcs", 2:ncol(rcs_basis)))
  # Attach to design
  for (k in seq_len(ncol(rcs_basis))) {
    design$variables[[colnames(rcs_basis)[k]]] <- rcs_basis[, k]
  }
  rcs_terms <- paste(colnames(rcs_basis), collapse = " + ")

  fml_str <- sprintf("%s ~ %s + %s", outcome, rcs_terms, adj)
  fam_obj <- if (fam == "binomial") quasibinomial() else gaussian()
  fit_full <- tryCatch(svyglm(as.formula(fml_str), design = design, family = fam_obj),
                       error = function(e) NULL)
  if (is.null(fit_full)) return(NULL)

  # Linear-only model (just rcs1) for likelihood ratio comparison
  fml_lin <- sprintf("%s ~ %s + %s", outcome, colnames(rcs_basis)[1], adj)
  fit_lin <- tryCatch(svyglm(as.formula(fml_lin), design = design, family = fam_obj),
                      error = function(e) NULL)
  # Null model (just adj)
  fml_null <- sprintf("%s ~ %s", outcome, adj)
  fit_null <- tryCatch(svyglm(as.formula(fml_null), design = design, family = fam_obj),
                       error = function(e) NULL)

  # Use survey::regTermTest for testing rcs terms
  # Non-linear: test rcs2... (=0)
  nonlin_terms <- colnames(rcs_basis)[-1]
  rt_nl <- tryCatch(regTermTest(fit_full, as.formula(paste("~", paste(nonlin_terms, collapse=" + ")))),
                    error = function(e) NULL)
  p_nonlin <- if (!is.null(rt_nl)) rt_nl$p else NA_real_

  # Linear (overall): test all rcs terms (=0)
  all_terms <- colnames(rcs_basis)
  rt_all <- tryCatch(regTermTest(fit_full, as.formula(paste("~", paste(all_terms, collapse=" + ")))),
                     error = function(e) NULL)
  p_overall <- if (!is.null(rt_all)) rt_all$p else NA_real_

  # Linear-only: p of rcs1 in linear model
  p_linear <- NA_real_
  if (!is.null(fit_lin)) {
    cf <- summary(fit_lin)$coefficients
    rcs1 <- colnames(rcs_basis)[1]
    if (rcs1 %in% rownames(cf)) p_linear <- cf[rcs1, "Pr(>|t|)"]
  }

  # Prediction grid for plotting
  xs <- seq(quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE),
            length.out = 100)
  rcs_pred <- Hmisc::rcspline.eval(xs, knots = knots, inclx = TRUE)
  colnames(rcs_pred) <- colnames(rcs_basis)
  # Get coefficients for rcs terms only
  coefs <- coef(fit_full)
  if (any(is.na(coefs))) return(NULL)
  # rcs basis x coef → log-odds or beta
  log_or <- as.numeric(rcs_pred %*% coefs[colnames(rcs_pred)])
  # Center at median x (reference)
  ref_idx <- which.min(abs(xs - median(x, na.rm = TRUE)))
  log_or_centered <- log_or - log_or[ref_idx]
  # SE via vcov - only rcs terms
  V <- vcov(fit_full)
  V_rcs <- V[colnames(rcs_pred), colnames(rcs_pred)]
  # Contrast vs reference: predict_i - predict_ref
  C <- sweep(rcs_pred, 2, rcs_pred[ref_idx, ], "-")
  se_contrast <- sqrt(rowSums((C %*% V_rcs) * C))

  pred_df <- data.frame(
    metabolite = sub("_log2$","", log_var),
    log_x = xs,
    est = log_or_centered,
    lo  = log_or_centered - 1.96 * se_contrast,
    hi  = log_or_centered + 1.96 * se_contrast,
    fam = fam,
    outcome = outcome,
    stringsAsFactors = FALSE
  )

  list(
    summary = data.frame(
      metabolite = sub("_log2$","", log_var),
      outcome = outcome,
      fam = fam,
      p_overall = p_overall,
      p_linear = p_linear,
      p_nonlinear = p_nonlin,
      n_used = sum(!is.na(x)),
      stringsAsFactors = FALSE),
    prediction = pred_df
  )
}

# ----------------------------------------------------------
# Run
# ----------------------------------------------------------
summaries <- list()
preds <- list()

# Outcome A: HOMA-IR continuous (homa_ir_log, gaussian)
for (lv in phth_log_cols) {
  cat(sprintf("  RCS: %s → homa_ir_log (gaussian)\n", lv))
  r <- run_rcs(design_main, lv, "homa_ir_log", fam = "gaussian")
  if (!is.null(r)) { summaries[[length(summaries)+1]] <- r$summary
                     preds[[length(preds)+1]] <- r$prediction }
}

# Outcome B: IR binary (ir_binary, quasibinomial)
for (lv in phth_log_cols) {
  cat(sprintf("  RCS: %s → ir_binary (binomial)\n", lv))
  r <- run_rcs(design_main, lv, "ir_binary", fam = "binomial")
  if (!is.null(r)) { summaries[[length(summaries)+1]] <- r$summary
                     preds[[length(preds)+1]] <- r$prediction }
}

summary_df <- do.call(rbind, summaries)
pred_df <- do.call(rbind, preds)

# Save summary
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(summary_df, "output/tables/rcs_phth_homa.csv", row.names = FALSE)
cat("\n[OK] output/tables/rcs_phth_homa.csv\n")

cat("\n--- summary preview ---\n"); print(summary_df, row.names = FALSE)

# ----------------------------------------------------------
# Figure: 8-panel grid
# ----------------------------------------------------------
# Plot log-HOMA outcome (gaussian) — clean continuous Y
pred_homa <- pred_df %>% filter(outcome == "homa_ir_log")

p_homa <- ggplot(pred_homa, aes(x = log_x, y = est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = .25) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ metabolite, scales = "free", ncol = 4) +
  labs(x = "log2 (urinary phthalate metabolite, ng/mL)",
       y = "Centered β for log HOMA-IR (95% CI)",
       title = "Restricted cubic spline: Phthalate metabolites and HOMA-IR (M2 adjusted)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
# W16 SA-3 fix: 200 -> 300 DPI (BMC final production floor)
ggsave("output/figures/rcs_phth_homa.png", p_homa,
       width = 12, height = 6, dpi = 300, bg = "white")
cat("[OK] output/figures/rcs_phth_homa.png (HOMA continuous, 300 DPI)\n")

# Plot ir_binary as well (log-odds)
pred_irb <- pred_df %>% filter(outcome == "ir_binary")
p_irb <- ggplot(pred_irb, aes(x = log_x, y = est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkorange", alpha = .25) +
  geom_line(color = "darkorange", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ metabolite, scales = "free", ncol = 4) +
  labs(x = "log2 (urinary phthalate metabolite, ng/mL)",
       y = "Centered log-odds for IR (HOMA≥2.5) (95% CI)",
       title = "RCS: Phthalate metabolites and IR binary (HOMA≥2.5, M2 adjusted)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# W16 SA-3 fix: 200 -> 300 DPI
ggsave("output/figures/rcs_phth_ir_binary.png", p_irb,
       width = 12, height = 6, dpi = 300, bg = "white")
cat("[OK] output/figures/rcs_phth_ir_binary.png (IR binary, 300 DPI)\n")

cat("\n========================================\n")
cat("RCS done.\n")
cat("========================================\n")
