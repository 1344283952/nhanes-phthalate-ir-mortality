# ============================================
# 009 / scripts/13_iptw.R
# Inverse Probability of Treatment Weighting (Austin 2011 Stat Med
# DOI 10.1002/sim.4067; Robins 2000 Epidemiology
# DOI 10.1097/00001648-200009000-00011)
# for High vs Low Phthalate (Q4 vs Q1-Q3 Σ-DEHP) → IR / mortality
#
# 输入: data/processed/nhanes_final.RData
# 输出: output/tables/iptw_phth_ir.csv  (ATE on IR binary + HOMA-IR continuous)
#       output/tables/iptw_phth_mort.csv (ATE on all-cause + CM mortality, Cox HR)
#       output/tables/iptw_phth_balance.csv (cobalt balance table SMD pre/post)
#       output/figures/iptw_love_plot.png (Love plot SMD pre/post)
#       output/tables/iptw_phth_results.RData (raw weight + fit bundle)
#
# Treatment definition:
#   high_dehp = 1 if sum_dehp_mol falls in Q4 of distribution (top quartile)
#   high_dehp = 0 if Q1-Q3 (reference)
# Covariates (predictors of treatment + confounders):
#   age, sex_male, race, education, pir, bmi, waist, smoke_ever,
#   hypertension, kcal_day, fish_freq_30d
#   Note: 'drink' completely missing in 03_clean → not used
#
# Method: WeightIt::weightit (method = "ps", logistic regression PS)
# Diagnostics:
#   - cobalt::bal.tab (SMD threshold 0.1)
#   - cobalt::love.plot (PNG saved)
# Outcome models (doubly robust):
#   - IR binary       → svyglm(family=quasibinomial) with IPTW
#   - HOMA-IR (log)   → svyglm(family=gaussian) with IPTW
#   - Mortality       → svycoxph(Surv(permth, event)) with IPTW
# Seed: 20260524
# ============================================

set.seed(20260524)

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
  library(survival)
  library(WeightIt)
  library(cobalt)
  library(ggplot2)
})
options(survey.lonely.psu = "adjust")

cat("========================================\n")
cat("009 / 13_iptw — IPTW Phthalate Q4 vs Q1-Q3 → IR + mortality\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final n=%d ; df_mort n=%d (all=%d / CM=%d)\n",
            nrow(nhanes_final), nrow(df_mort),
            sum(df_mort$mort_allcause, na.rm = TRUE),
            sum(df_mort$mort_cm, na.rm = TRUE)))

# ---------------------------------------------------------------
# Step 0: Construct treatment (high_dehp = Q4 of sum_dehp_mol)
# ---------------------------------------------------------------
q_breaks <- quantile(nhanes_final$sum_dehp_mol, probs = c(0, 0.25, 0.5, 0.75, 1.0),
                     na.rm = TRUE, names = FALSE)
cat(sprintf("Σ-DEHP-mol quartile cutoffs: Q1=%.4g, Q2=%.4g, Q3=%.4g, Q4=%.4g, max=%.4g\n",
            q_breaks[1], q_breaks[2], q_breaks[3], q_breaks[4], q_breaks[5]))

nhanes_final$high_dehp <- ifelse(!is.na(nhanes_final$sum_dehp_mol) &
                                   nhanes_final$sum_dehp_mol > q_breaks[4],
                                 1L, 0L)
df_mort$high_dehp <- nhanes_final$high_dehp[match(df_mort$SEQN, nhanes_final$SEQN)]

# smoke_ever 0/1 from factor
nhanes_final$smoke_ever <- ifelse(!is.na(nhanes_final$smoke) &
                                    nhanes_final$smoke == "Ever", 1L, 0L)
df_mort$smoke_ever      <- ifelse(!is.na(df_mort$smoke) &
                                    df_mort$smoke == "Ever", 1L, 0L)

# Median impute on covariates for stability
for (cv in c("waist","kcal_day","fish_freq_30d")) {
  if (cv %in% names(nhanes_final)) {
    med <- median(nhanes_final[[cv]], na.rm = TRUE)
    nhanes_final[[cv]][is.na(nhanes_final[[cv]])] <- med
  }
  if (cv %in% names(df_mort)) {
    medm <- median(df_mort[[cv]], na.rm = TRUE)
    df_mort[[cv]][is.na(df_mort[[cv]])] <- medm
  }
}

cat(sprintf("Treated (Q4 high_dehp): n=%d (%.1f%%) ; Control (Q1-Q3): n=%d\n",
            sum(nhanes_final$high_dehp == 1, na.rm = TRUE),
            100 * mean(nhanes_final$high_dehp == 1, na.rm = TRUE),
            sum(nhanes_final$high_dehp == 0, na.rm = TRUE)))

# ---------------------------------------------------------------
# Step 1: Treatment model covariates
# ---------------------------------------------------------------
ps_covs <- c("age", "sex_male", "race", "education", "pir",
             "bmi", "waist", "smoke_ever", "hypertension",
             "kcal_day", "fish_freq_30d")

ps_covs <- ps_covs[ps_covs %in% names(nhanes_final)]
cat(sprintf("PS covariates available: %s\n", paste(ps_covs, collapse = ", ")))

# Complete-case for IR analysis
ir_keep <- c("SEQN", "high_dehp", ps_covs, "ir_binary", "homa_ir", "homa_ir_log",
             "wt_pooled", "SDMVPSU", "SDMVSTRA")
df_ir <- nhanes_final %>%
  select(any_of(ir_keep)) %>%
  filter(if_all(all_of(c("high_dehp", ps_covs, "ir_binary", "homa_ir_log",
                          "wt_pooled")), ~ !is.na(.)))
cat(sprintf("IR analytic n=%d (Q4=%d, Q1-Q3=%d ; IR cases=%d)\n",
            nrow(df_ir), sum(df_ir$high_dehp == 1), sum(df_ir$high_dehp == 0),
            sum(df_ir$ir_binary == 1)))

# Mortality
mort_keep <- c("SEQN", "high_dehp", ps_covs, "mort_allcause", "mort_cm", "permth",
               "wt_pooled", "SDMVPSU", "SDMVSTRA")
df_m <- df_mort %>%
  select(any_of(mort_keep)) %>%
  filter(if_all(all_of(c("high_dehp", ps_covs, "mort_allcause", "permth",
                          "wt_pooled")), ~ !is.na(.)))
cat(sprintf("Mortality analytic n=%d (events=%d / CM=%d)\n",
            nrow(df_m), sum(df_m$mort_allcause), sum(df_m$mort_cm, na.rm = TRUE)))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ---------------------------------------------------------------
# Step 2: Estimate IPTW (ATE estimand) on IR cohort
# W16 R-Stats C5 fix: GBM propensity (method="gbm") + stabilize=TRUE + 99% trim
# Manuscript §2.5 + Figure 6 declare "generalised-boosted-model + stabilised + 99% cap".
# Falls back to logistic ("ps") if gbm package unavailable.
# ---------------------------------------------------------------
cat("\n[1/4] WeightIt::weightit (method='gbm', ATE, stabilize=TRUE) on IR cohort ...\n")
fml_ps <- as.formula(paste("high_dehp ~", paste(ps_covs, collapse = " + ")))

w_ir <- tryCatch({
  WeightIt::weightit(
    formula  = fml_ps,
    data     = df_ir,
    method   = "gbm",
    estimand = "ATE",
    stabilize = TRUE
  )
}, error = function(e) {
  cat("WARN gbm IR failed (", conditionMessage(e), ") → fallback to method='ps' logistic\n")
  tryCatch({
    WeightIt::weightit(
      formula  = fml_ps,
      data     = df_ir,
      method   = "ps",
      estimand = "ATE",
      stabilize = TRUE
    )
  }, error = function(e2) {
    cat("ERROR weightit IR fallback:", conditionMessage(e2), "\n"); NULL
  })
})

# W16 R-Stats C5 fix: 99% trim
if (!is.null(w_ir)) {
  w_ir <- tryCatch(
    WeightIt::trim(w_ir, at = .99),
    error = function(e) {
      cat("WARN trim 99% failed: ", conditionMessage(e), " — using untrimmed\n")
      w_ir
    }
  )
}

if (!is.null(w_ir)) {
  df_ir$iptw <- w_ir$weights
  # Combine with NHANES survey weight for doubly-robust population estimand
  df_ir$comb_w <- df_ir$iptw * df_ir$wt_pooled

  # ---- Balance diagnostics ----
  bt <- cobalt::bal.tab(w_ir, un = TRUE, m.threshold = 0.1)
  bal_df <- as.data.frame(bt$Balance)
  bal_df$variable <- rownames(bal_df)
  write.csv(bal_df, "output/tables/iptw_phth_balance.csv", row.names = FALSE)
  cat("→ output/tables/iptw_phth_balance.csv\n")

  # ---- Love plot ----
  lp <- tryCatch({
    cobalt::love.plot(w_ir, threshold = 0.1, abs = TRUE,
                      var.order = "unadjusted",
                      title = "Σ-DEHP Q4 vs Q1-Q3 IPTW balance (SMD)")
  }, error = function(e) {
    cat("WARN love.plot failed:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(lp)) {
    ggsave("output/figures/iptw_love_plot.png", plot = lp,
           width = 7, height = 5, dpi = 200)
    cat("→ output/figures/iptw_love_plot.png\n")
  }
}

# ---------------------------------------------------------------
# Step 3: Outcome models (doubly robust IPTW + survey weight)
# ---------------------------------------------------------------
ir_results <- list()
if (!is.null(w_ir)) {
  cat("\n[2/4] svyglm IPTW → IR binary (logistic) ...\n")
  des_ir <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                       weights = ~comb_w, data = df_ir, nest = TRUE)

  m_ir_bin <- tryCatch({
    survey::svyglm(ir_binary ~ high_dehp, design = des_ir,
                   family = quasibinomial())
  }, error = function(e) {
    cat("ERROR svyglm IR binary:", conditionMessage(e), "\n"); NULL
  })

  if (!is.null(m_ir_bin)) {
    co <- summary(m_ir_bin)$coefficients
    or <- exp(co["high_dehp", 1])
    lcl <- exp(co["high_dehp", 1] - 1.96 * co["high_dehp", 2])
    ucl <- exp(co["high_dehp", 1] + 1.96 * co["high_dehp", 2])
    ir_results[["ir_binary"]] <- data.frame(
      outcome = "IR binary (HOMA>=2.5)",
      treatment = "Σ-DEHP Q4 vs Q1-Q3",
      n_treated = sum(df_ir$high_dehp == 1),
      n_control = sum(df_ir$high_dehp == 0),
      scale = "OR",
      estimate = or, lcl = lcl, ucl = ucl,
      p = co["high_dehp", 4],
      effect_str = sprintf("OR=%.3f (%.3f-%.3f), p=%.4g", or, lcl, ucl, co["high_dehp", 4]),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  IR binary: %s\n", ir_results[["ir_binary"]]$effect_str))
  }

  cat("\n[3/4] svyglm IPTW → HOMA-IR (log, continuous) ...\n")
  m_ir_cont <- tryCatch({
    survey::svyglm(homa_ir_log ~ high_dehp, design = des_ir,
                   family = gaussian())
  }, error = function(e) {
    cat("ERROR svyglm HOMA cont:", conditionMessage(e), "\n"); NULL
  })

  if (!is.null(m_ir_cont)) {
    co <- summary(m_ir_cont)$coefficients
    beta <- co["high_dehp", 1]
    se   <- co["high_dehp", 2]
    lcl <- beta - 1.96 * se
    ucl <- beta + 1.96 * se
    ir_results[["homa_cont"]] <- data.frame(
      outcome = "HOMA-IR (log)",
      treatment = "Σ-DEHP Q4 vs Q1-Q3",
      n_treated = sum(df_ir$high_dehp == 1),
      n_control = sum(df_ir$high_dehp == 0),
      scale = "log-scale beta",
      estimate = beta, lcl = lcl, ucl = ucl,
      p = co["high_dehp", 4],
      effect_str = sprintf("beta=%.4f (%.4f to %.4f), p=%.4g",
                            beta, lcl, ucl, co["high_dehp", 4]),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  HOMA-IR log: %s\n", ir_results[["homa_cont"]]$effect_str))
  }
}

ir_df <- if (length(ir_results) > 0) dplyr::bind_rows(ir_results) else data.frame()
if (nrow(ir_df) > 0) {
  write.csv(ir_df, "output/tables/iptw_phth_ir.csv", row.names = FALSE)
  cat("\n→ output/tables/iptw_phth_ir.csv\n")
}

# ---------------------------------------------------------------
# Step 4: IPTW on mortality cohort + svycoxph
# ---------------------------------------------------------------
mort_results <- list()
# W16 R-Stats C5 fix: GBM + stabilize=TRUE + 99% trim for mortality cohort too
w_m <- tryCatch({
  WeightIt::weightit(
    formula  = fml_ps,
    data     = df_m,
    method   = "gbm",
    estimand = "ATE",
    stabilize = TRUE
  )
}, error = function(e) {
  cat("WARN gbm mort failed → fallback to method='ps' logistic stabilize=TRUE\n")
  tryCatch({
    WeightIt::weightit(
      formula  = fml_ps,
      data     = df_m,
      method   = "ps",
      estimand = "ATE",
      stabilize = TRUE
    )
  }, error = function(e2) {
    cat("ERROR weightit mort fallback:", conditionMessage(e2), "\n"); NULL
  })
})
if (!is.null(w_m)) {
  w_m <- tryCatch(
    WeightIt::trim(w_m, at = .99),
    error = function(e) {
      cat("WARN trim mort 99% failed — using untrimmed\n")
      w_m
    }
  )
}

if (!is.null(w_m)) {
  df_m$iptw <- w_m$weights
  df_m$comb_w <- df_m$iptw * df_m$wt_pooled
  des_m <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                      weights = ~comb_w, data = df_m, nest = TRUE)

  cat("\n[4/4] svycoxph IPTW → mortality (all-cause + CM) ...\n")
  m_all <- tryCatch({
    survey::svycoxph(Surv(permth, mort_allcause) ~ high_dehp, design = des_m)
  }, error = function(e) {
    cat("ERROR svycoxph all-cause:", conditionMessage(e), "\n"); NULL
  })
  m_cm <- tryCatch({
    survey::svycoxph(Surv(permth, mort_cm) ~ high_dehp, design = des_m)
  }, error = function(e) {
    cat("ERROR svycoxph CM:", conditionMessage(e), "\n"); NULL
  })

  add_cox_row <- function(m, lbl, events) {
    if (is.null(m)) return(NULL)
    co <- summary(m)$coefficients
    if (!"high_dehp" %in% rownames(co)) return(NULL)
    hr <- exp(co["high_dehp", 1])
    se <- co["high_dehp", 2]
    lcl <- exp(co["high_dehp", 1] - 1.96 * se)
    ucl <- exp(co["high_dehp", 1] + 1.96 * se)
    data.frame(
      outcome = lbl,
      treatment = "Σ-DEHP Q4 vs Q1-Q3",
      events = events,
      n_treated = sum(df_m$high_dehp == 1),
      n_control = sum(df_m$high_dehp == 0),
      scale = "HR",
      estimate = hr, lcl = lcl, ucl = ucl,
      p = co["high_dehp", ncol(co)],
      effect_str = sprintf("HR=%.3f (%.3f-%.3f), p=%.4g",
                           hr, lcl, ucl, co["high_dehp", ncol(co)]),
      stringsAsFactors = FALSE
    )
  }
  mort_results[["all"]] <- add_cox_row(m_all, "All-cause mortality",
                                         sum(df_m$mort_allcause))
  mort_results[["cm"]]  <- add_cox_row(m_cm,  "Cardiometabolic mortality",
                                         sum(df_m$mort_cm, na.rm = TRUE))
}

mort_df <- if (length(mort_results) > 0) dplyr::bind_rows(mort_results) else data.frame()
if (nrow(mort_df) > 0) {
  write.csv(mort_df, "output/tables/iptw_phth_mort.csv", row.names = FALSE)
  cat("\n→ output/tables/iptw_phth_mort.csv\n")
  for (i in seq_len(nrow(mort_df))) {
    cat(sprintf("  %s: %s\n", mort_df$outcome[i], mort_df$effect_str[i]))
  }
}

# ---------------------------------------------------------------
# Save bundle
# ---------------------------------------------------------------
save(w_ir, w_m, ir_df, mort_df,
     file = "output/tables/iptw_phth_results.RData")

cat("\n保存:\n")
cat("  output/tables/iptw_phth_ir.csv\n")
cat("  output/tables/iptw_phth_mort.csv\n")
cat("  output/tables/iptw_phth_balance.csv\n")
cat("  output/figures/iptw_love_plot.png\n")
cat("  output/tables/iptw_phth_results.RData\n")
cat("\nDONE 13_iptw.R\n")
