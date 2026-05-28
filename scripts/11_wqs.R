# ============================================
# 009 / scripts/11_wqs.R
# WQS regression (Carrico 2015 J Agric Biol Environ Stat,
# DOI 10.1007/s13253-014-0180-3) for Phthalate mixture × IR (continuous + binary)
#
# 输入: data/processed/nhanes_final.RData
# 输出: output/tables/wqs_phth_homa.csv (HOMA-IR continuous)
#       output/tables/wqs_phth_ir_binary.csv (IR binary)
#       output/tables/wqs_phth_summary.csv (compact summary across all 4 fits)
#       output/tables/wqs_phth_weights.csv (component weights long-format)
#
# Mixture (8 Phthalate metabolite z-scores) — match 10_qgcomp.R
# Direction:
#   - Positive (b1_pos=TRUE): default, hypothesis = Phthalate → ↑ IR
#   - Negative (b1_pos=FALSE): sensitivity (some metabolites may protective via mole-weighted balance)
# Train/Valid split: 50/50 (validation=0.5)
# Bootstrap: b=100 on training
# Seed: 20260524
# Note: gWQS automatically standardizes via quantile (q=4) — feed raw z-scores OK
# ============================================

set.seed(20260524)

suppressPackageStartupMessages({
  library(dplyr)
  library(gWQS)
})

cat("========================================\n")
cat("009 / 11_wqs — Phthalate mixture × IR (continuous + binary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final loaded: n=%d ; IR cases=%d (%.1f%%)\n",
            nrow(nhanes_final),
            sum(nhanes_final$ir_binary == 1, na.rm = TRUE),
            100 * mean(nhanes_final$ir_binary == 1, na.rm = TRUE)))

# ---------------------------------------------------------------
# Step 1: Mixture + covariates
# ---------------------------------------------------------------
mixture <- c("URXMEP_z", "URXMBP_z", "URXMIB_z", "URXMZP_z",
             "URXMHP_z", "URXMHH_z", "URXMOH_z", "URXECP_z")

# Construct smoke_ever 0/1 if missing
if (!"smoke_ever" %in% names(nhanes_final)) {
  nhanes_final$smoke_ever <- ifelse(!is.na(nhanes_final$smoke) &
                                      nhanes_final$smoke == "Ever", 1L, 0L)
}

# Median impute cotinine_log for the half missing
if ("cotinine_log" %in% names(nhanes_final)) {
  median_cot <- median(nhanes_final$cotinine_log, na.rm = TRUE)
  nhanes_final$cotinine_log <- ifelse(is.na(nhanes_final$cotinine_log),
                                      median_cot, nhanes_final$cotinine_log)
}

cov_set <- c("age", "sex_male", "race", "education", "pir",
             "bmi", "waist", "smoke_ever", "cotinine_log")

core_keep <- c("SEQN", mixture, cov_set, "homa_ir_log", "ir_binary", "wt_pooled")
d <- nhanes_final %>%
  select(any_of(core_keep)) %>%
  filter(if_all(all_of(c(mixture, "age", "sex_male", "race", "education",
                          "pir", "bmi", "homa_ir_log", "ir_binary", "wt_pooled")),
                ~ !is.na(.)))
# Median-impute waist if a few NA; default smoke_ever 0
for (cv in c("waist")) {
  d[[cv]][is.na(d[[cv]])] <- median(d[[cv]], na.rm = TRUE)
}
d$smoke_ever[is.na(d$smoke_ever)] <- 0L

cat(sprintf("Analytic n: %d ; IR cases = %d\n", nrow(d), sum(d$ir_binary == 1)))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

# ---------------------------------------------------------------
# Step 2: run_wqs helper — supports gaussian/binomial × pos/neg
# ---------------------------------------------------------------
run_wqs <- function(outcome_var, family_str, direction) {
  is_pos <- (direction == "positive")
  cov_str <- paste(cov_set, collapse = " + ")
  fml <- as.formula(paste(outcome_var, "~ wqs +", cov_str))

  fit <- tryCatch(
    gWQS::gwqs(
      formula    = fml,
      mix_name   = mixture,
      data       = d,
      q          = 4,
      validation = 0.5,
      b          = 100,
      b1_pos     = is_pos,
      b1_constr  = TRUE,
      family     = family_str,
      seed       = 20260524,
      plan_strategy = "sequential"
    ),
    error = function(e) {
      cat("  WQS failed:", conditionMessage(e), "\n"); NULL
    })

  if (is.null(fit)) return(NULL)

  # Extract regression coefficient for the wqs index
  fs <- summary(fit$fit)
  co <- fs$coefficients
  if (!"wqs" %in% rownames(co)) {
    cat("  WARN: 'wqs' coef row missing\n"); return(NULL)
  }
  beta <- as.numeric(co["wqs", 1])
  se   <- as.numeric(co["wqs", 2])
  zval <- as.numeric(co["wqs", 3])
  pval <- as.numeric(co["wqs", 4])
  # 95% CI from beta ± 1.96 SE on linear-predictor scale
  ci_lcl <- beta - 1.96 * se
  ci_ucl <- beta + 1.96 * se

  # Component weights (gWQS returns mean of b bootstrap fits)
  fw <- fit$final_weights
  w_df <- data.frame(
    outcome   = outcome_var,
    direction = direction,
    metabolite = rownames(fw),
    weight    = fw$mean_weight,
    stringsAsFactors = FALSE)

  list(
    outcome = outcome_var, direction = direction,
    family = family_str,
    beta = beta, se = se, z = zval, p = pval,
    ci_lcl = ci_lcl, ci_ucl = ci_ucl,
    weights = w_df,
    n = nobs(fit$fit)
  )
}

# ---------------------------------------------------------------
# Step 3: run 4 fits: 2 outcomes × 2 directions
# ---------------------------------------------------------------
fits <- list(
  list(var = "homa_ir_log", fam = "gaussian", lbl = "HOMA-IR (log)"),
  list(var = "ir_binary",   fam = "binomial", lbl = "IR binary (HOMA>=2.5)")
)
directions <- c("positive", "negative")

all_results <- list()
all_weights <- list()
for (oc in fits) {
  cat(sprintf("\n--- %s [%s] ---\n", oc$lbl, oc$fam))
  for (dir in directions) {
    cat(sprintf("  direction = %s ...\n", dir))
    r <- run_wqs(oc$var, oc$fam, dir)
    if (!is.null(r)) {
      key <- paste0(oc$var, "_", dir)
      all_results[[key]] <- r
      all_weights[[key]] <- r$weights
      cat(sprintf("    beta=%.4f SE=%.4f p=%.4g  (n=%d)\n",
                  r$beta, r$se, r$p, r$n))
    }
  }
}

# ---------------------------------------------------------------
# Step 4: Per-outcome CSV (main = positive direction; neg = sensitivity)
# ---------------------------------------------------------------
make_outcome_csv <- function(outcome_var, scale_lbl, csv_path) {
  pos_key <- paste0(outcome_var, "_positive")
  neg_key <- paste0(outcome_var, "_negative")
  rows <- list()
  for (k in c(pos_key, neg_key)) {
    if (!is.null(all_results[[k]])) {
      r <- all_results[[k]]
      row <- data.frame(
        outcome = r$outcome, direction = r$direction,
        family  = r$family,  scale = scale_lbl,
        beta = r$beta, se = r$se,
        ci_lcl = r$ci_lcl, ci_ucl = r$ci_ucl,
        z = r$z, p = r$p, n = r$n,
        stringsAsFactors = FALSE
      )
      if (r$family == "binomial") {
        row$or <- exp(r$beta)
        row$or_lcl <- exp(r$ci_lcl)
        row$or_ucl <- exp(r$ci_ucl)
        row$effect_str <- sprintf("OR=%.3f (%.3f-%.3f), p=%.4g",
                                  row$or, row$or_lcl, row$or_ucl, row$p)
      } else {
        row$effect_str <- sprintf("beta=%.4f (%.4f to %.4f), p=%.4g",
                                  row$beta, row$ci_lcl, row$ci_ucl, row$p)
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  if (length(rows) > 0) {
    out_df <- dplyr::bind_rows(rows)
    write.csv(out_df, csv_path, row.names = FALSE)
    cat(sprintf("→ %s\n", csv_path))
    for (i in seq_len(nrow(out_df))) {
      cat(sprintf("  [%s] %s\n", out_df$direction[i], out_df$effect_str[i]))
    }
  }
}

make_outcome_csv("homa_ir_log", "log-scale", "output/tables/wqs_phth_homa.csv")
make_outcome_csv("ir_binary",   "OR",        "output/tables/wqs_phth_ir_binary.csv")

# ---------------------------------------------------------------
# Step 5: Combined summary + weights long-format
# ---------------------------------------------------------------
summary_rows <- lapply(all_results, function(r)
  data.frame(outcome = r$outcome, direction = r$direction, family = r$family,
             beta = r$beta, se = r$se, ci_lcl = r$ci_lcl, ci_ucl = r$ci_ucl,
             p = r$p, n = r$n, stringsAsFactors = FALSE))
summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else data.frame()
if (nrow(summary_df) > 0) {
  summary_df$p_BH <- p.adjust(summary_df$p, method = "BH")
  write.csv(summary_df, "output/tables/wqs_phth_summary.csv", row.names = FALSE)
  cat("\n→ output/tables/wqs_phth_summary.csv\n")
  print(summary_df)
}

weights_df <- if (length(all_weights) > 0) dplyr::bind_rows(all_weights) else data.frame()
if (nrow(weights_df) > 0) {
  write.csv(weights_df, "output/tables/wqs_phth_weights.csv", row.names = FALSE)
  cat("→ output/tables/wqs_phth_weights.csv\n")
}

# Save raw fits
save(all_results, summary_df, weights_df,
     file = "output/tables/wqs_phth_results.RData")

cat("\nDONE 11_wqs.R\n")
