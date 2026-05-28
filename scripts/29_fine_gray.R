# ============================================
# 009 / scripts/29_fine_gray.R
# Fine-Gray competing-risk model for cardiometabolic-specific mortality
#
# W16 Round 1 reset — R-Stats C2 fix:
# Manuscript §2.6 declares "Fine-Gray competing-risk models for cardiometabolic-
# specific mortality". Previously absent in code. Implement using
# cmprsk::crr (Fine-Gray subdistribution hazard).
#
# 输入: data/processed/nhanes_final.RData
# 输出: output/tables/fine_gray_cm_mortality.csv
#       output/tables/fine_gray_retrodesign_warning.csv (Type S/M for 64-68 events)
#
# Methodology (Fine-Gray 1999 JASA):
#   - mort_cm = 1: event of interest (cardiometabolic death, UCOD 1/5/7)
#   - mort_other = 1: competing event (death from other cause: cancer, accident, etc.)
#   - censored = neither
# Subdistribution HR (sHR) with 95% CI + Wald P, for Σ-DEHP-z + 8 metabolites
# ============================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(cmprsk)
})

cat("========================================\n")
cat("009 / 29_fine_gray.R — Fine-Gray CM mortality (R-Stats C2 fix)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

# ----------------------------------------------------------
# Build competing-risk outcome
# 0 = censored, 1 = CM death (event of interest), 2 = other death (competing)
# ----------------------------------------------------------
df_fg <- df_mort %>%
  dplyr::filter(!is.na(permth), permth > 0,
                !is.na(mort_allcause), !is.na(mort_cm))

# Status: 0/1/2
df_fg$fg_status <- with(df_fg, ifelse(
  mort_allcause == 0, 0L,
  ifelse(mort_cm == 1, 1L, 2L)
))

cat(sprintf("Fine-Gray cohort: N=%d\n", nrow(df_fg)))
cat(sprintf("  Censored (status=0):     n=%d\n", sum(df_fg$fg_status == 0)))
cat(sprintf("  CM deaths (status=1):    n=%d (event of interest)\n", sum(df_fg$fg_status == 1)))
cat(sprintf("  Other deaths (status=2): n=%d (competing)\n", sum(df_fg$fg_status == 2)))

# ----------------------------------------------------------
# Exposures
# ----------------------------------------------------------
phth_exposures <- c(
  "URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
  "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z",
  "sum_dehp_mol_z","sum_hmw_z","sum_lmw_z"
)
phth_exposures <- intersect(phth_exposures, names(df_fg))

# Cast factors for model matrix
for (v in c("race","education","smoke")) {
  if (v %in% names(df_fg) && !is.factor(df_fg[[v]])) {
    df_fg[[v]] <- factor(df_fg[[v]])
  }
}

# M2-style covariates (pre-exposure + adiposity, matching 07_cox structure)
cov_vars <- c("age","sex_male","race","education","pir","bmi","smoke","hypertension")
cov_vars <- intersect(cov_vars, names(df_fg))

# ----------------------------------------------------------
# Helper: run Fine-Gray for one exposure
# ----------------------------------------------------------
run_fine_gray <- function(df, exposure, cov_vars) {
  keep <- c("permth","fg_status",exposure,cov_vars)
  d <- df[, keep, drop = FALSE]
  d <- d[complete.cases(d), ]
  if (nrow(d) < 30) return(NULL)
  # Build design matrix
  fml_cov <- as.formula(paste("~", paste(c(exposure, cov_vars), collapse = " + ")))
  mm <- tryCatch(model.matrix(fml_cov, data = d)[, -1, drop = FALSE],
                 error = function(e) NULL)
  if (is.null(mm)) return(NULL)

  fit <- tryCatch(
    cmprsk::crr(ftime = d$permth, fstatus = d$fg_status, cov1 = mm,
                failcode = 1, cencode = 0),
    error = function(e) { cat(sprintf("  [err] %s: %s\n", exposure, conditionMessage(e))); NULL }
  )
  if (is.null(fit)) return(NULL)

  # Find the row corresponding to the exposure
  cf_names <- names(fit$coef)
  exp_idx <- which(cf_names == exposure)
  if (length(exp_idx) == 0) {
    # Try with mm-aware naming
    exp_idx <- which(grepl(paste0("^", exposure), cf_names))
  }
  if (length(exp_idx) == 0) return(NULL)
  exp_idx <- exp_idx[1]

  beta <- fit$coef[exp_idx]
  se   <- sqrt(diag(fit$var)[exp_idx])
  shr  <- exp(beta)
  lo   <- exp(beta - 1.96 * se)
  hi   <- exp(beta + 1.96 * se)
  z    <- beta / se
  p    <- 2 * pnorm(-abs(z))

  data.frame(
    exposure = exposure,
    sHR = shr,
    lo  = lo,
    hi  = hi,
    p   = p,
    n_obs = nrow(d),
    n_cm_event = sum(d$fg_status == 1),
    n_other_event = sum(d$fg_status == 2),
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# Run Fine-Gray for all exposures
# ----------------------------------------------------------
res_list <- list()
for (exp in phth_exposures) {
  cat(sprintf("  Fine-Gray: %s ...\n", exp))
  r <- run_fine_gray(df_fg, exp, cov_vars)
  if (!is.null(r)) res_list[[length(res_list)+1]] <- r
}

res_df <- if (length(res_list) > 0) do.call(rbind, res_list) else data.frame()

if (nrow(res_df) > 0) {
  res_df$sHR <- round(res_df$sHR, 4)
  res_df$lo  <- round(res_df$lo, 4)
  res_df$hi  <- round(res_df$hi, 4)
  res_df$p   <- signif(res_df$p, 4)
}

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(res_df, "output/tables/fine_gray_cm_mortality.csv", row.names = FALSE)
cat(sprintf("\n[OK] output/tables/fine_gray_cm_mortality.csv (rows=%d)\n", nrow(res_df)))

if (nrow(res_df) > 0) {
  cat("\n--- Fine-Gray CM mortality preview ---\n")
  print(res_df, row.names = FALSE)
}

# ----------------------------------------------------------
# Retrodesign warning (Type S/M for 64-68 CM events)
# Manuscript §3.8 acknowledges CM events count is in underpowered zone.
# This is a redundant standalone calculation; 32_retrodesign.R runs the full set.
# ----------------------------------------------------------
n_cm <- sum(df_fg$fg_status == 1)
n_total <- nrow(df_fg)
cat(sprintf("\n[retrodesign] CM events = %d / %d (%.1f%%); EPV (sample-based) = %.1f\n",
            n_cm, n_total, 100*n_cm/n_total, n_cm / length(cov_vars)))
cat(sprintf("[retrodesign] WARNING: events < 100 → Cox / Fine-Gray inference is in the\n"))
cat(sprintf("              underpowered zone (Peduzzi 1996 EPV rule; Gelman-Carlin 2014).\n"))
cat(sprintf("              See 32_retrodesign.R for formal Type S / Type M estimates.\n"))

write.csv(data.frame(
  metric = c("n_total","n_cm_event","n_other_event","epv_sample_based","cov_count"),
  value  = c(n_total, n_cm, sum(df_fg$fg_status==2), round(n_cm/length(cov_vars), 2), length(cov_vars))
), "output/tables/fine_gray_retrodesign_warning.csv", row.names = FALSE)

cat("\n========================================\n")
cat("Fine-Gray done.\n")
cat("========================================\n")
