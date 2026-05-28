# ============================================
# 009 / scripts/30_reri_ap_si.R
# RERI / AP / SI additive-scale interaction (VanderWeele-Knol 2014 Epidemiol Methods)
#
# W16 Round 1 reset — R-Causal CI-6 fix:
# Public-health attributable-fraction language requires ADDITIVE-scale interaction
# (RERI, AP, SI). Manuscript Discussion §4 ¶5 invokes "adiposity reduction reverses
# approximately half of phthalate-attributable IR component" — that claim is only
# legitimate on the additive scale.
#
# Targets: Σ-DEHP × adiposity (waist tertile high) → IR binary
# Sources:
#   VanderWeele 2014 Epidemiol Methods 3(1):33-72; Knol-VanderWeele 2012 Int J Epidemiol
#   RERI = RR11 - RR10 - RR01 + 1
#   AP   = RERI / RR11   (Attributable Proportion)
#   SI   = (RR11 - 1) / (RR10 + RR01 - 2)  (Synergy Index)
#
# Bootstrap CI (1000 reps) for RERI per VanderWeele 2014 recommendation
# ============================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 30_reri_ap_si.R — Additive-scale interaction (R-Causal CI-6 fix)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
options(survey.lonely.psu = "adjust")

# ----------------------------------------------------------
# Build binary exposure & binary effect modifier
# A = high Σ-DEHP  = top tertile (1) vs bottom 2 tertiles (0)
# M = high adiposity = top tertile waist (1) vs bottom 2 tertiles (0)
# ----------------------------------------------------------
d <- nhanes_final
d <- d[complete.cases(d[, c("sum_dehp_mol_z","waist","ir_binary","age","sex_male","race",
                            "education","pir","smoke","cotinine_log")]), ]

dehp_t <- quantile(d$sum_dehp_mol_z, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
d$A <- as.integer(d$sum_dehp_mol_z > dehp_t[3])

waist_t <- quantile(d$waist, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
d$M <- as.integer(d$waist > waist_t[3])

# Joint exposure category
d$AM <- with(d, paste0("A", A, "M", M))
n_by_cell <- table(d$AM, d$ir_binary)
cat("Cell counts (A {0,1} × M {0,1} × IR {0,1}):\n")
print(n_by_cell)

# Adjustment
cov_vars <- c("age","sex_male","race","education","pir","smoke","cotinine_log")

# ----------------------------------------------------------
# Fit weighted logistic with A × M interaction + cov
# ----------------------------------------------------------
des_d <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                   weights = ~wt_pooled, data = d, nest = TRUE)

fml_main <- as.formula(paste("ir_binary ~ A * M +", paste(cov_vars, collapse = " + ")))
fit_main <- survey::svyglm(fml_main, design = des_d, family = quasibinomial())

cf <- coef(fit_main)
vcv <- vcov(fit_main)

# Recover OR_A only (M=0), OR_M only (A=0), OR_AM joint
b_A  <- cf["A"]
b_M  <- cf["M"]
b_AM <- cf["A:M"]

OR10 <- exp(b_A)          # A=1, M=0 vs A=0, M=0
OR01 <- exp(b_M)          # A=0, M=1 vs A=0, M=0
OR11 <- exp(b_A + b_M + b_AM)  # A=1, M=1 vs A=0, M=0

# RERI = OR11 - OR10 - OR01 + 1
RERI <- OR11 - OR10 - OR01 + 1
AP   <- RERI / OR11
SI   <- (OR11 - 1) / (OR10 + OR01 - 2)

cat(sprintf("\nMain-table point estimates:\n"))
cat(sprintf("  OR10 (A=1,M=0) = %.4f\n", OR10))
cat(sprintf("  OR01 (A=0,M=1) = %.4f\n", OR01))
cat(sprintf("  OR11 (A=1,M=1) = %.4f\n", OR11))
cat(sprintf("  RERI = %.4f\n", RERI))
cat(sprintf("  AP   = %.4f\n", AP))
cat(sprintf("  SI   = %.4f\n", SI))

# ----------------------------------------------------------
# Bootstrap CI for RERI/AP/SI (500 reps, complete-case)
# ----------------------------------------------------------
cat("\nBootstrap 500 reps for RERI/AP/SI CI ...\n")

boot_one <- function(d_boot) {
  fit <- tryCatch(glm(fml_main, data = d_boot, family = quasibinomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, 3))
  cf <- coef(fit)
  if (!all(c("A","M","A:M") %in% names(cf))) return(rep(NA, 3))
  bA <- cf["A"]; bM <- cf["M"]; bAM <- cf["A:M"]
  or10 <- exp(bA); or01 <- exp(bM); or11 <- exp(bA + bM + bAM)
  reri <- or11 - or10 - or01 + 1
  ap   <- reri / or11
  si   <- (or11 - 1) / (or10 + or01 - 2)
  c(RERI = reri, AP = ap, SI = si)
}

n_boot <- 500L
boot_res <- matrix(NA_real_, nrow = n_boot, ncol = 3,
                   dimnames = list(NULL, c("RERI","AP","SI")))
for (i in seq_len(n_boot)) {
  idx <- sample.int(nrow(d), replace = TRUE)
  boot_res[i, ] <- tryCatch(boot_one(d[idx, ]),
                            error = function(e) rep(NA, 3))
}

# Percentile CI
ci_reri <- quantile(boot_res[, "RERI"], probs = c(.025, .975), na.rm = TRUE)
ci_ap   <- quantile(boot_res[, "AP"],   probs = c(.025, .975), na.rm = TRUE)
ci_si   <- quantile(boot_res[, "SI"],   probs = c(.025, .975), na.rm = TRUE)

# Approximate P (test against null = 0 for RERI, AP; null = 1 for SI)
p_reri <- 2 * min(mean(boot_res[, "RERI"] >= 0, na.rm = TRUE),
                  mean(boot_res[, "RERI"] <= 0, na.rm = TRUE))
p_ap   <- 2 * min(mean(boot_res[, "AP"]   >= 0, na.rm = TRUE),
                  mean(boot_res[, "AP"]   <= 0, na.rm = TRUE))
p_si   <- 2 * min(mean(boot_res[, "SI"]   >= 1, na.rm = TRUE),
                  mean(boot_res[, "SI"]   <= 1, na.rm = TRUE))

cat("\n--- RERI / AP / SI (VanderWeele-Knol 2014, additive scale) ---\n")
cat(sprintf("  RERI = %.4f (95%% boot CI %.4f, %.4f), P ≈ %.4f  (null = 0)\n",
            RERI, ci_reri[1], ci_reri[2], p_reri))
cat(sprintf("  AP   = %.4f (95%% boot CI %.4f, %.4f), P ≈ %.4f  (null = 0)\n",
            AP, ci_ap[1], ci_ap[2], p_ap))
cat(sprintf("  SI   = %.4f (95%% boot CI %.4f, %.4f), P ≈ %.4f  (null = 1)\n",
            SI, ci_si[1], ci_si[2], p_si))

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------
res_df <- data.frame(
  measure = c("RERI","AP","SI"),
  null    = c(0, 0, 1),
  point   = c(RERI, AP, SI),
  ci_lcl  = c(ci_reri[1], ci_ap[1], ci_si[1]),
  ci_ucl  = c(ci_reri[2], ci_ap[2], ci_si[2]),
  p_value = c(p_reri, p_ap, p_si),
  scale   = "additive (VanderWeele-Knol 2014)",
  interp  = c(
    "Relative excess risk due to interaction (>0 = positive additive)",
    "Proportion of OR11 attributable to interaction",
    "Synergy index (>1 = positive additive)"
  ),
  stringsAsFactors = FALSE
)
res_df$point  <- round(res_df$point, 4)
res_df$ci_lcl <- round(res_df$ci_lcl, 4)
res_df$ci_ucl <- round(res_df$ci_ucl, 4)
res_df$p_value <- signif(res_df$p_value, 4)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(res_df, "output/tables/reri_ap_si.csv", row.names = FALSE)

cat(sprintf("\n[OK] output/tables/reri_ap_si.csv (n_boot=%d, contrasts: A=Σ-DEHP top-tertile, M=waist top-tertile)\n",
            n_boot))
cat("========================================\n")
cat("RERI / AP / SI done.\n")
cat("========================================\n")
