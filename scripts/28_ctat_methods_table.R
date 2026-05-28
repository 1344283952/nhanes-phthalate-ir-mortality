# ============================================
# 009 / scripts/28_ctat_methods_table.R
# CTAT 4-quadrant methods table (Hernan 2020 Causal Inference book;
# Diabetes Care / J Hep supplementary standard)
#
# Four quadrants per method:
#   (a) Estimand definition
#   (b) Identification assumption
#   (c) Estimation method
#   (d) Sensitivity analysis
#
# Methods covered:
#   1. BKMR (Bayesian kernel machine regression)
#   2. qgcomp (quantile g-computation)
#   3. WQS (weighted quantile sum)
#   4. CMAverse (4-way mediation decomposition)
#   5. IPTW (inverse probability of treatment weighting)
#   6. Bayesian g-computation
#   7. Negative control outcomes
#
# Output: output/tables/ctat_methods_table.csv (supplementary table S1)
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

cat("========================================\n")
cat("009 / 28 CTAT 4-quadrant methods table\n")
cat("========================================\n\n")

# ---- Construct rows ----
ctat_rows <- list(

  list(method = "BKMR (Bayesian kernel machine regression)",
       reference = "Bobb 2015 Biostatistics",
       estimand = paste0(
         "Joint exposure-response surface E[log(HOMA-IR) | Z = z*, X = x] ",
         "across all phthalate metabolites simultaneously; ",
         "estimands include single-exposure effect, posterior inclusion probability (PIP), ",
         "bivariate interaction, and overall mixture effect at population quantiles."),
       identification = paste0(
         "No unmeasured confounders given X; no exposure-response model misspecification (relaxed by ",
         "Gaussian kernel); positivity across observed mixture support."),
       estimation = paste0(
         "Markov chain Monte Carlo (2 chains x 10,000 iter; component-wise variable selection; ",
         "knots K=100; checkpointed)."),
       sensitivity = paste0(
         "Knots K=50 sensitivity; alternative h-prior; rolling-3 chain history; ",
         "compare with qgcomp / WQS triangulation.")
  ),

  list(method = "qgcomp (quantile g-computation)",
       reference = "Keil 2020 EHP",
       estimand = paste0(
         "Joint effect of 1-unit increase in all phthalate quantiles simultaneously ",
         "on log(HOMA-IR) (psi); plus signed component-weights."),
       identification = paste0(
         "Same as BKMR (no unmeasured confounders, additive linear form within quantile blocks). ",
         "Does NOT require positive monotonicity, but assumes linear within quantile."),
       estimation = paste0(
         "Quantize each phthalate to deciles; fit linear regression on summed quantile ",
         "with bootstrap (B=2000) for inference."),
       sensitivity = paste0(
         "Number of quantiles (4 vs 10); signed coefficient direction match; ",
         "ground in null permutation.")
  ),

  list(method = "WQS (weighted quantile sum)",
       reference = "Carrico 2015 J Agric Biol Environ Stat",
       estimand = paste0(
         "Directional mixture index w * z (w summing to 1) and effect on log(HOMA-IR); ",
         "yields a normalised weight per phthalate."),
       identification = paste0(
         "Directional homogeneity (all phthalates point same way); ",
         "same confounding assumption as qgcomp."),
       estimation = paste0(
         "40/60 train/validation split; B=200 bootstrap; positive- and negative-directional WQS."),
       sensitivity = paste0(
         "Random hold-out repetition; positive vs negative direction comparison; ",
         "compare top-weighted compound with BKMR PIP.")
  ),

  list(method = "CMAverse (4-way mediation decomposition)",
       reference = "Shi B 2021 Bioinformatics; VanderWeele 2014 Epi Methods",
       estimand = paste0(
         "Total effect decomposition: ",
         "(i) controlled direct effect CDE, ",
         "(ii) reference interaction INT_ref, ",
         "(iii) mediated interaction INT_med, ",
         "(iv) pure indirect effect PIE; ",
         "sequential mediators HSCRP -> adiposity -> cardiometabolic profile."),
       identification = paste0(
         "VanderWeele 4 assumptions: no confounding for (1) E-Y, (2) M-Y, (3) E-M, ",
         "(4) no exposure-induced M-Y confounding (relaxed for sequential)."),
       estimation = paste0(
         "Parametric outcome and mediator regressions with B=500 bootstrap. ",
         "Exposure: high_dehp Q4 vs Q1-Q3."),
       sensitivity = paste0(
         "E-value for each pathway; alternative mediator order; ",
         "include/exclude HOMA-IR proxy.")
  ),

  list(method = "IPTW (inverse probability of treatment weighting)",
       reference = "Robins 2000 Epidemiology",
       estimand = paste0(
         "Average treatment effect on the treated (ATT) / overall ATE of high_dehp on IR ",
         "after standardising covariate distribution between exposed and unexposed."),
       identification = paste0(
         "No unmeasured confounders; positivity (probability of high_dehp 0<p<1 for each X); ",
         "correct propensity model specification."),
       estimation = paste0(
         "Propensity score from logistic regression on M2 covariates; ",
         "stabilised IPW; Cobalt SMD < 0.1 balance check; ",
         "weighted outcome model svyglm."),
       sensitivity = paste0(
         "Trim extreme weights (1st/99th centile); SMD threshold sensitivity (0.05 vs 0.1); ",
         "alternative propensity model (Random Forest).")
  ),

  list(method = "Bayesian g-computation",
       reference = "Keil 2017 EpiMethods",
       estimand = paste0(
         "Counterfactual mean outcome E[Y^{do(high_dehp=1)}] - E[Y^{do(high_dehp=0)}] ",
         "averaged across observed covariate distribution (ATE on risk-difference and RR scale)."),
       identification = paste0(
         "Standard g-formula assumptions: no unmeasured confounders, positivity, ",
         "consistency, correct outcome model."),
       estimation = paste0(
         "Bayesian logistic outcome model (rstanarm, weakly informative N(0,2.5) priors, ",
         "HMC 2 chains x 4000 iter, warmup 1000); ",
         "posterior predictive at high_dehp = 0 and 1 for each subject; ",
         "mean over subjects per posterior draw."),
       sensitivity = paste0(
         "ROPE [-0.01, +0.01] for ATE; alternative prior scale N(0,5); ",
         "comparison with frequentist g-comp and IPTW.")
  ),

  list(method = "Negative control outcomes",
       reference = "Lipsitch 2010 Epidemiology",
       estimand = paste0(
         "Conditional E[NCO | high_dehp, X] - should be null when phthalate is causal for IR ",
         "but not for NCO (e.g., hemoglobin, RBC count, MCV, platelets)."),
       identification = paste0(
         "Shared confounding between phthalate and IR / phthalate and NCO; ",
         "no causal link from phthalate -> NCO."),
       estimation = paste0(
         "Weighted linear regression of NCO on high_dehp + M2 covariates ",
         "using nhanes pooled survey design."),
       sensitivity = paste0(
         "Multiple NCO biological systems (red cell line); ",
         "interpret PASS = p>0.05 as supporting absence of unmeasured-confounder bias.")
  ),

  list(method = "Probabilistic Bias Analysis (episensr)",
       reference = "Lash 2009 book / Fox 2021 IJE",
       estimand = paste0(
         "Bias-adjusted OR distribution under specified ranges of ",
         "(a) exposure misclassification, (b) outcome misclassification, ",
         "(c) unmeasured confounder strength."),
       identification = paste0(
         "Trapezoidal prior distributions on sensitivity / specificity / confounder RR; ",
         "non-differential misclassification within type."),
       estimation = paste0(
         "100,000 Monte Carlo iterations; per-iteration bias-adjusted 2x2 reconstruction."),
       sensitivity = paste0(
         "Range of prior parameter trapezoids; compare adj OR to crude OR; ",
         "report % iterations crossing null.")
  ),

  list(method = "Multiverse / Specification curve",
       reference = "Simonsohn 2020 Nat Hum Behav; Steegen 2016 Persp Psych Sci",
       estimand = paste0(
         "Distribution of effect estimates across many defensible specifications ",
         "(exposure sum, cutoff, outcome scale, covariate set, subset)."),
       identification = paste0(
         "Each specification individually defensible; collective summary reflects ",
         "robustness to analyst degrees of freedom."),
       estimation = paste0(
         "3 exposures x 3 cutoffs x 2 outcomes x 2 model specs x 4 subsets = 144 specs; ",
         "each fit with weighted survey design."),
       sensitivity = paste0(
         "Report % positive direction, % significantly positive, ",
         "median estimate per outcome; spec curve plot.")
  )
)

ctat_df <- bind_rows(lapply(ctat_rows, as.data.frame, stringsAsFactors = FALSE))

# Make wider-format CTAT 4-quadrant table per Hernan layout
ctat_long <- ctat_df %>%
  pivot_longer(cols = c(estimand, identification, estimation, sensitivity),
               names_to = "quadrant", values_to = "content") %>%
  mutate(quadrant_label = recode(quadrant,
    estimand       = "(a) Estimand definition",
    identification = "(b) Identification assumption",
    estimation     = "(c) Estimation method",
    sensitivity    = "(d) Sensitivity analysis"))

# Wide (one row per method, 4 quadrant columns) for supplementary S1
ctat_wide <- ctat_df %>%
  select(method, reference, estimand, identification, estimation, sensitivity)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

# Wide table (publication-ready)
write.csv(ctat_wide, "output/tables/ctat_methods_table.csv", row.names = FALSE)

# Long table (machine-readable)
write.csv(ctat_long, "output/tables/ctat_methods_table_long.csv", row.names = FALSE)

cat(sprintf("Methods covered: %d\n", nrow(ctat_df)))
print(ctat_df$method)

cat("\nSaved: output/tables/ctat_methods_table.csv\n")
cat("Saved: output/tables/ctat_methods_table_long.csv\n")

cat("\n========================================\n")
cat(sprintf("28 CTAT 4-quadrant methods table complete (%d methods x 4 quadrants)\n",
            nrow(ctat_df)))
cat("========================================\n")
