# ============================================
# 009 / scripts/25_bayesian_gcomp.R
# Bayesian g-computation (Keil 2017 EpiMethods; rstanarm/brms)
# Phthalate Sigma-DEHP Q4 vs Q1 -> IR binary
#
# Design:
#   Treatment: high_dehp = sum_dehp_mol_q4 (1 if Q4, 0 if Q1)
#   Outcome:   ir_binary (binomial, probit link by Keil convention)
#   Adjustment: M2 covariate set
#   Sampler: HMC 4 chains x 4000 iter (warmup 1000) — W16 R-Stats C4 fix
#   Inference: counterfactual prediction
#   Reports: ATE, RR, RD with 95% credible interval + ROPE
#
# Output:
#   output/tables/bayesian_gcomp_phth.csv
#   output/figures/bayesian_posterior_plot.png
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(rstanarm); library(posterior)
  library(bayesplot)
})

set.seed(20260523)
options(mc.cores = min(parallel::detectCores(), 4))  # W16 R-Stats C4 fix: 4 chains parallel

cat("========================================\n")
cat("009 / 25 Bayesian g-computation (Sigma-DEHP Q4 vs Q1 -> IR binary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

# ---- Build analytic dataset ----
df <- nhanes_final

q_breaks <- quantile(df$sum_dehp_mol, probs = seq(0, 1, 0.25), na.rm = TRUE)
df$dehp_q <- cut(df$sum_dehp_mol, breaks = q_breaks, include.lowest = TRUE,
                 labels = c("Q1", "Q2", "Q3", "Q4"))

# Subset to Q1 and Q4 for contrast definition
df_q14 <- df %>% filter(dehp_q %in% c("Q1", "Q4")) %>%
  mutate(high_dehp = as.integer(dehp_q == "Q4"))

cat(sprintf("Q1 n=%d, Q4 n=%d (Q1+Q4 total n=%d)\n",
            sum(df_q14$high_dehp == 0),
            sum(df_q14$high_dehp == 1),
            nrow(df_q14)))

# Covariate vars (M2 set)
cov_vars <- c("age", "RIAGENDR", "race", "education", "pir",
              "bmi", "waist", "smoke")
keep_vars <- c("ir_binary", "high_dehp", cov_vars, "wt_pooled")

dat <- df_q14 %>% select(any_of(keep_vars)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

dat$RIAGENDR  <- factor(dat$RIAGENDR)
dat$race      <- factor(dat$race)
dat$education <- factor(dat$education)
dat$smoke     <- factor(dat$smoke)

cat(sprintf("Analytic n (complete cases)      : %d\n", nrow(dat)))
cat(sprintf("IR cases                          : %d (%.1f%%)\n",
            sum(dat$ir_binary), 100*mean(dat$ir_binary)))

# ---- Fit Bayesian outcome model (probit link, IR binary) ----
fmla <- ir_binary ~ high_dehp + age + RIAGENDR + race + education +
                    pir + bmi + waist + smoke

t0 <- Sys.time()
cat("\nFitting Bayesian logistic regression (rstanarm; probit-like via logit)...\n")
fit <- stan_glm(
  formula = fmla,
  data    = dat,
  family  = binomial(link = "logit"),
  prior          = normal(0, 2.5),
  prior_intercept = normal(0, 2.5),
  chains  = 4L,        # W16 R-Stats C4 fix: 4 chains (Methods §2.5 declares 4)
  iter    = 4000L,
  warmup  = 1000L,
  refresh = 0,
  seed    = 20260523L
)
elapsed_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Fit complete in %.1fs\n", elapsed_fit))

# ---- Counterfactual prediction (g-computation) ----
# Strategy: posterior_epred under high_dehp=1 vs high_dehp=0 for each obs.
# ATE = mean( pred1 - pred0 ) over rows; aggregate across posterior draws.

dat_0 <- dat; dat_0$high_dehp <- 0L
dat_1 <- dat; dat_1$high_dehp <- 1L

cat("Generating posterior predictions under counterfactuals...\n")
ep_0 <- posterior_epred(fit, newdata = dat_0)
ep_1 <- posterior_epred(fit, newdata = dat_1)

# Per-draw averages
risk_0 <- rowMeans(ep_0)     # vector length = ndraws (after warmup)
risk_1 <- rowMeans(ep_1)

ate_draws <- risk_1 - risk_0
rr_draws  <- risk_1 / risk_0
rd_draws  <- risk_1 - risk_0
# For OR (cf. ATE on prob scale): not standard for g-comp; we report RR/RD/ATE

# ---- Summaries ----
q_ci <- function(x) c(median = median(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975))

s_risk0 <- q_ci(risk_0)
s_risk1 <- q_ci(risk_1)
s_ate   <- q_ci(ate_draws)
s_rr    <- q_ci(rr_draws)
s_rd    <- q_ci(rd_draws)

# ROPE: practically equivalent zone for ATE = [-0.01, 0.01] (1% absolute risk diff)
rope_low  <- -0.01
rope_high <-  0.01
pct_in_rope <- mean(ate_draws >= rope_low & ate_draws <= rope_high) * 100
pct_above_rope <- mean(ate_draws > rope_high) * 100
pct_below_rope <- mean(ate_draws < rope_low)  * 100
prob_ate_gt0   <- mean(ate_draws > 0) * 100

cat("\n=== Bayesian g-computation results ===\n")
cat(sprintf("  E[Y|do(high_dehp=0)]: %.3f (95%%CrI %.3f-%.3f)\n",
            s_risk0[1], s_risk0[2], s_risk0[3]))
cat(sprintf("  E[Y|do(high_dehp=1)]: %.3f (95%%CrI %.3f-%.3f)\n",
            s_risk1[1], s_risk1[2], s_risk1[3]))
cat(sprintf("  ATE (RD)            : %.4f (95%%CrI %.4f-%.4f)\n",
            s_ate[1], s_ate[2], s_ate[3]))
cat(sprintf("  RR                  : %.3f (95%%CrI %.3f-%.3f)\n",
            s_rr[1], s_rr[2], s_rr[3]))
cat(sprintf("  P(ATE > 0)          : %.1f%%\n", prob_ate_gt0))
cat(sprintf("  ROPE [-0.01,0.01]   : %.1f%% in / %.1f%% above / %.1f%% below\n",
            pct_in_rope, pct_above_rope, pct_below_rope))

# ---- Save table ----
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

res_df <- data.frame(
  estimand = c("E[Y|do(Q1)]", "E[Y|do(Q4)]", "ATE (RD)", "RR", "P(ATE>0)",
               "ROPE [-0.01,0.01] in %", "ROPE above %", "ROPE below %"),
  value    = c(round(s_risk0[1], 3), round(s_risk1[1], 3), round(s_ate[1], 4),
               round(s_rr[1], 3), round(prob_ate_gt0, 1),
               round(pct_in_rope, 1), round(pct_above_rope, 1), round(pct_below_rope, 1)),
  ci_lo    = c(round(s_risk0[2], 3), round(s_risk1[2], 3), round(s_ate[2], 4),
               round(s_rr[2], 3), NA, NA, NA, NA),
  ci_hi    = c(round(s_risk0[3], 3), round(s_risk1[3], 3), round(s_ate[3], 4),
               round(s_rr[3], 3), NA, NA, NA, NA)
)
write.csv(res_df, "output/tables/bayesian_gcomp_phth.csv", row.names = FALSE)
cat("\nSaved: output/tables/bayesian_gcomp_phth.csv\n")

# ---- Posterior plots ----
plot_df <- data.frame(
  draw      = seq_along(ate_draws),
  ate       = as.numeric(ate_draws),
  rr        = as.numeric(rr_draws),
  risk_0    = as.numeric(risk_0),
  risk_1    = as.numeric(risk_1)
)

# ATE density with ROPE
p_ate <- ggplot(plot_df, aes(x = ate)) +
  geom_density(fill = "#1b4965", alpha = 0.5, color = "#1b4965") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(rope_low, rope_high), linetype = "dotted", color = "red") +
  geom_vline(xintercept = s_ate[1], color = "#d7263d", linewidth = 0.8) +
  labs(title = sprintf("Posterior of ATE (Sigma-DEHP Q4 vs Q1 on IR)"),
       subtitle = sprintf("Median = %.4f (95%%CrI %.4f, %.4f); P(ATE>0) = %.1f%%",
                          s_ate[1], s_ate[2], s_ate[3], prob_ate_gt0),
       x = "ATE on risk-difference scale",
       y = "Posterior density") +
  theme_minimal(base_size = 11)

p_rr <- ggplot(plot_df, aes(x = rr)) +
  geom_density(fill = "#7570b3", alpha = 0.5, color = "#7570b3") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = s_rr[1], color = "#d7263d", linewidth = 0.8) +
  labs(title = "Posterior of Risk Ratio",
       subtitle = sprintf("Median = %.3f (95%%CrI %.3f, %.3f)",
                          s_rr[1], s_rr[2], s_rr[3]),
       x = "Risk Ratio (Q4 vs Q1)",
       y = "Posterior density") +
  theme_minimal(base_size = 11)

png("output/figures/bayesian_posterior_plot.png",
    width = 1300, height = 600, res = 150)
gridExtra::grid.arrange(p_ate, p_rr, nrow = 1)
dev.off()

cat("Saved: output/figures/bayesian_posterior_plot.png\n")

# ---- Save full posterior summary RData (for downstream MR comparison) ----
save(fit, risk_0, risk_1, ate_draws, rr_draws, s_ate, s_rr,
     file = "data/processed/bayesian_gcomp_009.RData")

cat("\n========================================\n")
cat("25 Bayesian g-computation complete\n")
cat(sprintf("  -> ATE = %.4f (95%%CrI %.4f, %.4f)\n", s_ate[1], s_ate[2], s_ate[3]))
cat(sprintf("  -> RR  = %.3f (95%%CrI %.3f, %.3f)\n", s_rr[1], s_rr[2], s_rr[3]))
cat(sprintf("  -> P(ATE>0) = %.1f%%\n", prob_ate_gt0))
cat("========================================\n")
