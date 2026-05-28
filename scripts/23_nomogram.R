# ============================================
# 009 / scripts/23_nomogram.R
# Nomogram (rms package, Harrell 2015)
# 输入: data/processed/tripod_ir_models.RData
# 输出: output/figures/nomogram_phth_ir.png + nomogram_phth_ir.pdf
#       output/tables/nomogram_points_table.csv
#
# 包含: 8 Phthalate metabolites + 7 key covariates (age, sex, race, edu, pir, bmi, waist)
# 输出: nomogram + 5-tier risk strata (very-low / low / medium / high / very-high)
# ============================================

suppressPackageStartupMessages({
  library(dplyr); library(rms)
})

set.seed(20260523)

cat("========================================\n")
cat("009 / 23 Nomogram + 5-tier risk strata\n")
cat("========================================\n\n")

load("data/processed/tripod_ir_models.RData")
cat(sprintf("Cohort N = %d, IR=1: %d\n", nrow(df), sum(y==1)))

# 7 key covariates + 8 Phthalate metabolites (task requirement)
key_cov <- c("age","sex_male","race","education","pir","bmi","waist")
phth_8  <- c("URXMEP_z","URXMBP_z","URXMIB_z","URXMZP_z",
             "URXMHP_z","URXMHH_z","URXMOH_z","URXECP_z")
preds <- c(key_cov, phth_8)

# W16 SA-3 fix: Attach human-readable labels via Hmisc::label() so nomogram
# shows "Age (years) / Male sex / MEP (z-score) ..." instead of raw column names.
suppressPackageStartupMessages(library(Hmisc))
Hmisc::label(df$age)      <- "Age (years)"
Hmisc::label(df$sex_male) <- "Male sex"
Hmisc::label(df$race)     <- "Race/Ethnicity"
Hmisc::label(df$education) <- "Education"
Hmisc::label(df$pir)      <- "PIR"
Hmisc::label(df$bmi)      <- "BMI (kg/m2)"
Hmisc::label(df$waist)    <- "Waist (cm)"
Hmisc::label(df$URXMEP_z) <- "MEP (z-score)"
Hmisc::label(df$URXMBP_z) <- "MnBP (z-score)"
Hmisc::label(df$URXMIB_z) <- "MiBP (z-score)"
Hmisc::label(df$URXMZP_z) <- "MBzP (z-score)"
Hmisc::label(df$URXMHP_z) <- "MEHP (z-score)"
Hmisc::label(df$URXMHH_z) <- "MEHHP (z-score)"
Hmisc::label(df$URXMOH_z) <- "MEOHP (z-score)"
Hmisc::label(df$URXECP_z) <- "MECPP (z-score)"

# W16 Wave 3 SA-A6 fix B: Rename factor levels BEFORE datadist so we can keep
# abbrev = FALSE on nomogram() and avoid abbreviate()'s vowel-drop bug that
# turns "High school" into "Hghs". We use readable short names that fit the
# nomogram axis width, then expand them in the figure footer legend.
race_short <- c("Non-Hispanic White" = "N-HW",
                "Non-Hispanic Black" = "N-HB",
                "Mexican American"   = "MxcA",
                "Other Hispanic"     = "OthH",
                "Other Race"         = "OthR")
edu_short  <- c("Less than HS"     = "LessHS",
                "High school"      = "HighS",
                "College or above" = "ColPlus")
df$race      <- factor(race_short[as.character(df$race)],
                       levels = race_short)
df$education <- factor(edu_short[as.character(df$education)],
                       levels = edu_short)
# Re-attach labels lost by factor() rebuild
Hmisc::label(df$race)      <- "Race/Ethnicity"
Hmisc::label(df$education) <- "Education"

# rms 要求 datadist
dd <- rms::datadist(df[, preds])
options(datadist = "dd")

# W16 Wave 3 SA-A6 fix A: Clip nomogram axes to clinically plausible ranges.
# rms::datadist defaults to the full data range, which for BMI yields 10 -> 140
# (extreme outliers) and for waist 60 -> 170. NHANES adult plausible bounds are
# BMI 18-50 (covers ~99% of US adults) and waist 60-160 cm. This shortens the
# nomogram axes dramatically and removes the implausible 140 BMI upper bound
# called out by R-Figure v2 round-2 review.
dd$limits["Low:prediction",  "bmi"] <- 18
dd$limits["High:prediction", "bmi"] <- 50
dd$limits["Low",             "bmi"] <- 18
dd$limits["High",            "bmi"] <- 50
dd$limits["Low:prediction",  "waist"] <- 60
dd$limits["High:prediction", "waist"] <- 160
dd$limits["Low",             "waist"] <- 60
dd$limits["High",            "waist"] <- 160
cat("\n[axis clip] BMI -> [18, 50]; Waist -> [60, 160]\n")

f_nom <- as.formula(paste("ir_binary ~", paste(preds, collapse = " + ")))
fit_nom <- rms::lrm(f_nom, data = df, x = TRUE, y = TRUE)
cat("\nlrm fit summary:\n")
print(fit_nom)

# Discrimination summary
disc <- fit_nom$stats
cat(sprintf("\nC-statistic = %.3f\n", disc["C"]))

# ============================================
# Nomogram (W16 SA-3 fix: abbrev = TRUE for race + education to prevent overlap)
# rms 不支持 force-label-on-separate-rows for multi-level factors when many
# levels have similar score positions; the practical fix is abbreviation +
# include legend below + wider canvas.
# ============================================
nom <- rms::nomogram(
  fit_nom,
  fun = plogis,
  fun.at = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80),
  funlabel = "Predicted IR risk",
  abbrev = FALSE  # W16 Wave 3 SA-A6 fix B: factor levels were manually renamed
                  # above (race -> N-HW/N-HB/MxcA/OthH/OthR; education ->
                  # LessHS/HighS/ColPlus) so we do NOT need rms abbreviate()'s
                  # vowel-drop algo (which previously turned "High school" into
                  # "Hghs"). Custom short names are guaranteed readable.
)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# W16 SA-3 fix: race axis labels non-overlapping via custom short names
# (W16 Wave 3 SA-A6 fix B: factor levels manually renamed; legend now uses
# the rename map directly instead of rms::abbreviate()).
#   - widen overall canvas (14 -> 16 inch)
#   - xfrac 0.30 -> 0.35 (more room for labels)
#   - cex.var 0.7 -> 0.85 (legible)
#   - cex.axis 0.6 -> 0.75
#   - Legend below shows full names
legend_text  <- c(
  paste0("Race: ",
         paste0(race_short, " = ", names(race_short), collapse = "; ")),
  paste0("Education: ",
         paste0(edu_short, " = ", names(edu_short), collapse = "; "))
)

png("output/figures/nomogram_phth_ir.png",
    width = 16, height = 10, units = "in", res = 300)
par(mar = c(8, 8, 3, 3))  # bottom margin bumped to 8 for legend
plot(nom, xfrac = 0.40, lmgp = 0.25,
     cex.axis = 0.75, cex.var = 0.9, label.every = 1)
title(main = "Nomogram: Phthalate-augmented IR risk (NHANES 2005-2018)",
      sub = sprintf("C-statistic = %.3f; 5-tier risk stratification very-low to very-high", disc["C"]),
      cex.main = 1.3)
mtext(legend_text, side = 1, line = c(5, 6), at = 0.05, adj = 0,
      cex = 0.75, font = 3)
dev.off()

pdf("output/figures/nomogram_phth_ir.pdf", width = 16, height = 10)
par(mar = c(8, 8, 3, 3))
plot(nom, xfrac = 0.40, lmgp = 0.25,
     cex.axis = 0.75, cex.var = 0.9, label.every = 1)
title(main = "Nomogram: Phthalate-augmented IR risk (NHANES 2005-2018)",
      sub = sprintf("C-statistic = %.3f; 5-tier risk stratification very-low to very-high", disc["C"]),
      cex.main = 1.3)
mtext(legend_text, side = 1, line = c(5, 6), at = 0.05, adj = 0,
      cex = 0.75, font = 3)
dev.off()

# ============================================
# Points table (long format for paper Supp)
# nomogram() returns list; each component has $points + reference vals
# ============================================
pt_rows <- list()
for (vn in names(nom)) {
  comp <- nom[[vn]]
  if (is.list(comp) && !is.null(comp$points)) {
    refs <- comp[[1]]
    pts  <- comp$points
    L <- min(length(refs), length(pts))
    if (L > 0) {
      pt_rows[[vn]] <- data.frame(
        variable = vn,
        ref_value = as.character(refs[seq_len(L)]),
        points = as.numeric(pts[seq_len(L)])
      )
    }
  }
}
pt_df <- do.call(rbind, pt_rows)
rownames(pt_df) <- NULL

# Total-points → predicted prob mapping (5 tiers)
df$lin_pred <- predict(fit_nom, type = "lp")
df$pred_prob <- plogis(df$lin_pred)
qs <- quantile(df$pred_prob, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
df$risk_tier <- cut(df$pred_prob, breaks = qs, include.lowest = TRUE,
                    labels = c("very-low","low","medium","high","very-high"))
risk_tbl <- df %>% group_by(risk_tier) %>%
  summarise(n = n(),
            ir_observed = sum(ir_binary == 1),
            ir_rate = mean(ir_binary == 1),
            mean_pred = mean(pred_prob),
            min_pred = min(pred_prob),
            max_pred = max(pred_prob),
            .groups = "drop")
cat("\n--- Risk tier 验证 ---\n")
print(risk_tbl, digits = 3)

write.csv(pt_df,    "output/tables/nomogram_points_table.csv", row.names = FALSE)
write.csv(risk_tbl, "output/tables/nomogram_risk_tiers.csv",  row.names = FALSE)

cat("\n输出: output/figures/nomogram_phth_ir.png\n")
cat("输出: output/figures/nomogram_phth_ir.pdf\n")
cat("输出: output/tables/nomogram_points_table.csv\n")
cat("输出: output/tables/nomogram_risk_tiers.csv\n")
cat("========================================\n23 完成\n========================================\n")
