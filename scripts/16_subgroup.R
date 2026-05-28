# ============================================
# 009 / scripts/16_subgroup.R
# 亚组分析 + p-interaction (主 outcome IR binary)
#
# 9 亚组:
#   1. age (20-39 / 40-59 / >=60)
#   2. sex (M / F)
#   3. race (5 levels)
#   4. education (3 levels)
#   5. BMI (<25 / 25-30 / >=30)
#   6. smoke (Never / Ever)
#   7. drink (Yes / No)           # 数据全 NA -> 跳过 + 记日志
#   8. hypertension (Yes / No)
#   9. postmenopausal (女性 Yes/No) # 现数据全 0 -> 跳过 + 记日志
#
# 主暴露: sum_dehp_mol_z (Σ-DEHP) — Plan M13-α 主 phthalate mixture proxy
# Model: M2 完全调整 svyglm + p-interaction via Wald test on exposure × strata
#
# 输出: output/tables/subgroup_ir.csv
#        output/figures/forest_subgroup.png
# ============================================

suppressPackageStartupMessages({
  library(survey); library(dplyr); library(broom); library(ggplot2); library(purrr)
})

cat("========================================\n")
cat("009 / 16_subgroup.R: Subgroup + p-interaction (IR binary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

# Primary exposure (per SD increase, log2-transformed z-score)
PRIMARY_EXP <- "sum_dehp_mol_z"
if (!PRIMARY_EXP %in% names(nhanes_final)) {
  stop("Primary exposure ", PRIMARY_EXP, " not found in nhanes_final")
}

# Build BMI category if missing
if (!"bmi_cat" %in% names(nhanes_final)) {
  nhanes_final$bmi_cat <- cut(nhanes_final$bmi,
                              breaks = c(-Inf, 25, 30, Inf),
                              labels = c("<25","25-30",">=30"),
                              right = FALSE)
}

# Hypertension factor for interaction term
nhanes_final$htn_fac <- factor(ifelse(nhanes_final$hypertension == 1, "Yes", "No"),
                               levels = c("No","Yes"))

# Sex factor
nhanes_final$sex_fac <- factor(ifelse(nhanes_final$sex_male == 1, "Male", "Female"),
                               levels = c("Male","Female"))

# Postmenopausal factor (current data: all 0 -> skip)
n_postmen_yes <- sum(nhanes_final$postmenopausal == 1, na.rm = TRUE)
n_postmen_no  <- sum(nhanes_final$postmenopausal == 0, na.rm = TRUE)
cat(sprintf("Postmenopausal: Yes=%d / No=%d (sparse, see note below)\n",
            n_postmen_yes, n_postmen_no))

# Drink: 检查
drink_table <- table(nhanes_final$drink, useNA = "ifany")
cat("Drink table:\n"); print(drink_table)

# Re-build design with these new columns
design_main2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                          weights = ~wt_pooled, data = nhanes_final, nest = TRUE)

# Base covariate adjust (M2) — drop 当前分层变量
cov_all <- c("age","sex_male","race","education","pir","bmi","waist","smoke","hypertension")

# Helper: 在 sub-design 上跑 svyglm 取 OR
run_subgroup_glm <- function(sub_design, exp_var, drop_cov = NULL, label = "") {
  cov_use <- setdiff(cov_all, drop_cov)
  f <- as.formula(paste("ir_binary ~", exp_var, "+",
                        paste(cov_use, collapse = " + ")))
  out <- tryCatch({
    m <- svyglm(f, design = sub_design, family = quasibinomial())
    n_sub <- nrow(sub_design$variables)
    tt <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::filter(term == exp_var)
    if (nrow(tt) == 0) return(NULL)
    data.frame(
      level = label,
      n = n_sub,
      OR = tt$estimate,
      CI_low = tt$conf.low,
      CI_high = tt$conf.high,
      p_value = tt$p.value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("  [err]", label, ":", conditionMessage(e), "\n")
    data.frame(level=label, n=NA, OR=NA, CI_low=NA, CI_high=NA, p_value=NA,
               stringsAsFactors=FALSE)
  })
  out
}

# Helper: p-interaction via svyglm with exposure × strata, Wald on interaction term
get_pinter <- function(strata_var, exp_var, design = design_main2) {
  cov_use <- setdiff(cov_all, strata_var)
  f_int <- as.formula(paste("ir_binary ~", exp_var, "*", strata_var, "+",
                            paste(cov_use, collapse = " + ")))
  tryCatch({
    m <- svyglm(f_int, design = design, family = quasibinomial())
    # regTermTest on interaction
    int_term <- paste0(exp_var, ":", strata_var)
    pt <- regTermTest(m, as.formula(paste0("~ ", int_term)))
    as.numeric(pt$p)
  }, error = function(e) {
    cat("  [pinter err]", strata_var, ":", conditionMessage(e), "\n")
    NA_real_
  })
}

# ------------------------------------------------------------------
# 1. age_group
# ------------------------------------------------------------------
cat("\n--- 1/9 age_group ---\n")
sg1 <- map_dfr(levels(nhanes_final$age_group), function(lv) {
  sd <- subset(design_main2, age_group == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "age", label = lv)
  r$strata <- "age_group"
  r
})
p1 <- get_pinter("age_group", PRIMARY_EXP)

# ------------------------------------------------------------------
# 2. sex
# ------------------------------------------------------------------
cat("\n--- 2/9 sex ---\n")
sg2 <- map_dfr(levels(nhanes_final$sex_fac), function(lv) {
  sd <- subset(design_main2, sex_fac == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "sex_male", label = lv)
  r$strata <- "sex"
  r
})
p2 <- get_pinter("sex_fac", PRIMARY_EXP)

# ------------------------------------------------------------------
# 3. race
# ------------------------------------------------------------------
cat("\n--- 3/9 race ---\n")
sg3 <- map_dfr(levels(nhanes_final$race), function(lv) {
  sd <- subset(design_main2, race == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "race", label = lv)
  r$strata <- "race"
  r
})
p3 <- get_pinter("race", PRIMARY_EXP)

# ------------------------------------------------------------------
# 4. education
# ------------------------------------------------------------------
cat("\n--- 4/9 education ---\n")
sg4 <- map_dfr(levels(nhanes_final$education), function(lv) {
  sd <- subset(design_main2, education == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "education", label = lv)
  r$strata <- "education"
  r
})
p4 <- get_pinter("education", PRIMARY_EXP)

# ------------------------------------------------------------------
# 5. bmi_cat
# ------------------------------------------------------------------
cat("\n--- 5/9 bmi ---\n")
sg5 <- map_dfr(levels(nhanes_final$bmi_cat), function(lv) {
  sd <- subset(design_main2, bmi_cat == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = c("bmi","waist"), label = lv)
  r$strata <- "bmi"
  r
})
p5 <- get_pinter("bmi_cat", PRIMARY_EXP)

# ------------------------------------------------------------------
# 6. smoke
# ------------------------------------------------------------------
cat("\n--- 6/9 smoke ---\n")
sg6 <- map_dfr(levels(nhanes_final$smoke), function(lv) {
  sd <- subset(design_main2, smoke == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "smoke", label = lv)
  r$strata <- "smoke"
  r
})
p6 <- get_pinter("smoke", PRIMARY_EXP)

# ------------------------------------------------------------------
# 7. drink — 检查可用性
# ------------------------------------------------------------------
cat("\n--- 7/9 drink ---\n")
drink_nonNA <- sum(!is.na(nhanes_final$drink))
if (drink_nonNA >= 100 && length(unique(na.omit(nhanes_final$drink))) >= 2) {
  sg7 <- map_dfr(levels(nhanes_final$drink), function(lv) {
    sd <- subset(design_main2, drink == lv)
    r <- run_subgroup_glm(sd, PRIMARY_EXP, label = lv)
    r$strata <- "drink"
    r
  })
  p7 <- get_pinter("drink", PRIMARY_EXP)
} else {
  cat(sprintf("  [skip] drink (non-NA=%d, levels avail too few) — 数据缺失\n", drink_nonNA))
  sg7 <- data.frame(level = "Yes/No", n = NA, OR = NA, CI_low = NA, CI_high = NA,
                    p_value = NA, strata = "drink", stringsAsFactors = FALSE)
  p7 <- NA_real_
}

# ------------------------------------------------------------------
# 8. hypertension
# ------------------------------------------------------------------
cat("\n--- 8/9 hypertension ---\n")
sg8 <- map_dfr(levels(nhanes_final$htn_fac), function(lv) {
  sd <- subset(design_main2, htn_fac == lv)
  r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "hypertension", label = lv)
  r$strata <- "hypertension"
  r
})
p8 <- get_pinter("htn_fac", PRIMARY_EXP)

# ------------------------------------------------------------------
# 9. postmenopausal — 检查可用性
# ------------------------------------------------------------------
cat("\n--- 9/9 postmenopausal ---\n")
if (n_postmen_yes >= 30 && n_postmen_no >= 30) {
  nhanes_final$postmen_fac <- factor(ifelse(nhanes_final$postmenopausal == 1, "Yes", "No"),
                                     levels = c("No","Yes"))
  # W16 fix: rebuild design with postmen_fac column attached
  design_main2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                            weights = ~wt_pooled, data = nhanes_final, nest = TRUE)
  # Female-only subgroup
  design_f <- subset(design_main2, sex_male == 0)
  sg9 <- map_dfr(c("No","Yes"), function(lv) {
    sd <- subset(design_f, postmen_fac == lv)
    r <- run_subgroup_glm(sd, PRIMARY_EXP, drop_cov = "sex_male", label = lv)
    r$strata <- "postmenopausal (Female)"
    r
  })
  p9 <- tryCatch({
    f_int <- as.formula(paste("ir_binary ~", PRIMARY_EXP, "* postmen_fac +",
                              paste(setdiff(cov_all, "sex_male"), collapse = " + ")))
    m <- svyglm(f_int, design = design_f, family = quasibinomial())
    pt <- regTermTest(m, as.formula(paste0("~ ", PRIMARY_EXP, ":postmen_fac")))
    as.numeric(pt$p)
  }, error = function(e) NA_real_)
} else {
  cat(sprintf("  [skip] postmenopausal (Yes=%d / No=%d, 数据上游清洗将所有女性标 0 -> 全是 No, 跳过)\n",
              n_postmen_yes, n_postmen_no))
  sg9 <- data.frame(level = "Yes/No", n = NA, OR = NA, CI_low = NA, CI_high = NA,
                    p_value = NA, strata = "postmenopausal (Female)", stringsAsFactors = FALSE)
  p9 <- NA_real_
}

# ------------------------------------------------------------------
# 汇总
# ------------------------------------------------------------------
all_subgroups <- bind_rows(sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9) %>%
  dplyr::select(strata, level, n, OR, CI_low, CI_high, p_value)

p_inter_df <- data.frame(
  strata = c("age_group","sex","race","education","bmi","smoke",
             "drink","hypertension","postmenopausal (Female)"),
  p_interaction = c(p1, p2, p3, p4, p5, p6, p7, p8, p9),
  stringsAsFactors = FALSE
)
all_subgroups <- all_subgroups %>% left_join(p_inter_df, by = "strata")

cat("\n--- Full subgroup table ---\n")
print(all_subgroups, row.names = FALSE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(all_subgroups, "output/tables/subgroup_ir.csv", row.names = FALSE)
cat("\nSaved: output/tables/subgroup_ir.csv\n")

# ------------------------------------------------------------------
# Forest plot
# ------------------------------------------------------------------
cat("\n--- Forest plot ---\n")
plot_df <- all_subgroups %>%
  dplyr::filter(!is.na(OR), is.finite(OR), is.finite(CI_low), is.finite(CI_high)) %>%
  mutate(
    label = sprintf("%s: %s", strata, level),
    row_id = row_number()
  )

if (nrow(plot_df) > 0) {
  pinter_label <- p_inter_df %>%
    mutate(p_lab = ifelse(is.na(p_interaction), "p-int: NA",
                          sprintf("p-int: %.3f", p_interaction)))
  plot_df <- plot_df %>% left_join(pinter_label, by = "strata")

  fp <- ggplot(plot_df, aes(x = OR, y = reorder(label, -row_id))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, color = "#2C3E50") +
    geom_point(size = 3, color = "#E74C3C") +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
    labs(x = "OR per SD increase in log2(Sum-DEHP)  [95% CI, log10 scale]",
         y = NULL,
         title = "Subgroup OR for Insulin Resistance (HOMA-IR >= 2.5)",
         subtitle = paste0("Primary exposure: Sum-DEHP z-score. Model M2 (age + sex + race + ",
                           "edu + PIR + BMI + waist + smoke + HTN)")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9))

  if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
  ggsave("output/figures/forest_subgroup.png", fp,
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved: output/figures/forest_subgroup.png\n")
} else {
  cat("[warn] No valid rows for forest plot\n")
}

cat("\n========================================\n")
cat("Subgroup analysis complete.\n")
cat("Notes:\n")
cat("  - drink subgroup: skipped (ALQ111 全 NA in cleaned data, no drink categories preserved)\n")
cat("  - postmenopausal subgroup: skipped (postmenopausal var 全 0 in cleaned data, RHQ060 logic needs revisit)\n")
cat("========================================\n")
