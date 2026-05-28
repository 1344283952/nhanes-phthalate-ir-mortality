# ============================================
# 009 / scripts/04_survey_design.R
# 复杂抽样设计 (NHANES weighted, pooled D-J per NCHS Series 2 No. 190)
# 输入: data/processed/nhanes_final.RData
# 输出: data/processed/nhanes_design.RData
#
# Plan M13-α: 4 个 design 对应 4 个 cohort stack:
#   design_main     → Stack 1 主分析 (N=2,239)
#   design_mortality → Stack 3 Mortality eligible (N=2,238)
#   design_t2d_prog  → Stack 2 T2D progression (N=1,247 prediabetes)
#   design_phth_pfas → Stack 4 Phthalate+PFAS (N=1,082)
# ============================================

library(survey); library(dplyr)

cat("========================================\n")
cat("009 / 复杂抽样设计 (4 stack, Plan M13-α)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
options(survey.lonely.psu = "adjust")

# ----------------------------------------------------------
# Design 1: Stack 1 主分析 (PHTHTE D-J + fasting + 排除)
# ----------------------------------------------------------
design_main <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~wt_pooled, data = nhanes_final,
  nest = TRUE
)
cat(sprintf("design_main (Stack 1 主): n=%d, weights=wt_pooled\n", nrow(nhanes_final)))

# Stack 3 Mortality (subset of Stack 1)
design_mortality <- subset(design_main, !is.na(ELIGSTAT) & ELIGSTAT == 1)
cat(sprintf("design_mortality (Stack 3): n=%d\n", nrow(design_mortality$variables)))

# ----------------------------------------------------------
# Design 2: Stack 2 T2D progression (prediabetes subset)
# ----------------------------------------------------------
if (exists("df_t2d_prog") && nrow(df_t2d_prog) > 0) {
  design_t2d_prog <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = ~wt_pooled, data = df_t2d_prog,
    nest = TRUE
  )
  cat(sprintf("design_t2d_prog (Stack 2 prediabetes): n=%d\n", nrow(df_t2d_prog)))
} else {
  design_t2d_prog <- NULL
}

# ----------------------------------------------------------
# Design 3: Stack 4 Phthalate+PFAS (cross-program synergy with 007)
# ----------------------------------------------------------
if (exists("df_phth_pfas") && nrow(df_phth_pfas) > 0) {
  design_phth_pfas <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = ~wt_pooled, data = df_phth_pfas,
    nest = TRUE
  )
  cat(sprintf("design_phth_pfas (Stack 4 Phthalate+PFAS): n=%d\n", nrow(df_phth_pfas)))
} else {
  design_phth_pfas <- NULL
}

# ----------------------------------------------------------
# W16 R-NHANES MAJOR 4 fix: Survey degrees of freedom report
# Per NCHS Series 2 No. 190 §6 + Lumley 2010 §4, complex-survey CIs should
# use t with df = nPSU - nStrata. degf() reports this for each design.
# ----------------------------------------------------------
cat("\n--- Survey degrees of freedom (per design) ---\n")
degf_rows <- list()

safe_degf <- function(d, label) {
  if (is.null(d)) return(NULL)
  v <- tryCatch(survey::degf(d), error = function(e) NA_integer_)
  cat(sprintf("  %-25s degf = %s\n", label, as.character(v)))
  data.frame(design = label, n_rows = nrow(d$variables), degf = v,
             stringsAsFactors = FALSE)
}

degf_rows[[1]] <- safe_degf(design_main,      "design_main")
degf_rows[[2]] <- safe_degf(design_mortality, "design_mortality")
degf_rows[[3]] <- safe_degf(design_t2d_prog,  "design_t2d_prog")
degf_rows[[4]] <- safe_degf(design_phth_pfas, "design_phth_pfas")

degf_df <- dplyr::bind_rows(Filter(Negate(is.null), degf_rows))
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(degf_df, "output/tables/survey_degf.csv", row.names = FALSE)
cat("→ output/tables/survey_degf.csv\n")

# Save
save(design_main, design_mortality, design_t2d_prog, design_phth_pfas,
     file = "data/processed/nhanes_design.RData")
cat("\n已保存 data/processed/nhanes_design.RData (4 design)\n")
cat("========================================\n")
