# ============================================
# 009 / scripts/31_multiplicity_padjust.R
# Benjamini-Hochberg FDR correction for all manuscript P-values
#
# W16 Round 1 reset — R-Stats M1/M2/M3 fix:
# Manuscript runs ~66 single-metabolite × outcome tests, 9 subgroup
# interaction tests, 8 sensitivity tests — none multiplicity-corrected.
# Implement family-level BH (Benjamini-Hochberg 1995) FDR at q=0.05.
#
# Output: output/tables/multiplicity_corrected.csv
#         (P_raw + P_BH + significance status by family)
# ============================================

suppressPackageStartupMessages({
  library(dplyr)
})

cat("========================================\n")
cat("009 / 31_multiplicity_padjust.R — BH FDR (R-Stats M1/M2/M3 fix)\n")
cat("========================================\n\n")

# ----------------------------------------------------------
# Family 1: 11 single-metab × 2 outcome × 3 model = ≤ 66 tests
# Source: output/tables/logistic_ir_binary.csv + linear_homa_ir.csv
# ----------------------------------------------------------
fam1_rows <- list()

f1_logit <- tryCatch(read.csv("output/tables/logistic_ir_binary.csv",
                              stringsAsFactors = FALSE),
                     error = function(e) NULL)
if (!is.null(f1_logit) && "p_per_SD" %in% names(f1_logit)) {
  d <- f1_logit %>% dplyr::select(exposure, outcome, model, P_raw = p_per_SD)
  d$family <- "Single-metab × IR_binary × 3 models"
  fam1_rows[[length(fam1_rows)+1]] <- d
}

f1_homa <- tryCatch(read.csv("output/tables/linear_homa_ir.csv",
                             stringsAsFactors = FALSE),
                    error = function(e) NULL)
if (!is.null(f1_homa) && "p_per_SD" %in% names(f1_homa)) {
  d <- f1_homa %>% dplyr::select(exposure, outcome, model, P_raw = p_per_SD)
  d$family <- "Single-metab × log-HOMA × 3 models"
  fam1_rows[[length(fam1_rows)+1]] <- d
}

f1_asian <- tryCatch(read.csv("output/tables/logistic_ir_binary_asian.csv",
                              stringsAsFactors = FALSE),
                     error = function(e) NULL)
if (!is.null(f1_asian) && "p_per_SD" %in% names(f1_asian)) {
  d <- f1_asian %>% dplyr::select(exposure, outcome, model, P_raw = p_per_SD)
  d$family <- "Single-metab × IR_binary_asian × 3 models"
  fam1_rows[[length(fam1_rows)+1]] <- d
}

# ----------------------------------------------------------
# Family 2: 9 subgroup interaction tests
# Source: output/tables/subgroup_ir.csv (if exists)
# ----------------------------------------------------------
f2 <- tryCatch(read.csv("output/tables/subgroup_ir.csv",
                        stringsAsFactors = FALSE),
               error = function(e) NULL)
fam2 <- if (!is.null(f2)) {
  pcol <- intersect(c("p","p_interaction","P","P_interaction","p_value"), names(f2))[1]
  if (is.na(pcol)) NULL else {
    sg <- f2 %>% dplyr::filter(!is.na(.data[[pcol]]))
    sgrows <- data.frame(
      exposure = if ("exposure" %in% names(sg)) sg$exposure else NA,
      outcome  = "IR_binary",
      model    = if ("subgroup" %in% names(sg)) sg$subgroup else
                 if ("level" %in% names(sg)) sg$level else
                 if ("group" %in% names(sg)) sg$group else NA,
      P_raw    = sg[[pcol]],
      family   = "Subgroup interactions (9 strata)",
      stringsAsFactors = FALSE
    )
    sgrows
  }
} else NULL

# ----------------------------------------------------------
# Family 3: Sensitivity analyses (8 sets S1-S8)
# Source: output/tables/sensitivity_8sets.csv
# ----------------------------------------------------------
f3 <- tryCatch(read.csv("output/tables/sensitivity_8sets.csv",
                        stringsAsFactors = FALSE),
               error = function(e) NULL)
fam3 <- if (!is.null(f3) && "p_value" %in% names(f3)) {
  d <- data.frame(
    exposure = "Σ-DEHP-z (primary)",
    outcome  = f3$outcome,
    model    = f3$scenario,
    P_raw    = f3$p_value,
    family   = "Sensitivity analyses (S1-S8)",
    stringsAsFactors = FALSE
  )
  d[!is.na(d$P_raw), ]
} else NULL

# ----------------------------------------------------------
# Combine + BH adjust within each family
# ----------------------------------------------------------
all_rows <- list()
if (length(fam1_rows) > 0) all_rows <- c(all_rows, fam1_rows)
if (!is.null(fam2)) all_rows[[length(all_rows)+1]] <- fam2
if (!is.null(fam3)) all_rows[[length(all_rows)+1]] <- fam3

if (length(all_rows) == 0) {
  cat("[WARN] No source CSVs found for multiplicity correction\n")
  write.csv(data.frame(note = "No source CSVs available"),
            "output/tables/multiplicity_corrected.csv", row.names = FALSE)
} else {
  combined <- dplyr::bind_rows(all_rows)
  combined <- combined[!is.na(combined$P_raw), ]

  # BH adjust per family
  out_list <- list()
  for (fam in unique(combined$family)) {
    sub <- combined[combined$family == fam, ]
    sub$rank <- rank(sub$P_raw, ties.method = "min")
    sub$n_tests_in_family <- nrow(sub)
    sub$P_BH <- p.adjust(sub$P_raw, method = "BH")
    sub$BH_q05_sig <- sub$P_BH < 0.05
    out_list[[fam]] <- sub
  }
  out_df <- dplyr::bind_rows(out_list)

  # Round
  out_df$P_raw <- signif(out_df$P_raw, 4)
  out_df$P_BH  <- signif(out_df$P_BH, 4)

  cat(sprintf("Multiplicity-corrected table: %d total tests across %d families\n",
              nrow(out_df), length(unique(out_df$family))))

  by_fam <- out_df %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(n_tests = dplyr::n(),
                     n_raw_p_lt_05 = sum(P_raw < 0.05),
                     n_bh_q_lt_05  = sum(BH_q05_sig),
                     .groups = "drop")
  cat("\n--- Summary by family ---\n")
  print(by_fam)

  cat("\n--- Top 15 by adjusted P ---\n")
  print(head(out_df[order(out_df$P_BH), ], 15), row.names = FALSE)

  if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
  write.csv(out_df, "output/tables/multiplicity_corrected.csv", row.names = FALSE)
  cat(sprintf("\n[OK] output/tables/multiplicity_corrected.csv (rows=%d)\n",
              nrow(out_df)))
}

cat("\n========================================\n")
cat("BH FDR adjustment done.\n")
cat("========================================\n")
