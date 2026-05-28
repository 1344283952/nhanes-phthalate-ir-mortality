# ============================================
# 009 / scripts/03_clean_data.R
# 输入: data/processed/nhanes_raw_merged.RData (含 PHTHTE D-J + INS/GLU fasting + PFAS for Stack 4)
# 输出: data/processed/nhanes_final.RData
#       (nhanes_final + df_mort + df_t2d_prog + df_phth_pfas + scale_attrs)
#
# Plan M13-α 关键 logic:
#   清洗加项: HEPB/HEPC 排除 + 既存糖尿病排除 + Cotinine + Fish freq + 蛋白质 + HSCRP
#   Phthalate 处理: LOD/√2 + log2 + z-score + Σ-DEHP/Σ-HMW/Σ-LMW + creatinine 调整 + creatinine-adj 双套
#   IR outcome:
#     - HOMA-IR continuous = LBXIN * LBXGLU / 405 (fasting only)
#     - IR binary HOMA ≥ 2.5 (McAuley 2001 ATP III)
#     - 敏感性 cutoff HOMA ≥ 3.6 (Yamada 2013 Asian-adapted)
#   Fasting subsample 严格: PHAFSTHR ≥ 8.5 h
#   既存糖尿病排除: DIQ010==1 OR HbA1c≥6.5 OR diabetes med (IR 关注 *前期*)
#   MetALD/FIB-4 as co-mediator (复用 007)
#
# Cohort 切法:
#   Stack 1 主: PHTHTE ∩ fasting ∩ ≥20 + 排除, N ≈ 1,974 (W10 quick 实查)
#   Stack 2 T2D progression: prediabetes (HbA1c 5.7-6.4 + fasting glucose 100-125)
#   Stack 3 Mortality eligible: Stack 1 ∩ LMF
#   Stack 4 Phthalate+PFAS: Stack 1 ∩ PFAS subsample (~800 估)
#
# 复用 templates/_shared/{fib4.R, pooled_weight.R, var_aliases.R, cycle_config.R,
#                          diabetes_define.R, covariates_config.R}
#
# ===========================================================
# W16 Round 1 reset — Wave 1 SA-1 fixes (R-NHANES + R-Stats):
#
# A1. Weight chain (R-NHANES C1):
#     Stack 1 cohort = MEC ∩ PHTHTE ∩ Fasting. Per NCHS Series 2 No. 190 §4
#     (Least Common Denominator rule), the appropriate base weight is the
#     **most restrictive subsample weight** = WTSAF2YR (fasting subsample),
#     NOT WTMEC2YR (MEC exam) nor WTSB2YR (PHTHTE subsample, less restrictive).
#     Implementation: pooled_saf_weight() (templates/_shared/pooled_weight.R L58-66)
#     applies WTSAF2YR × (cycle_years / total_years) cycle-by-cycle.
#
# A2. PERMTH_INT cap (R-NHANES C2):
#     NCHS LMF 2019 maximum theoretical follow-up for cycle D (2005-2006)
#     = 15 yr = 180 mo; cycle J (2017-2018) max ≈ 36 mo. Any PERMTH_INT > 200 mo
#     is biologically impossible for D-J cohort + likely SEQN merge contamination
#     or unstripped NCHS sentinel (9997 = LTFU). Cap permth at 200 mo as ceiling.
#
# A3. Cycle A/B/C exclusion documentation (R-NHANES C3):
#     Cycles A (1999-2000), B (2001-2002), C (2003-2004) excluded because:
#     - LBXIN serum-insulin assay protocol changed at cycle D (RIA → ECLIA,
#       Roche Diagnostics; NCHS Lab Procedure Manual 2005-2006 revision),
#       introducing assay heterogeneity incompatible with pooled HOMA-IR.
#     - BIOPRO module field coverage inconsistent pre-2005 (templates/_shared/
#       cycle_config.R: biochem panel changed).
#     - PHTHTE D-J 是 NHANES standard 集 (creatinine adjustment URXUCR fully
#       available D-onward).
#     - PLOS-Bio 2025 (Suchak et al., DOI 10.1371/journal.pbio.3003152)
#       cycle 完整性 audit 要求 cycle exclusion 显式 justify.
#
# A4. Σ-DEHP-mol weights (R-NHANES MAJOR 2):
#     Use PubChem 4-digit decimals (MEHP 278.34, MEHHP 294.34, MEOHP 292.32,
#     MECPP 308.32) not integer rounded weights (278/294/292/308).
# ===========================================================

library(dplyr); library(tidyr); library(purrr)

source("../../templates/_shared/fib4.R")
source("../../templates/_shared/pooled_weight.R")
source("../../templates/_shared/var_aliases.R")
source("../../templates/_shared/cycle_config.R")
source("../../templates/_shared/diabetes_define.R")
source("../../templates/_shared/covariates_config.R")

set.seed(20260523)  # 009 启动日 seed

cat("========================================\n")
cat("009 / 清洗 Phthalate × IR + Mortality cohort (Plan M13-α)\n")
cat("========================================\n\n")

load("data/processed/nhanes_raw_merged.RData")
cat(sprintf("原始合并: %d 行 x %d 列\n\n", nrow(nhanes_all), ncol(nhanes_all)))

flow <- list()
log_flow <- function(label, n) {
  cat(sprintf("  [流程] %-70s n = %d\n", label, n))
  flow[[length(flow) + 1]] <<- data.frame(step = label, n = n, stringsAsFactors = FALSE)
}
log_flow("原始合并 (9 cycles, C-J + P_)", nrow(nhanes_all))

# Util
na_codes <- function(x, codes) ifelse(x %in% codes, NA, x)
zero_to_na <- function(x) ifelse(!is.na(x) & x == 0, NA, x)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Step 0: 跨周期变量名 normalize
nhanes_all <- apply_var_aliases(nhanes_all)
cat("apply_var_aliases: canonical 列已生成\n\n")

# ----------------------------------------------------------
# Step 1: 入选 — age ≥ 20 + 非妊娠
# ----------------------------------------------------------
df <- nhanes_all %>% filter(RIDAGEYR >= 20)
log_flow("Age >= 20", nrow(df))
df <- df %>% filter(is.na(RIDEXPRG) | RIDEXPRG != 1)
log_flow("+ not pregnant", nrow(df))

# ----------------------------------------------------------
# Step 2: PHTHTE 6 核心 metabolites 任一非缺 (Stack 1 主入选)
# NHANES Phthalate 字段: URXMEP / URXMBP / URXMIB / URXMZP / URXMHH / URXMHP / URXMOH / URXECP 等
# ----------------------------------------------------------
phth_cols_core <- c("URXMEP",    # MEP (DEP metabolite, cosmetics)
                    "URXMBP",    # MnBP (DBP metabolite, PVC)
                    "URXMIB",    # MiBP (DiBP)
                    "URXMZP",    # MBzP (BBzP, flooring)
                    "URXMHP",    # MEHP (DEHP M1)
                    "URXMHH",    # MEHHP (DEHP M2)
                    "URXMOH",    # MEOHP (DEHP M3)
                    "URXECP",    # MECPP (DEHP M4)
                    "URXMCOH",   # MCOH
                    "URXMNP")    # MNP
phth_cols_avail <- intersect(phth_cols_core, names(df))
phth_cols_avail <- phth_cols_avail[sapply(phth_cols_avail, function(c) sum(!is.na(df[[c]])) > 0)]
cat(sprintf("Phthalate metabolites 可用: %s\n", paste(phth_cols_avail, collapse = ", ")))

df <- df %>% filter(rowSums(!is.na(across(all_of(phth_cols_avail)))) > 0)
log_flow(sprintf("+ Phthalate 核心 (%d 化合物) 任一非缺", length(phth_cols_avail)), nrow(df))

# ----------------------------------------------------------
# Step 3: Phthalate LOD/√2 substitution + log2-transform + z-score
# ----------------------------------------------------------
for (col in phth_cols_avail) {
  vals <- df[[col]]
  if (sum(!is.na(vals)) < 50) next
  lod_est <- min(vals[vals > 0], na.rm = TRUE)
  imp_val <- lod_est / sqrt(2)
  df[[paste0(col, "_imp")]] <- ifelse(is.na(vals) | vals <= 0, imp_val, vals)
  df[[paste0(col, "_log2")]] <- log2(df[[paste0(col, "_imp")]])
  df[[paste0(col, "_z")]] <- as.numeric(scale(df[[paste0(col, "_log2")]]))
}

# Σ-DEHP (4 metabolites mole-weighted): MEHP + MEHHP + MEOHP + MECPP
dehp_metab_cols <- intersect(c("URXMHP_imp","URXMHH_imp","URXMOH_imp","URXECP_imp"), names(df))
if (length(dehp_metab_cols) >= 3) {
  # PubChem 4-digit mole weights (R-NHANES MAJOR 2 fix W16; CIDs reconciled W18 via live PubChem REST API):
  #   MEHP   278.34 g/mol (PubChem CID 20393)
  #   MEHHP  294.34 g/mol (PubChem CID 170295)
  #   MEOHP  292.32 g/mol (PubChem CID 119096)
  #   MECPP  308.32 g/mol (PubChem CID 148386)
  df$sum_dehp_mol <- with(df,
    (`URXMHP_imp` %||% 0) / 278.34 +
    (`URXMHH_imp` %||% 0) / 294.34 +
    (`URXMOH_imp` %||% 0) / 292.32 +
    (`URXECP_imp` %||% 0) / 308.32)
  df$sum_dehp_mol_log2 <- log2(pmax(df$sum_dehp_mol, 0.001))
  df$sum_dehp_mol_z <- as.numeric(scale(df$sum_dehp_mol_log2))
}

# Σ-HMW (high MW): MEHP/MEHHP/MEOHP/MECPP/MCOH/MNP
hmw_cols <- intersect(c("URXMHP_imp","URXMHH_imp","URXMOH_imp","URXECP_imp","URXMCOH_imp","URXMNP_imp"), names(df))
if (length(hmw_cols) >= 3) {
  df$sum_hmw <- rowSums(df[, hmw_cols], na.rm = TRUE)
  df$sum_hmw_log2 <- log2(pmax(df$sum_hmw, 0.01))
  df$sum_hmw_z <- as.numeric(scale(df$sum_hmw_log2))
}

# Σ-LMW (low MW): MEP/MBP/MIB/MZP
lmw_cols <- intersect(c("URXMEP_imp","URXMBP_imp","URXMIB_imp","URXMZP_imp"), names(df))
if (length(lmw_cols) >= 3) {
  df$sum_lmw <- rowSums(df[, lmw_cols], na.rm = TRUE)
  df$sum_lmw_log2 <- log2(pmax(df$sum_lmw, 0.01))
  df$sum_lmw_z <- as.numeric(scale(df$sum_lmw_log2))
}

# HMW/LMW ratio (exposure source proxy)
if ("sum_hmw" %in% names(df) && "sum_lmw" %in% names(df)) {
  df$hmw_lmw_ratio <- df$sum_hmw / df$sum_lmw
  df$hmw_lmw_ratio[!is.finite(df$hmw_lmw_ratio)] <- NA
}

cat(sprintf("\nPhthalate 衍生暴露: sum_dehp_mol (n=%d), sum_hmw (n=%d), sum_lmw (n=%d)\n",
            sum(!is.na(df$sum_dehp_mol %||% rep(NA, nrow(df)))),
            sum(!is.na(df$sum_hmw %||% rep(NA, nrow(df)))),
            sum(!is.na(df$sum_lmw %||% rep(NA, nrow(df))))))

# Creatinine adjustment (URXUCR available?)
if ("URXUCR" %in% names(df)) {
  # Creatinine-adjusted Phthalate (per gram creatinine)
  for (col in phth_cols_avail) {
    df[[paste0(col, "_cr")]] <- ifelse(!is.na(df[[col]]) & !is.na(df$URXUCR) & df$URXUCR > 0,
                                       df[[col]] / df$URXUCR * 100,  # per 100 mg/dL creatinine
                                       NA)
  }
  cat(sprintf("Creatinine adjustment 完成 (URXUCR avail N=%d)\n", sum(!is.na(df$URXUCR))))
}

# ----------------------------------------------------------
# Step 4: Fasting subsample (PHAFSTHR ≥ 8.5 h, M13 关键 filter)
# ----------------------------------------------------------
df$fasting_hours <- with(df,
  ifelse(!is.na(PHAFSTHR), PHAFSTHR + ifelse(!is.na(PHAFSTMN), PHAFSTMN/60, 0), NA))
df_fasting <- df %>% filter(fasting_hours >= 8.5)
log_flow(sprintf("+ Fasting subsample (PHAFSTHR >= 8.5 h)"), nrow(df_fasting))
df <- df_fasting

# ----------------------------------------------------------
# Step 5: HOMA-IR 完整 (LBXIN + LBXGLU)
# ----------------------------------------------------------
df <- df %>% filter(!is.na(LBXIN), !is.na(LBXGLU))
log_flow("+ HOMA-IR ready (LBXIN + LBXGLU)", nrow(df))

df$homa_ir <- (df$LBXIN * df$LBXGLU) / 405
df$homa_ir_log <- log(pmax(df$homa_ir, 0.01))
df$ir_binary <- as.integer(df$homa_ir >= 2.5)        # McAuley 2001 ATP III
df$ir_binary_asian <- as.integer(df$homa_ir >= 3.6)  # Yamada 2013 Asian-adapted (sensitivity)

cat(sprintf("\nHOMA-IR computed: median=%.2f, IR binary (≥2.5)=%d (%.1f%%), IR binary (≥3.6)=%d (%.1f%%)\n",
            median(df$homa_ir, na.rm = TRUE),
            sum(df$ir_binary == 1, na.rm = TRUE),
            100 * mean(df$ir_binary == 1, na.rm = TRUE),
            sum(df$ir_binary_asian == 1, na.rm = TRUE),
            100 * mean(df$ir_binary_asian == 1, na.rm = TRUE)))

# ----------------------------------------------------------
# Step 6: HEPB/HEPC 排除
# ----------------------------------------------------------
df$hep_excl <- (df$hbsag %in% 1) | (df$hcv %in% 1)
df$hep_excl[is.na(df$hep_excl)] <- FALSE
n_hep_pos <- sum(df$hep_excl, na.rm = TRUE)
df <- df %>% filter(!hep_excl)
log_flow(sprintf("+ 排除 HEPB/HEPC 阳性 (n=%d)", n_hep_pos), nrow(df))

# ----------------------------------------------------------
# Step 7: 既存糖尿病排除 (M13 关键 — IR 关注 *前期*)
# DIQ010==1 (self-report) OR HbA1c ≥ 6.5 OR diabetes med
# ----------------------------------------------------------
df$DIQ010 <- na_codes(df$DIQ010, c(7, 9))
df$diabetes_self <- as.integer(!is.na(df$DIQ010) & df$DIQ010 == 1)
df$hba1c_high <- as.integer(!is.na(df$LBXGH) & df$LBXGH >= 6.5)

# Diabetes med (from rx_all RXQ_RX merge)
if (exists("rx_all") && nrow(rx_all) > 0 && "RXDDRUG" %in% names(rx_all)) {
  rx_dm <- detect_antidiabetic_rx(rx_all)
  df <- df %>% left_join(rx_dm, by = "SEQN") %>%
    mutate(rx_diabetes_yes = ifelse(is.na(rx_diabetes_yes), FALSE, rx_diabetes_yes))
}
df$diabetes_existing <- with(df,
  as.integer(diabetes_self == 1 | hba1c_high == 1 | (rx_diabetes_yes %||% FALSE)))

n_dm <- sum(df$diabetes_existing == 1, na.rm = TRUE)
df <- df %>% filter(diabetes_existing != 1 | is.na(diabetes_existing))
log_flow(sprintf("+ 排除既存糖尿病 (DIQ010==1 OR HbA1c>=6.5 OR Rx, n=%d)", n_dm), nrow(df))

# ----------------------------------------------------------
# Step 8: Active cancer 排除
# ----------------------------------------------------------
df$MCQ220 <- na_codes(df$MCQ220, c(7, 9))
n_cancer <- sum(df$MCQ220 == 1, na.rm = TRUE)
df <- df %>% filter(!(MCQ220 %in% 1))
log_flow(sprintf("+ 排除 self-reported cancer (n=%d)", n_cancer), nrow(df))

# ----------------------------------------------------------
# Step 9: 血压 + sex_male + other 协变量
# ----------------------------------------------------------
df <- df %>% mutate(
  BPXSY1 = zero_to_na(BPXSY1), BPXSY2 = zero_to_na(BPXSY2),
  BPXSY3 = zero_to_na(BPXSY3), BPXSY4 = zero_to_na(BPXSY4),
  BPXDI1 = zero_to_na(BPXDI1), BPXDI2 = zero_to_na(BPXDI2),
  BPXDI3 = zero_to_na(BPXDI3), BPXDI4 = zero_to_na(BPXDI4)
) %>% rowwise() %>% mutate(
  sbp = mean(c(BPXSY2, BPXSY3, BPXSY4), na.rm = TRUE),
  dbp = mean(c(BPXDI2, BPXDI3, BPXDI4), na.rm = TRUE)
) %>% ungroup() %>% mutate(
  sbp = ifelse(is.nan(sbp), NA, sbp),
  dbp = ifelse(is.nan(dbp), NA, dbp),
  sex_male = as.integer(RIAGENDR == 1),
  albumin_gdl = if ("alb_gl" %in% names(df)) alb_gl / 10 else NA_real_
)

# ----------------------------------------------------------
# Step 10: FIB-4 (作 co-mediator, 复用 007)
# ----------------------------------------------------------
if (all(c("ast_unl","alt_unl","LBXPLTSI") %in% names(df))) {
  df$fib4 <- calc_fib4(df$RIDAGEYR, df$ast_unl, df$alt_unl, df$LBXPLTSI)
  df$fib4_log <- log(pmax(df$fib4, 0.01))
  df$fib4_advanced <- as.integer(df$fib4 >= 2.67)
}

# ----------------------------------------------------------
# Step 11: Mortality outcome
# ----------------------------------------------------------
df <- df %>% mutate(
  mort_allcause = ifelse(!is.na(MORTSTAT), as.integer(MORTSTAT == 1), NA),
  mort_cm = ifelse(!is.na(MORTSTAT),
                   as.integer(MORTSTAT == 1 & UCOD_LEADING %in% c(1, 5, 7)), NA),
  # DM-cause specific mortality (M13 独占)
  mort_dm = ifelse(!is.na(MORTSTAT),
                   as.integer(MORTSTAT == 1 & UCOD_LEADING == 7), NA),
  # 优先 PERMTH_INT (range 0-180 months 正常, NCHS LMF 2019 standard);
  # PERMTH_EXM range 0-80 异常 (truncated, sub-agent A 2026-05-23 实查发现)
  permth = ifelse(!is.na(PERMTH_INT), PERMTH_INT, PERMTH_EXM)
)

# W16 Round 1 reset — R-NHANES C2 fix:
# Cap PERMTH at 200 mo (biological ceiling for 2005-2018 cohort).
# Cycle D (2005) max LMF 2019 follow-up = 180 mo; J (2017) = 36 mo.
# Values > 200 indicate SEQN merge contamination or unstripped NCHS sentinel
# (9997 = LTFU); set to NA so they drop from Cox cohort.
n_permth_outlier <- sum(!is.na(df$permth) & df$permth > 200)
if (n_permth_outlier > 0) {
  cat(sprintf("[R-NHANES C2 fix] PERMTH > 200 mo outliers detected: n=%d (biological ceiling for D-J cohort = 180); set to NA\n",
              n_permth_outlier))
}
df$permth <- ifelse(!is.na(df$permth) & df$permth > 200, NA_integer_, df$permth)
cat(sprintf("[diagnostic] PERMTH range after cap: [%s, %s] mo\n",
            ifelse(any(!is.na(df$permth)), as.character(min(df$permth, na.rm = TRUE)), "NA"),
            ifelse(any(!is.na(df$permth)), as.character(max(df$permth, na.rm = TRUE)), "NA")))

# ----------------------------------------------------------
# Step 12: 协变量 (复用 007 + Phthalate-specific)
# ----------------------------------------------------------
df <- df %>% mutate(
  DMDEDUC2 = na_codes(DMDEDUC2, c(7, 9)),
  DMDMARTL = na_codes(DMDMARTL, c(77, 99)),
  SMQ020   = na_codes(SMQ020, c(7, 9)),
  # ALQ101 是 NHANES D-J 全 cycle 跨周期通用字段 (lifetime ever-drinker)
  # ALQ111 只 2017-2018+ 有 (CDC 后期重命名), 早期 cycle 全 NA
  # Sub-agent A/B/C 2026-05-23 实查 ALQ111 全 NA 教训, 改 ALQ101 fallback
  ALQ101   = if ("ALQ101" %in% names(df)) na_codes(ALQ101, c(7, 9, 77, 99)) else NA,
  ALQ111   = if ("ALQ111" %in% names(df)) na_codes(ALQ111, c(7, 9)) else NA
)
df <- df %>% mutate(
  age = RIDAGEYR,
  age_group = factor(case_when(
    age >= 20 & age < 40 ~ "20-39",
    age >= 40 & age < 60 ~ "40-59",
    age >= 60            ~ ">=60"
  ), levels = c("20-39","40-59",">=60")),
  race = factor(recode(as.character(RIDRETH1),
                       "1"="Mexican American","2"="Other Hispanic",
                       "3"="Non-Hispanic White","4"="Non-Hispanic Black",
                       "5"="Other Race"),
                levels = c("Non-Hispanic White","Non-Hispanic Black",
                           "Mexican American","Other Hispanic","Other Race")),
  education = factor(case_when(
    DMDEDUC2 %in% c(1,2) ~ "Less than HS",
    DMDEDUC2 == 3        ~ "High school",
    DMDEDUC2 %in% c(4,5) ~ "College or above"
  ), levels = c("Less than HS","High school","College or above")),
  pir = INDFMPIR,
  bmi = BMXBMI,
  waist = BMXWAIST,
  smoke = factor(case_when(
    SMQ020 == 2 ~ "Never", SMQ020 == 1 ~ "Ever"
  ), levels = c("Never","Ever")),
  drink = factor(case_when(
    # 优先用 ALQ101 (跨周期全 D-J), ALQ111 fallback
    !is.na(ALQ101) & ALQ101 == 1 ~ "Yes",
    !is.na(ALQ101) & ALQ101 == 2 ~ "No",
    !is.na(ALQ111) & ALQ111 == 1 ~ "Yes",
    !is.na(ALQ111) & ALQ111 == 2 ~ "No"
  ), levels = c("No","Yes")),
  htn_med = as.integer(!is.na(BPQ020) & BPQ020 == 1),
  hypertension = htn_med,
  cotinine_log = ifelse(!is.na(LBXCOT) & LBXCOT > 0, log(LBXCOT), NA),
  smoke_objective = case_when(
    !is.na(LBXCOT) & LBXCOT >= 10  ~ "Active",
    !is.na(LBXCOT) & LBXCOT >= 0.05 ~ "Passive",
    !is.na(LBXCOT) & LBXCOT < 0.05 ~ "Non",
    TRUE ~ NA_character_
  )
)
if ("DR1TKCAL" %in% names(df)) df$kcal_day <- df$DR1TKCAL
if ("DR1TPROT" %in% names(df)) df$protein_g <- df$DR1TPROT

# Fish freq
fish_cols <- intersect(c("DRD340","DRD350","DRD360","DRD370"), names(df))
if (length(fish_cols) > 0) {
  for (c in fish_cols) df[[c]] <- na_codes(df[[c]], c(7,9,77,99,777,999))
  df$fish_freq_30d <- rowSums(df[, fish_cols, drop = FALSE], na.rm = TRUE)
}

# Postmenopausal
# Postmenopausal — RHQ060 是 age at menopause (continuous, NOT Y/N!)
# Sub-agent C 2026-05-23 实查 bug: 原用 RHQ060 %in% 1 全 NA / 0
# 修正: female + (RHQ060 valid age 10-80 OR age ≥ 50 default postmenopausal)
df$postmenopausal <- with(df, ifelse(
  sex_male == 0,
  ifelse(
    !is.na(RHQ060) & RHQ060 >= 10 & RHQ060 <= 80, 1L,    # 有明确 menopause age
    ifelse(age >= 50, 1L, 0L)                                # 否则 age ≥ 50 视为 postmenopausal
  ),
  NA_integer_                                                 # 男性 NA
))

# HSCRP inflammation 中介
df$hscrp_log <- ifelse(!is.na(df$LBXHSCRP) & df$LBXHSCRP > 0, log(df$LBXHSCRP), NA)

# ----------------------------------------------------------
# Step 13: Cohort split (Stack 1 主: PHTHTE D-J + fasting + 排除)
# ----------------------------------------------------------
cycle_years_main <- c(
  NHANES_2005_2006 = 2, NHANES_2007_2008 = 2,
  NHANES_2009_2010 = 2, NHANES_2011_2012 = 2,
  NHANES_2013_2014 = 2, NHANES_2015_2016 = 2,
  NHANES_2017_2018 = 2
)

df_main <- df %>% filter(cycle_tag %in% names(cycle_years_main))
log_flow("Stack 1 cohort (PHTHTE D-J + fasting + 排除)", nrow(df_main))

# Pooled weight — W16 Round 1 reset R-NHANES C1 fix:
#   Stack 1 = MEC ∩ PHTHTE ∩ Fasting. Per NCHS Series 2 No. 190 §4 LCD rule,
#   appropriate weight = most restrictive subsample = WTSAF2YR (fasting subsample).
#   pooled_saf_weight() applies WTSAF2YR × (cycle_years / total_years) cycle-by-cycle.
#   This replaces the previous pooled_mec_weight() which over-counted non-fasting.
if ("WTSAF2YR" %in% names(df_main)) {
  df_main$wt_pooled <- pooled_saf_weight(df_main, cycle_years_main)
  # Keep MEC pooled as descriptive backup (NOT used in design)
  df_main$wt_pooled_mec <- pooled_mec_weight(df_main, cycle_years_main)
  n_zero_saf <- sum(is.na(df_main$wt_pooled) | df_main$wt_pooled <= 0)
  if (n_zero_saf > 0) {
    cat(sprintf("[R-NHANES C1 fix] %d rows with WTSAF2YR ≤ 0 or NA — falling back to MEC weight\n", n_zero_saf))
    bad_idx <- is.na(df_main$wt_pooled) | df_main$wt_pooled <= 0
    df_main$wt_pooled[bad_idx] <- df_main$wt_pooled_mec[bad_idx]
  }
} else {
  cat("[WARN] WTSAF2YR not available → falling back to WTMEC2YR-based pooling (LCD violation)\n")
  df_main$wt_pooled <- pooled_mec_weight(df_main, cycle_years_main)
}
# W16 Round 3 SA-A3 cleanup (R-NHANES Round 2 M-7 fix):
# wt_diet_pooled was previously assigned here but never used downstream
# (design_main / design_mortality / design_t2d_prog / design_phth_pfas all
#  use wt_pooled = WTSAF2YR-pooled, not WTDRD1-pooled).
# Per Option A clean-removal recommendation: DROP the unused assignment.
# Limitation: kcal_day + protein_g enter as covariates in scripts/13_iptw.R
# (IPTW) and scripts/19_tripod_ai.R (TRIPOD-AI prediction); under the LCD
# rule these sub-cohorts (MEC ∩ PHTHTE ∩ Fasting ∩ DR1) would warrant a
# WTDRD1-pooled weight, but: (a) 13_iptw is a sensitivity probe not a
# primary inference and (b) 19_tripod_ai trains xgboost without survey
# weights (sample-conditional AUROC, per Methods §3.10). Documented in
# Methods §3.10 + §3.6 as a sensitivity/robustness limitation, not a
# primary-LCD violation. See _review_log/wave3_sa_a3_refs_scripts_csv_report.md.

# Core covariate strict
core_cov_strict <- c("age","race","education","pir","bmi","RIAGENDR","SDMVPSU","SDMVSTRA")
n_before <- nrow(df_main)
df_main <- df_main %>% filter(if_all(all_of(core_cov_strict), ~ !is.na(.)))
log_flow(sprintf("+ 核心协变量完整 (剔 %d)", n_before - nrow(df_main)), nrow(df_main))

# ----------------------------------------------------------
# Step 14: Stack 3 Mortality cohort
# ----------------------------------------------------------
df_mort <- df_main %>% filter(!is.na(ELIGSTAT), ELIGSTAT == 1)
log_flow("Stack 3 Mortality-linked cohort", nrow(df_mort))

# ----------------------------------------------------------
# Step 15: Stack 2 T2D progression (prediabetes subset)
# Prediabetes: HbA1c 5.7-6.4 OR fasting glucose 100-125
# ----------------------------------------------------------
df_t2d_prog <- df_main %>% filter(
  (LBXGH >= 5.7 & LBXGH < 6.5) | (LBXGLU >= 100 & LBXGLU < 126)
)
log_flow("Stack 2 Prediabetes subset (HbA1c 5.7-6.4 OR FBG 100-125)", nrow(df_t2d_prog))

# ----------------------------------------------------------
# Step 16: Stack 4 Phthalate+PFAS (与 007 PFAS subsample 交集)
# ----------------------------------------------------------
pfas_cols <- c("LBXPFOA","LBXPFOS","LBXPFNA","LBXPFHS","LBXPFDE","LBXMPAH")
pfas_avail <- intersect(pfas_cols, names(df_main))
if (length(pfas_avail) >= 4) {
  df_phth_pfas <- df_main %>% filter(rowSums(!is.na(across(all_of(pfas_avail)))) >= 4)
  log_flow("Stack 4 Phthalate+PFAS subset (PFAS 4+ 化合物 complete)", nrow(df_phth_pfas))
} else {
  df_phth_pfas <- data.frame()
}

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------
nhanes_final <- df_main

scale_attrs <- list(
  phth_cols = phth_cols_avail,
  pfas_cols_for_stack4 = pfas_avail,
  cycle_years_main = cycle_years_main
)

if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_final, df_mort, df_t2d_prog, df_phth_pfas,
     cycle_years_main, mort_all, rx_all, scale_attrs,
     file = "data/processed/nhanes_final.RData")

cat("\n========================================\n")
cat(sprintf("Stack 1 主分析 N = %d\n", nrow(nhanes_final)))
cat(sprintf("Stack 2 T2D progression 子集 N = %d\n", nrow(df_t2d_prog)))
cat(sprintf("Stack 3 Mortality eligible N = %d / 全因死亡 = %d / CM 死亡 = %d / DM 死亡 = %d\n",
            nrow(df_mort),
            sum(df_mort$mort_allcause, na.rm = TRUE),
            sum(df_mort$mort_cm, na.rm = TRUE),
            sum(df_mort$mort_dm, na.rm = TRUE)))
cat(sprintf("Stack 4 Phthalate+PFAS subset N = %d\n", nrow(df_phth_pfas)))
cat(sprintf("\nIR primary outcome:\n"))
cat(sprintf("  HOMA-IR median: %.2f\n", median(nhanes_final$homa_ir, na.rm = TRUE)))
cat(sprintf("  IR binary (HOMA>=2.5): %d / %d (%.1f%%)\n",
            sum(nhanes_final$ir_binary == 1, na.rm = TRUE), nrow(nhanes_final),
            100 * mean(nhanes_final$ir_binary == 1, na.rm = TRUE)))
cat(sprintf("  IR binary Asian (HOMA>=3.6): %d / %d (%.1f%%)\n",
            sum(nhanes_final$ir_binary_asian == 1, na.rm = TRUE), nrow(nhanes_final),
            100 * mean(nhanes_final$ir_binary_asian == 1, na.rm = TRUE)))
cat(sprintf("========================================\n"))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
flow_df <- do.call(rbind, flow)
write.csv(flow_df, "output/tables/flow_counts.csv", row.names = FALSE)
cat("已保存 data/processed/nhanes_final.RData + output/tables/flow_counts.csv\n")
