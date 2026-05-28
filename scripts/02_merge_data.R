# ============================================
# 009 / scripts/02_merge_data.R (复用 007 + 加 PHTHTE)
# 合并 NHANES PFAS C-J + P_PFAS + PHTHTE D-J + 26 模块 + NCHS LMF
# 输入:  data/raw/*.xpt (复用 007 raw) + data/raw/mortality/*.dat
# 输出:  data/processed/nhanes_raw_merged.RData
#         (含 nhanes_all + rx_all + mort_all + cycle_tag)
#
# 关键差异 vs 007:
#   - 加 PHTHTE 模块 (D-J 7 cycles; P_PHTHTE 未发布)
#   - PHTHTE 文件名一致 (PHTHTE_D 到 PHTHTE_J, 无漂移), 用 modules_uniform 自动 merge
#   - PFAS 仍合并 (作 Stack 4 cross-program synergy)
# ============================================

library(haven); library(dplyr); library(purrr)

cat("========================================\n")
cat("009 / 合并 PFAS C-J + P_ + PHTHTE D-J + 26 模块 + 死亡链接\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"

# --------------------------------------------------
# Step 1: 读所有 .xpt
# --------------------------------------------------
xpt_files <- list.files(raw_dir, pattern = "\\.xpt$", full.names = TRUE,
                        ignore.case = TRUE)
cat(sprintf("找到 %d 个 .xpt 文件\n\n", length(xpt_files)))

read_safe <- function(p) {
  tryCatch(read_xpt(p),
           error = function(e) {
             cat(sprintf("  [读失败] %s : %s\n", basename(p), e$message))
             NULL
           })
}

data_list <- map(xpt_files, read_safe)
names(data_list) <- gsub("\\.xpt$", "", basename(xpt_files), ignore.case = TRUE)
data_list <- data_list[!sapply(data_list, is.null)]
cat(sprintf("成功读入 %d 个数据帧\n\n", length(data_list)))

# --------------------------------------------------
# Step 2: 周期 + 后缀 + PFAS 文件名漂移
# --------------------------------------------------
cycles <- data.frame(
  cycle_tag = c("NHANES_2003_2004","NHANES_2005_2006","NHANES_2007_2008",
                "NHANES_2009_2010","NHANES_2011_2012","NHANES_2013_2014",
                "NHANES_2015_2016","NHANES_2017_2018","PrePandemic_2017_March2020"),
  suffix = c("_C","_D","_E","_F","_G","_H","_I","_J","P_"),
  year   = c(2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2017),
  stringsAsFactors = FALSE
)

# PFAS 跨周期文件名映射 (CDC 命名漂移)
pfas_map <- data.frame(
  suffix = c("_C","_D","_E","_F","_G","_H","_I","_J","P_"),
  file   = c("L24PFC_C","PFC_D","PFC_E","PFC_F","PFC_G",
             "PFAS_H","PFAS_I","PFAS_J","P_PFAS"),
  stringsAsFactors = FALSE
)

# 跨周期"宽"模块 (每人 1 行)
# 009 加 PHTHTE (D-J 7 cycles; 文件名一致 PHTHTE_D...PHTHTE_J 无漂移)
modules_uniform <- c(
  "DEMO","BMX","BPX","SMQ","ALQ","DIQ","BPQ","MCQ","PAQ","RHQ",
  "PBCD","IHGEM","UHM",
  "LUX","BIOPRO","CBC","HSCRP",
  "GHB","GLU","INS","HDL","TCHOL","TRIGLY",
  "DR1TOT","DR2TOT","VID","ALB_CR","COT",
  "HEPB_S","HEPC",
  "PHTHTE"  # 009: Phthalate metabolites (D-J 7 cycles, P_PHTHTE 未发布)
)
modules_long <- c("RXQ_RX")

# --------------------------------------------------
# Step 3: safe_left_join 工具
# --------------------------------------------------
safe_left_join <- function(x, y, key = "SEQN") {
  if (is.null(y) || nrow(y) == 0) return(x)
  if (!key %in% names(y)) return(x)
  dup <- intersect(setdiff(names(x), key), names(y))
  if (length(dup) > 0) y <- y[, setdiff(names(y), dup), drop = FALSE]
  left_join(x, y, by = key)
}

# --------------------------------------------------
# Step 4: 单周期合并
# --------------------------------------------------
merge_cycle <- function(cyc_row) {
  suffix <- cyc_row$suffix
  year <- cyc_row$year
  tag <- cyc_row$cycle_tag
  cat(sprintf("--- 周期 %s (%s) ---\n", tag, suffix))

  # DEMO 基表 (P_ 用 P_DEMO; std 用 DEMO_X)
  base_name <- if (suffix == "P_") "P_DEMO" else paste0("DEMO", suffix)
  if (!base_name %in% names(data_list)) {
    cat(sprintf("  [警告] 基表 %s 不存在, 跳过\n\n", base_name))
    return(NULL)
  }
  result <- data_list[[base_name]]
  result$cycle_tag <- tag
  result$cycle_year <- year
  result$is_prepandemic <- (suffix == "P_")

  # 合并 modules_uniform (skip DEMO 自身)
  for (mod in modules_uniform[-1]) {
    key <- if (suffix == "P_") paste0("P_", mod) else paste0(mod, suffix)
    if (key %in% names(data_list)) {
      result <- safe_left_join(result, data_list[[key]])
    }
  }

  # PFAS (文件名漂移): 用 pfas_map 找文件
  pfas_file <- pfas_map$file[pfas_map$suffix == suffix]
  if (length(pfas_file) > 0 && pfas_file %in% names(data_list)) {
    result <- safe_left_join(result, data_list[[pfas_file]])
    cat(sprintf("  + PFAS %s 已合并\n", pfas_file))
  } else {
    cat(sprintf("  [警告] PFAS 文件 %s 缺失\n", pfas_file))
  }

  cat(sprintf("  -> %d 行 x %d 列\n\n", nrow(result), ncol(result)))
  result
}

merged_list <- lapply(seq_len(nrow(cycles)), function(i) merge_cycle(cycles[i, ]))
merged_list <- merged_list[!sapply(merged_list, is.null)]
nhanes_all <- bind_rows(merged_list)

cat(sprintf("========================================\n"))
cat(sprintf("9 周期合并完成: %d 行 x %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))
cat(sprintf("========================================\n\n"))

# --------------------------------------------------
# Step 4.5: 长表 RXQ_RX (用药)
# --------------------------------------------------
cat("--- Step 4.5: 长表 RXQ_RX 单独合并 ---\n")
rx_list <- list()
for (suf in cycles$suffix) {
  key <- if (suf == "P_") "P_RXQ_RX" else paste0("RXQ_RX", suf)
  if (key %in% names(data_list)) rx_list[[key]] <- data_list[[key]]
}
rx_all <- if (length(rx_list) > 0) bind_rows(rx_list) else data.frame(SEQN = integer(0))
cat(sprintf("rx_all: %d 行 x %d 列\n\n", nrow(rx_all), ncol(rx_all)))

# --------------------------------------------------
# Step 5: 读 NCHS LMF (PUBLICID 1-14 含 SEQN 1-6, fixed-width fwf)
# --------------------------------------------------
cat("--- Step 5: 接入死亡链接 ---\n")
mort_files <- list.files(mort_dir, pattern = "\\.dat$", full.names = TRUE)
cat(sprintf("找到 %d 个 LMF 文件\n", length(mort_files)))

read_mort <- function(path) {
  widths <- c(
    SEQN = 6, PADDING1 = 8, ELIGSTAT = 1, MORTSTAT = 1,
    UCOD_LEADING = 3, DIABETES = 1, HYPERTEN = 1,
    DODQTR = 2, DODYEAR = 4, WGT_NEW = 8, SA_WGT_NEW = 8,
    PERMTH_INT = 3, PERMTH_EXM = 3
  )
  df <- tryCatch(
    read.fwf(path, widths = widths, header = FALSE,
             na.strings = c("", "."), stringsAsFactors = FALSE,
             col.names = names(widths)),
    error = function(e) { cat(sprintf("  [读失败] %s : %s\n", basename(path), e$message)); NULL }
  )
  if (is.null(df)) return(NULL)
  df$PADDING1 <- NULL
  df$SEQN <- as.integer(df$SEQN)
  for (col in c("ELIGSTAT","MORTSTAT","UCOD_LEADING","DIABETES","HYPERTEN",
                "PERMTH_INT","PERMTH_EXM")) {
    df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  }
  df
}
mort_all <- bind_rows(map(mort_files, read_mort))
cat(sprintf("死亡链接合并: %d 行\n", nrow(mort_all)))
cat(sprintf("  ELIGSTAT==1 (合格): %d\n", sum(mort_all$ELIGSTAT == 1, na.rm=TRUE)))
cat(sprintf("  MORTSTAT==1 (死亡): %d\n", sum(mort_all$MORTSTAT == 1, na.rm=TRUE)))

nhanes_all <- safe_left_join(nhanes_all, mort_all)
cat(sprintf("\n合并死亡后: %d 行 x %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))

# --------------------------------------------------
# 保存
# --------------------------------------------------
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_all, rx_all, mort_all, cycles,
     file = "data/processed/nhanes_raw_merged.RData")
cat("\n已保存 data/processed/nhanes_raw_merged.RData\n")
cat(sprintf("  含 nhanes_all (n=%d), rx_all (n=%d), mort_all (n=%d), cycles (9 周期 spec)\n",
            nrow(nhanes_all), nrow(rx_all), nrow(mort_all)))
cat("========================================\n")
