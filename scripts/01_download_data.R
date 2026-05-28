# ============================================
# 009 / scripts/01_download_data.R
# 下载 NHANES PHTHTE D-J (复用 007 已下) + INS + GLU 协变量 + NCHS LMF
#
# 论文: Phthalate × Insulin Resistance + Cardiometabolic Mortality
#       Plan M13-α 主投 Diabetes Care IF 16
# Cohort:
#   Stack 1 主: PHTHTE D-J ∩ fasting (PHAFSTHR ≥ 8.5) ∩ ≥20, N ≈ 1,974 (W10 quick 实查)
#   Stack 3 mortality: + NCHS LMF
#   Stack 4 联合: + PFAS (复用 007 cohort 交集 ~800)
#
# 关键: 007 raw 已含 PHTHTE_D-J + INS_D-J + GLU_D-J + PFAS + 协变量
#       009 复用 007 raw, P_PHTHTE 未发布 (007 W2 confirm 2026-05-19)
# ============================================

cat("========================================\n")
cat("009 / NHANES PHTHTE + INS + GLU + 协变量下载 (Plan M13-α)\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"
log_path <- "data/raw/_download.log"
if (!dir.exists(raw_dir))  dir.create(raw_dir,  recursive = TRUE)
if (!dir.exists(mort_dir)) dir.create(mort_dir, recursive = TRUE)

# --- 复用 007 已下文件 ---
shared_dir <- "../../projects/007_pfas_masld_mortality/data/raw"
if (dir.exists(shared_dir)) {
  cat("[REUSE] 复用 007 已下文件 (含 PHTHTE D-J + INS + GLU + mortality)...\n")
  shared_files <- list.files(shared_dir, pattern = "\\.xpt$", full.names = TRUE)
  copied <- 0
  for (f in shared_files) {
    fn <- basename(f)
    dest <- file.path(raw_dir, fn)
    if (!file.exists(dest)) { file.copy(f, dest, overwrite = FALSE); copied <- copied + 1 }
  }
  cat(sprintf("[REUSE] 已复用 %d / %d 个 007 xpt (%d 已存在跳过)\n",
              copied, length(shared_files), length(shared_files) - copied))

  shared_mort <- file.path(shared_dir, "mortality")
  if (dir.exists(shared_mort)) {
    mort_files <- list.files(shared_mort, pattern = "\\.dat$", full.names = TRUE)
    copied_m <- 0
    for (f in mort_files) {
      fn <- basename(f)
      dest <- file.path(mort_dir, fn)
      if (!file.exists(dest)) { file.copy(f, dest, overwrite = FALSE); copied_m <- copied_m + 1 }
    }
    cat(sprintf("[REUSE] 已复用 %d / %d 个 mortality .dat\n",
                copied_m, length(mort_files)))
  }
} else {
  cat("[INFO] 007 raw 目录不存在, 全量下载\n")
}

writeLines(c("# 009 NHANES download log", paste("# started:", Sys.time())), log_path)

n_xpt <- length(list.files(raw_dir, pattern = "\\.xpt$"))
n_mort <- length(list.files(mort_dir, pattern = "\\.dat$"))
cat(sprintf("\ndata/raw/         %d 个 .xpt\n", n_xpt))
cat(sprintf("data/raw/mortality/ %d 个 .dat\n", n_mort))
cat("\n009 下载完成 (复用 007 raw)\n")
cat("========================================\n")
