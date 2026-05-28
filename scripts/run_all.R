# ============================================
# run_all.R
# 一键执行全部分析脚本
# 用法: 在 projects/<NNN_xxx>/ 目录下执行
#   cd projects/<NNN_xxx>
#   Rscript scripts/run_all.R
# 所有相对路径 (scripts/, data/, output/) 以当前 cwd 为基准
# ============================================

cat("╔══════════════════════════════════════════╗\n")
cat("║  NHANES 论文复现 - 全流程执行            ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

start_time <- Sys.time()

scripts <- c(
  "scripts/01_download_data.R",
  "scripts/02_merge_data.R",
  "scripts/03_clean_data.R",
  "scripts/04_survey_design.R",
  "scripts/05_table1.R",
  "scripts/06_table2.R",
  "scripts/07_table3.R",
  "scripts/08_table4.R",
  "scripts/09_figures.R"
)

for (s in scripts) {
  cat(paste0("\n>>> 正在执行: ", s, "\n"))
  cat(paste0(rep("-", 50), collapse = ""), "\n")
  
  tryCatch({
    source(s, echo = FALSE)
    cat(paste0(">>> ", s, " 执行成功 ✓\n"))
  }, error = function(e) {
    cat(paste0(">>> ", s, " 执行失败 ✗\n"))
    cat(paste0("    错误: ", e$message, "\n"))
  })
}

end_time <- Sys.time()
elapsed <- round(difftime(end_time, start_time, units = "mins"), 1)

cat(paste0("\n╔══════════════════════════════════════════╗\n"))
cat(paste0("║  全部完成! 总耗时: ", elapsed, " 分钟\n"))
cat(paste0("╚══════════════════════════════════════════╝\n"))

cat("\n输出文件:\n")
cat("  output/tables/table1.csv\n")
cat("  output/tables/table2.csv\n")
cat("  output/tables/table3.csv\n")
cat("  output/tables/table4.csv\n")
cat("  output/figures/figure2.png\n")
