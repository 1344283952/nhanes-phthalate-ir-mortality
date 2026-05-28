# ============================================
# _make_delivery_pi.R (009 phthalate-IR-mortality)
# 拆分打两个独立 zip 给通讯作者：
#
#   1. 投稿主文档包.zip (~5MB)
#      仅含 投稿系统上传的文档 + 给通讯作者看的操作手册
#
#   2. GitHub上传包.zip (~50MB)
#      仅含 GitHub 仓库内容（scripts + data/processed + output + README.md）
#      通讯作者解压 → git init → push 即可
#
# 第 3 个完整包 (GitHub 上传包) 由 _make_delivery_v2.R 生成
# ============================================

cat("========================================\n")
cat("打包通讯作者友好分包 (主投 Diabetes Care)\n")
cat("========================================\n\n")

ts <- format(Sys.time(), "%Y%m%d_%H%M")

if (!requireNamespace("zip", quietly = TRUE)) {
  install.packages("zip", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# ============================================
# Pack 1: 投稿主文档包
# ============================================

pack1_name <- "投稿主文档包_009.zip"
top1 <- "投稿主文档包_009"

stage1 <- file.path(tempdir(), top1)
unlink(stage1, recursive = TRUE)
dir.create(stage1, recursive = TRUE)

cat(">>> 打包 1: 投稿主文档包\n\n")

# ============================================
# Pre-build: pandoc 自动渲染 .md → .docx (期刊投稿接受 .docx)
# 2026-05-18 W22 fix: 之前模板未自动渲染, PI 收到 .md 无法直接上传 EM 投稿系统
# 现在打包前先检查 .docx 是否存在/过期, 缺则用 pandoc 自动生成
# ============================================
pandoc_available <- nchar(Sys.which("pandoc")) > 0
if (!pandoc_available) {
  warning("pandoc 不在 PATH! 跳过 .md → .docx 自动渲染. PI 拿到 .md 后需手动转 docx.")
} else {
  # 项目根 .md 文件 (排除内部 _* / task / review log / archive)
  all_md <- list.files(pattern = "\\.md$", recursive = FALSE)
  exclude_pat <- c("^task\\.md$", "^_review_log", "^_BKMR_DAILY", "^_QA_pipeline", "^_m13_results_summary",
                   "^继续工作", "^_dag_spec", "^_figures_spec", "^_rebuttal",
                   "^manuscript_v0", "_pre_W2\\d", "_archive", "^_drafts")
  render_md <- all_md[!Reduce(`|`, lapply(exclude_pat, function(p) grepl(p, all_md)))]
  cat(sprintf("Pre-build: 检查 %d 个 .md 是否需渲染 .docx\n", length(render_md)))
  for (md_file in render_md) {
    docx_file <- sub("\\.md$", ".docx", md_file)
    needs_render <- !file.exists(docx_file) ||
                    (file.mtime(md_file) > file.mtime(docx_file))
    if (needs_render) {
      tryCatch({
        rc <- system2("pandoc",
                       args = c(shQuote(md_file), "-o", shQuote(docx_file),
                                "--from=gfm+yaml_metadata_block", "--to=docx"),
                       stdout = FALSE, stderr = FALSE)
        if (file.exists(docx_file)) cat(sprintf("  rendered %s\n", docx_file))
      }, error = function(e) {
        warning(sprintf("pandoc render %s failed: %s", md_file, e$message))
      })
    }
  }
  cat("\n")
}

# 顶层文档 — 极简版（7 个核心，全是投稿系统会用到的 + 操作手册）
# 注：references.bib 不打进来 —— Vancouver 已转完，PI 不需要 .bib
#
# W9-F (2026-05-23) fix: source .md 文件名 ≠ 期刊投稿期望的短名。
# 例如 source 是 manuscript_v1.md → render manuscript_v1.docx，
# 但投稿系统 / PI 期望 manuscript.docx (memory: feedback_manuscript_naming_for_submission).
# 用 rename_map 显式映射 source → delivered name，避免 silent skip。
rename_map <- list(
  "投稿操作指南.docx"             = "投稿操作指南.docx",          # ⭐ 第一份打开的：4 步走
  "manuscript_v1.docx"            = "manuscript.docx",            # 主文档
  "cover_letter_DiabetesCare.docx" = "cover_letter.docx",         # 投稿信
  "supplementary_information.docx" = "supplementary.docx",        # 补充材料
  "STROBE_checklist.docx"         = "STROBE_checklist.docx",      # STROBE 22 项（期刊要求附）
  "STROBE_MR_checklist.docx"      = "STROBE_MR_checklist.docx",   # STROBE-MR (NA per N/A declaration; required attachment)
  "AGReMA_checklist.docx"         = "AGReMA_checklist.docx",      # AGReMA mediation reporting
  "TRIPOD_AI_checklist.docx"      = "TRIPOD_AI_checklist.docx",   # TRIPOD-AI 27/27 prediction
  "SAGER_checklist.docx"          = "SAGER_checklist.docx",       # SAGER sex/gender reporting
  "suggested_reviewers.docx"      = "suggested_reviewers.docx",   # 5 位推荐审稿人
  "_osf_preregistration.docx"     = "OSF_preregistration.docx"    # OSF 预注册（粘贴到 OSF 网站）
)

for (src in names(rename_map)) {
  dst <- rename_map[[src]]
  if (file.exists(src)) {
    file.copy(src, file.path(stage1, dst), overwrite = TRUE)
    if (src == dst) {
      cat(sprintf("  + %s\n", src))
    } else {
      cat(sprintf("  + %s -> %s\n", src, dst))
    }
  } else {
    cat(sprintf("  - 缺 source: %s\n", src))
  }
}

# Figures (4 main + graphical abstract + 9 supplementary = 14 files PI uploads to ScholarOne)
fig_dir1 <- file.path(stage1, "figures")
dir.create(fig_dir1, showWarnings = FALSE)
# Main 4 + Graphical Abstract (DC 2026 required at submission)
main_figs <- c(
  "fig1_consort.png",
  "fig5_cmaverse_4way.png",
  "fig7_multiverse.png",
  "fig12_nomogram.png",
  "graphical_abstract.png"    # PNG only per pi_zip_single_figure_format
)
for (f in main_figs) {
  sp <- file.path("output", "figures", f)
  if (file.exists(sp)) {
    file.copy(sp, file.path(fig_dir1, f), overwrite = TRUE)
  }
}
# Supplementary Figures S1-S9 (per manuscript §7 recap + ScholarOne separate uploads)
# Rename to supp_sN_* prefix for clear PI ordering in EM system
supp_figs <- list(
  "supp_s1_dag.png"             = "fig2_dag.png",            # S1 DAG (ggdag)
  "supp_s2_rcs_dose_response.png" = "fig3_rcs_homa.png",     # S2 RCS dose-response (HOMA-IR primary)
  "supp_s3_wqs_weights.png"     = "fig4_wqs_weights.png",    # S3 WQS positive-direction weights
  "supp_s4_iptw_love.png"       = "fig6_iptw_balance.png",   # S4 IPTW Love plot covariate balance
  "supp_s5_subgroup_forest.png" = "fig8_subgroup_forest.png",# S5 9-subgroup × interaction forest
  "supp_s6_evalue.png"          = "fig9_evalue.png",         # S6 E-value point + CI bound
  "supp_s7_bayesian_posterior.png" = "fig10_bayesian_gcomp.png",  # S7 Bayesian g-comp posterior
  "supp_s8_dca_calibration.png" = "fig11_dca_calibration.png",# S8 DCA + calibration belt stitched
  "supp_s9_bkmr_pip.png"        = "fig13_bkmr_pip.png"       # S9 BKMR posterior inclusion probabilities (post-convergence real data)
)
for (dst_name in names(supp_figs)) {
  sp <- file.path("output", "figures", supp_figs[[dst_name]])
  if (file.exists(sp)) {
    file.copy(sp, file.path(fig_dir1, dst_name), overwrite = TRUE)
  } else {
    cat(sprintf("  ⚠️ supp source missing: %s\n", sp))
  }
}
cat(sprintf("  + figures/ (%d 张: 4 main + 1 graphical abstract + 9 supplementary S1-S9)\n", length(list.files(fig_dir1))))

# Tables (主文 7 张 XLSX)
tab_dir1 <- file.path(stage1, "tables")
dir.create(tab_dir1, showWarnings = FALSE)
main_tabs <- list.files("output/tables", pattern = "table[1-9].*\\.xlsx$", full.names = TRUE)
file.copy(main_tabs, tab_dir1, overwrite = TRUE)
cat(sprintf("  + tables/ (%d 张主表 xlsx)\n", length(main_tabs)))

# README.txt
readme1 <- c(
  "# 投稿主文档包",
  "",
  sprintf("打包时间: %s", format(Sys.time())),
  "目标期刊: Diabetes Care (American Diabetes Association, IF 16)",
  "",
  "## 文件清单",
  "- 投稿操作指南.docx       ⭐ 第一份打开",
  "- manuscript.docx          主稿（投稿系统上传）",
  "- cover_letter.docx        投稿信（投稿系统上传）",
  "- supplementary.docx       补充材料（投稿系统上传）",
  "- STROBE_checklist.docx    STROBE 22 项",
  "- suggested_reviewers.docx 5 位推荐审稿人（填表照抄）",
  "- OSF_preregistration.docx OSF 注册时粘贴",
  "- figures/                 4 张主图 (CONSORT + CMAverse 4-way + Multiverse + Nomogram)",
  "- tables/                  主表 xlsx",
  "",
  "## 第一步",
  "打开 投稿操作指南.docx，按 4 步走即可。"
)
writeLines(readme1, file.path(stage1, "README.txt"))
cat("  + README.txt\n")

cat("\n压缩 pack 1 ...\n")
old_wd <- getwd()
setwd(dirname(stage1))
zip::zip(zipfile = file.path(old_wd, pack1_name),
         files   = basename(stage1),
         recurse = TRUE,
         mode    = "cherry-pick")
setwd(old_wd)
unlink(stage1, recursive = TRUE)

zsize1 <- file.size(pack1_name)
cat(sprintf("\n[OK] pack 1: %s (%.1f MB)\n", pack1_name, zsize1 / 1024 / 1024))

# ============================================
# Pack 2: GitHub 上传包
# ============================================

pack2_name <- "GitHub上传包_009.zip"
top2 <- "GitHub上传包_009"

stage2 <- file.path(tempdir(), top2)
unlink(stage2, recursive = TRUE)
dir.create(stage2, recursive = TRUE)

cat("\n\n>>> 打包 2: GitHub 上传包\n\n")

# README.md (GitHub 仓库首页) — 唯一的顶层 doc
if (file.exists("README_github.md")) {
  file.copy("README_github.md", file.path(stage2, "README.md"), overwrite = TRUE)
  cat("  + README.md (GitHub 仓库首页)\n")
}
# 注：不包含 task.md (内部规划档) 和 references.bib (写稿用，复现无关)

# .gitignore
gitignore <- c(
  "# raw NHANES data (~700 MB, public domain — re-download via script 01)",
  "data/raw/",
  "",
  "# R / IDE",
  ".Rhistory",
  ".RData",
  ".Rproj.user/",
  "*.Rproj",
  "",
  "# OS",
  "Thumbs.db",
  ".DS_Store",
  "",
  "# build outputs",
  "*.docx",
  "*.zip"
)
writeLines(gitignore, file.path(stage2, ".gitignore"))
cat("  + .gitignore\n")

# LICENSE (MIT)
license <- c(
  "MIT License",
  "",
  "Copyright (c) 2026 Jie Li, Xiubo Sun, Jing Zhang, Lijie Zhai, Ling Yu / The Second Hospital of Jilin University",
  "",
  "Permission is hereby granted, free of charge, to any person obtaining a copy",
  "of this software and associated documentation files (the \"Software\"), to deal",
  "in the Software without restriction, including without limitation the rights",
  "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell",
  "copies of the Software, and to permit persons to whom the Software is",
  "furnished to do so, subject to the following conditions:",
  "",
  "The above copyright notice and this permission notice shall be included in",
  "all copies or substantial portions of the Software.",
  "",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND."
)
writeLines(license, file.path(stage2, "LICENSE"))
cat("  + LICENSE (MIT)\n")

# scripts/ — 论文复现需要的核心脚本 (analysis + helpers); 剔除打包/探索/QA 工具
sc_dir2 <- file.path(stage2, "scripts")
dir.create(sc_dir2, showWarnings = FALSE)
reproducibility_scripts <- c(
  "00_install_packages.R", "01_download_data.R", "02_merge_data.R",
  "03_clean_data.R", "04_survey_design.R", "05_table1.R",
  "06_bkmr_phth_homa_ir_checkpoint.R",
  "07_cox_mortality.R", "08_logistic_ir.R", "09_rcs.R",
  "10_qgcomp.R", "11_wqs.R", "12_cmaverse_4way.R", "13_iptw.R",
  "14_evalue.R", "15_dag.R", "16_subgroup.R", "17_sensitivity.R",
  "18_evalue_ci_bound.R", "19_tripod_ai.R", "20_nri_idi.R",
  "21_dca.R", "22_calibration.R", "23_nomogram.R",
  "24_multiverse.R", "25_bayesian_gcomp.R", "26_episensr.R",
  "27_lag_analysis.R", "28_ctat_methods_table.R",
  "29_fine_gray.R", "30_reri_ap_si.R", "31_multiplicity_padjust.R",
  "32_retrodesign.R", "33_figures_rework.R",
  "run_all.R", "_extract_numbers.R"
)
for (sc in reproducibility_scripts) {
  sp <- file.path("scripts", sc)
  if (file.exists(sp)) file.copy(sp, file.path(sc_dir2, sc), overwrite = TRUE)
}
cat(sprintf("  + scripts/ (%d R files — 仅论文复现需要的)\n",
            length(list.files(sc_dir2))))

# data/processed/
proc_dir2 <- file.path(stage2, "data", "processed")
dir.create(proc_dir2, showWarnings = FALSE, recursive = TRUE)
if (dir.exists("data/processed")) {
  rd <- list.files("data/processed", pattern = "\\.(RData|rds|csv)$", full.names = TRUE)
  # Exclude internal _* files + BKMR checkpoint .rds + BKMR results.RData (>100MB GitHub per-file hard limit)
  rd <- rd[!grepl("^_|/_", basename(rd))]
  rd <- rd[!grepl("checkpoint\\.rds$", basename(rd), ignore.case = TRUE)]
  rd <- rd[!grepl("bkmr.*results\\.RData$", basename(rd), ignore.case = TRUE)]  # 311 MB — users re-run 06_bkmr*.R (~10h) per feedback_github_zip_no_bkmr_rds
  for (f in rd) file.copy(f, file.path(proc_dir2, basename(f)), overwrite = TRUE)
  cat(sprintf("  + data/processed/ (%d files; BKMR checkpoint .rds + results.RData excluded — re-run scripts/06_bkmr*.R to regenerate)\n", length(rd)))
}
# Also filter output/tables for BKMR results.RData (>100MB hard limit)
# (handled inside copy_dir_clean below via additional pattern)
# 也建个 data/raw/ 空目录加 .gitkeep
raw_dir2 <- file.path(stage2, "data", "raw")
dir.create(raw_dir2, showWarnings = FALSE, recursive = TRUE)
writeLines(c(
  "# data/raw/ — NHANES XPT files",
  "",
  "Run `Rscript scripts/01_download_data.R` to populate this folder."
), file.path(raw_dir2, "README.md"))

# output/
copy_dir_clean <- function(src, dst_parent) {
  dst <- file.path(dst_parent, basename(src))
  dir.create(dst, showWarnings = FALSE, recursive = TRUE)
  items <- list.files(src, all.files = FALSE, no.. = TRUE, include.dirs = TRUE)
  items <- items[!grepl("^[._]", items)]  # exclude .-prefix (.git/.DS_Store) + _-prefix (_ai_tell_report.csv / _ceiling_report.csv / _consistency_report.csv / _review_log/ etc.) internal QA files from GitHub pack-2 zip
  items <- items[!grepl("\\.(tiff?|pdf)$", items, ignore.case = TRUE)]  # exclude .tiff and .pdf figures from GitHub zip (PNG variants are sufficient for reproducibility; TIFF retained in working directory for journal-quality submission)
  items <- items[!grepl("checkpoint\\.rds$", items, ignore.case = TRUE)]  # exclude BKMR checkpoint .rds (~310 MB) — exceeds GitHub 100MB per-file limit; users can re-run BKMR from scratch (10h) or request from authors
  items <- items[!grepl("bkmr.*results\\.RData$", items, ignore.case = TRUE)]  # exclude BKMR results.RData (311 MB) — exceeds GitHub 100MB; users re-run scripts/06_bkmr*.R (~10h)
  for (item in items) {
    sp <- file.path(src, item)
    if (dir.exists(sp)) copy_dir_clean(sp, dst)
    else file.copy(sp, file.path(dst, item), overwrite = TRUE)
  }
}
if (dir.exists("output")) {
  copy_dir_clean("output", stage2)
  n_tab <- length(list.files(file.path(stage2, "output", "tables")))
  n_fig <- length(list.files(file.path(stage2, "output", "figures")))
  cat(sprintf("  + output/tables/ (%d files)\n", n_tab))
  cat(sprintf("  + output/figures/ (%d files)\n", n_fig))
}

cat("\n压缩 pack 2 ...\n")
setwd(dirname(stage2))
zip::zip(zipfile = file.path(old_wd, pack2_name),
         files   = basename(stage2),
         recurse = TRUE,
         mode    = "cherry-pick")
setwd(old_wd)
unlink(stage2, recursive = TRUE)

zsize2 <- file.size(pack2_name)
cat(sprintf("\n[OK] pack 2: %s (%.1f MB)\n", pack2_name, zsize2 / 1024 / 1024))

# ============================================
# 总结
# ============================================

cat("\n========================================\n")
cat("[完成] 通讯作者友好分包齐\n")
cat("========================================\n")
cat(sprintf("  📄 %s  (%.1f MB)  ← 投稿系统上传 + 通讯作者看\n",
            pack1_name, zsize1 / 1024 / 1024))
cat(sprintf("  🐙 %s  (%.1f MB)  ← 解压后 git push\n",
            pack2_name, zsize2 / 1024 / 1024))
cat("\n通讯作者先看 投稿主文档包/投稿操作指南.docx\n")
cat("========================================\n")
