# ============================================
# 00_install_packages.R
# 首次运行：安装所有依赖的 R 包
# 只需运行一次
# ============================================

cat("========================================\n")
cat("正在安装 R 包，首次运行可能需要 5-10 分钟...\n")
cat("========================================\n\n")

# 定义需要安装的包
packages <- c(
  # Core data + reporting (Table 1-4 + figures)
  "tidyverse",    # dplyr/tidyr/ggplot2/purrr/stringr
  "survey",       # 复杂抽样设计与加权分析 (svydesign + svyglm + svycoxph)
  "haven",        # 读取 .xpt (SAS transport) 文件
  "broom",        # 回归模型 tidy 输出
  "tableone",     # 基线特征表
  "openxlsx",     # Excel 输出
  "DiagrammeR",   # 流程图 (Figure 1 CONSORT)
  "DiagrammeRsvg",# SVG export for DiagrammeR
  "rsvg",         # SVG → PNG/PDF
  "magick",       # 图像 LZW TIFF conversion
  "corrplot",     # 相关矩阵
  "gridExtra",    # multi-panel ggplot composition
  "cowplot",      # ggplot themes + plot_grid (graphical abstract)
  # Mixture estimators (Scripts 06, 10, 11)
  "bkmr",         # Bayesian Kernel Machine Regression (Bobb 2015)
  "qgcomp",       # quantile g-computation (Keil 2020)
  "gWQS",         # weighted quantile sum (Carrico 2015)
  # Causal mediation + IPTW (Scripts 12, 13)
  "CMAverse",     # 4-way decomposition (VanderWeele 2014; Shi 2021); install from GitHub: remotes::install_github("BS1125/CMAverse")
  "WeightIt",     # propensity weighting (Greifer 2024)
  "cobalt",       # covariate balance diagnostics
  "EValue",       # E-value (VanderWeele-Ding 2017)
  # Survival + prediction (Scripts 07, 19-23, 29)
  "survival",     # Cox + Schoenfeld
  "rms",          # nomogram + cph + lrm
  "Hmisc",        # rcs + ggplot helpers
  "cmprsk",       # Fine-Gray competing-risk
  "nricens",      # NRI + IDI
  "PredictABEL",  # discrimination/calibration
  "pROC",         # ROC + AUC
  "dcurves",      # decision-curve analysis
  "givitiR",      # calibration belt (Nattino 2014)
  "ResourceSelection", # Hosmer-Lemeshow test
  "xgboost",      # gradient boosting (Tier 4 TRIPOD-AI)
  # MICE imputation (Script 17 sensitivity S4)
  "mice",         # multiple imputation
  # Bayesian g-computation (Script 25)
  "rstanarm",     # Bayesian rstan front-end
  "rstan",        # backing
  "posterior",    # posterior summaries
  "bayesplot",    # MCMC diagnostics
  # DAG + multiverse + bias analysis (Scripts 15, 24, 26)
  "dagitty",      # DAG specification
  "ggdag",        # DAG ggplot rendering
  "specr",        # 144-spec multiverse (Steegen 2016)
  "episensr",     # probabilistic bias analysis (Lash 2009)
  # Misc utilities
  "stringr",      # string operations
  "httr2",        # API calls (Crossref / NCHS)
  "jsonlite"      # JSON parsing
)

# 检查并安装缺失的包
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("正在安装: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(paste0("已安装: ", pkg, "\n"))
  }
}

cat("\n========================================\n")
cat("所有 R 包安装完成！\n")
cat("========================================\n")
