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
  "tidyverse",    # 数据处理全家桶（dplyr, tidyr, ggplot2, purrr, stringr 等）
  "survey",       # 复杂抽样设计与加权分析
  "haven",        # 读取 .xpt (SAS transport) 文件
  "broom",        # 提取回归模型结果为 tidy 格式
  "tableone",     # 快速生成基线特征表 (Table 1)
  "openxlsx",     # 输出 Excel 文件
  "DiagrammeR",   # 绘制流程图 (Figure 1)
  "corrplot"      # 相关矩阵可视化（备用）
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
