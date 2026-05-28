# 工具脚本：从 output/tables 抠出关键数字，输出到控制台
# 用于撰写"结果段"
suppressPackageStartupMessages({library(dplyr)})

t3 <- read.csv("output/tables/table3_cox.csv")
t4 <- read.csv("output/tables/table4_logistic.csv")
pt3 <- read.csv("output/tables/table3_cox_ptrend.csv")
pt4 <- read.csv("output/tables/table4_logistic_ptrend.csv")
pv  <- read.csv("output/tables/predictive_value.csv")
md  <- read.csv("output/tables/mediation.csv")

cat("=== T3 (Cox) Q4 ALL ===\n")
print(t3[grepl("Q4", t3$term), c("exposure","status","model","HR_CI","p")], row.names = FALSE)

cat("\n=== T4 (Logit) Q4 Model 2 ===\n")
print(t4[grepl("Q4", t4$term) & t4$model == "Model2",
         c("exposure","outcome","OR_CI","p")], row.names = FALSE)

cat("\n=== T4 (Logit) Q4 Model 1 ===\n")
print(t4[grepl("Q4", t4$term) & t4$model == "Model1",
         c("exposure","outcome","OR_CI","p")], row.names = FALSE)

cat("\n=== P-trend (Cox) ===\n")
print(pt3, row.names = FALSE)
cat("\n=== P-trend (Logit) ===\n")
print(pt4, row.names = FALSE)

cat("\n=== Predictive value (C-stat + p) ===\n")
print(pv[, c("outcome","model","Cstat","p_C")], row.names = FALSE)

cat("\n=== 中介（Prop p<0.05） ===\n")
m2 <- md[!is.na(md$Prop_p) & md$Prop_p < 0.05,
         c("exposure","mediator","outcome","Prop_med","Prop_p")]
print(m2, row.names = FALSE)

cat("\n=== 事件数 ===\n")
load("data/processed/nhanes_design.RData")
cat(sprintf("N = %d\n", nrow(nhanes_final)))
cat(sprintf("Death all = %d, CVD death = %d\n",
            sum(nhanes_final$death_all == 1),
            sum(nhanes_final$death_cvd == 1)))
cat(sprintf("Total CVD=%d, CHF=%d, CHD=%d, Angina=%d, MI=%d, Stroke=%d\n",
            sum(nhanes_final$total_cvd_y == 1),
            sum(nhanes_final$chf_y == 1),
            sum(nhanes_final$chd_y == 1),
            sum(nhanes_final$angina_y == 1),
            sum(nhanes_final$mi_y == 1),
            sum(nhanes_final$stroke_y == 1)))
