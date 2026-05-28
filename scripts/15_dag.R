# ============================================
# 009 / scripts/15_dag.R
# DAG (Directed Acyclic Graph) for Phthalate -> IR -> Mortality
#
# 用 ggdag + dagitty:
#   - Exposure: Phthalate
#   - Outcome (primary): IR (HOMA-IR)
#   - Outcome (secondary): All-cause mortality
#   - Covariates: age, sex, race, SES (edu+pir), smoking, BMI/waist, diet,
#                 hypertension
#   - Get minimal adjustment set via dagitty::adjustmentSets()
#
# 输出: output/figures/dag_phth_ir_mortality.png
#        output/figures/dag_adjustment_set.txt
# ============================================

suppressPackageStartupMessages({
  library(dagitty); library(ggdag); library(ggplot2)
})

cat("========================================\n")
cat("009 / 15_dag.R: DAG + minimal adjustment set\n")
cat("========================================\n\n")

# ------------------------------------------------------------------
# Build DAG (dagitty syntax)
# ------------------------------------------------------------------
# 节点角色:
#   E (exposure)    = Phthalate
#   O (outcome)     = IR (primary)
#   O2 (outcome)    = Mortality (secondary)
#   Mediators       = BMI/waist (adiposity), hypertension
#   Confounders     = age, sex, race, SES (edu+pir), smoking, diet
dag <- dagitty('
dag {
  Phthalate [exposure]
  IR [outcome]
  Mortality [outcome]
  age [pos="0,0"]
  sex [pos="0,1"]
  race [pos="0,2"]
  SES [pos="0,3"]
  smoking [pos="0,4"]
  diet [pos="0,5"]
  Phthalate [pos="1,2.5"]
  BMI_waist [pos="2,1.5"]
  hypertension [pos="2,3.5"]
  IR [pos="3,2.5"]
  Mortality [pos="4,2.5"]

  age -> Phthalate
  age -> IR
  age -> Mortality
  age -> BMI_waist
  age -> hypertension
  sex -> Phthalate
  sex -> IR
  sex -> Mortality
  sex -> BMI_waist
  race -> Phthalate
  race -> IR
  race -> Mortality
  race -> SES
  SES -> Phthalate
  SES -> IR
  SES -> Mortality
  SES -> diet
  SES -> smoking
  smoking -> Phthalate
  smoking -> IR
  smoking -> Mortality
  smoking -> hypertension
  diet -> Phthalate
  diet -> BMI_waist
  diet -> IR
  Phthalate -> BMI_waist
  Phthalate -> hypertension
  Phthalate -> IR
  Phthalate -> Mortality
  BMI_waist -> IR
  BMI_waist -> hypertension
  BMI_waist -> Mortality
  hypertension -> IR
  hypertension -> Mortality
  IR -> Mortality
}
')

# ------------------------------------------------------------------
# Minimal adjustment sets
# ------------------------------------------------------------------
cat("\n--- Minimal adjustment set: Phthalate -> IR ---\n")
adj_ir <- adjustmentSets(dag, exposure = "Phthalate", outcome = "IR",
                         type = "minimal", effect = "total")
print(adj_ir)

cat("\n--- Minimal adjustment set: Phthalate -> Mortality (total effect) ---\n")
adj_mort <- adjustmentSets(dag, exposure = "Phthalate", outcome = "Mortality",
                           type = "minimal", effect = "total")
print(adj_mort)

cat("\n--- Minimal adjustment set: Phthalate -> Mortality (direct effect) ---\n")
adj_mort_direct <- adjustmentSets(dag, exposure = "Phthalate", outcome = "Mortality",
                                  type = "minimal", effect = "direct")
print(adj_mort_direct)

# Implied conditional independencies (Q: do we satisfy DAG?)
cat("\n--- Implied conditional independencies (first 8) ---\n")
ici <- impliedConditionalIndependencies(dag)
print(head(ici, 8))

# ------------------------------------------------------------------
# Plot via ggdag
# ------------------------------------------------------------------
cat("\n--- Plotting DAG ---\n")
tdag <- tidy_dagitty(dag)

p <- ggdag(tdag, text = FALSE) +
  geom_dag_point(aes(color = name), size = 18, alpha = 0.7) +
  geom_dag_text(color = "black", size = 3) +
  scale_color_manual(values = c(
    Phthalate    = "#E74C3C",
    IR           = "#3498DB",
    Mortality    = "#9B59B6",
    BMI_waist    = "#F39C12",
    hypertension = "#F39C12",
    age          = "#7F8C8D",
    sex          = "#7F8C8D",
    race         = "#7F8C8D",
    SES          = "#7F8C8D",
    smoking      = "#7F8C8D",
    diet         = "#7F8C8D"
  ), guide = "none") +
  labs(title = "DAG: Phthalate -> Insulin Resistance -> Mortality",
       subtitle = "Red = exposure | Blue/Purple = outcomes | Orange = mediators | Grey = confounders") +
  theme_dag() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/dag_phth_ir_mortality.png", p,
       width = 11, height = 7, dpi = 300, bg = "white")
cat("Figure saved: output/figures/dag_phth_ir_mortality.png\n")

# ------------------------------------------------------------------
# Adjustment-set text file
# ------------------------------------------------------------------
adj_lines <- c(
  "=== DAG: Phthalate -> IR -> Mortality ===",
  "",
  "Minimal adjustment set (Phthalate -> IR, total effect):",
  sprintf("  { %s }", paste(unlist(adj_ir), collapse = ", ")),
  "",
  "Minimal adjustment set (Phthalate -> Mortality, total effect):",
  sprintf("  { %s }", paste(unlist(adj_mort), collapse = ", ")),
  "",
  "Minimal adjustment set (Phthalate -> Mortality, direct effect, blocking BMI/HTN/IR mediator paths):",
  sprintf("  { %s }", paste(unlist(adj_mort_direct), collapse = ", ")),
  "",
  "Mediators (NOT in minimal adjustment set for total effect):",
  "  BMI_waist, hypertension, IR (for mortality)",
  "",
  "Confounders (in minimal adjustment set):",
  "  age, sex, race, SES (education + PIR), smoking, diet",
  "",
  "Note: 009 manuscript reports total-effect models adjusting for the minimal",
  "      confounder set + adiposity-blocking direct-effect models for mortality.",
  "      Mediation decomposition via CMAverse 4-way is reported separately."
)
writeLines(adj_lines, "output/figures/dag_adjustment_set.txt")
cat("Adjustment-set text saved: output/figures/dag_adjustment_set.txt\n")

cat("\n========================================\n")
cat("DAG complete.\n")
cat("========================================\n")
