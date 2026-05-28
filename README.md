# Urinary phthalate metabolites and adiposity-mediated insulin resistance (NHANES 2005–2018)

Reproducible analytic pipeline for the manuscript:

> **Urinary phthalate metabolites and adiposity-mediated insulin resistance: a multi-method triangulation analysis in NHANES 2005–2018**
> Submitted to *Diabetes Care*.

## What this repository contains

The full reproducible R pipeline for our analysis of ten urinary phthalate metabolites, adiposity-mediated insulin resistance (HOMA-IR), and all-cause + cardiometabolic mortality in **2,239 fasting US adults** from seven NHANES cycles (2005-2018, D-J), linked to the 2019 NCHS Public-Use Linked Mortality File.

- **Primary exposures**: ten urinary phthalate monoester metabolites (MEP, MnBP, MiBP, MBzP, MEHP, MEHHP, MEOHP, MECPP, MCNP, MCOP); primary mixture = molar-sum Σ-DEHP (4 DEHP metabolites)
- **Primary outcome**: insulin resistance (HOMA-IR ≥ 2.5, McAuley threshold)
- **Secondary outcome**: all-cause mortality (Linked Mortality File 2019, 141 events post 200-month PERMTH cap) + cardiometabolic mortality (52 events, treated as underpowered)
- **Mediator**: adiposity composite (BMI-*z* + waist-*z*)
- **Methods**: seven-framework triangulation — CMAverse 4-way (VanderWeele 2014), qgcomp (Keil 2020), WQS (Carrico 2015), BKMR (Bobb 2015/2018), IPTW (GBM + stabilised), Bayesian g-computation (rstanarm HMC), 144-spec multiverse (Steegen 2016)

## Key findings (pre-registered primary + secondary)

**Pre-registered primary**:
- CMAverse 4-way Q4-vs-Q1 total-effect OR for IR: **1.453 (95% CI 1.224, 1.777, *P* < 0.001)**, with **52.5% mediated by adiposity** (Rpnie OR 1.170, *P* < 0.001)
- Six of seven analytical frameworks pointed in the positive direction (qgcomp null due to sub-weight cancellation)
- Single-metabolite drivers: MEHP RCS *P*-nonlinear = 0.0005; MiBP per-SD OR 1.22 (*P*-BH = 0.011)

**Pre-registered secondary**:
- All-cause mortality Σ-DEHP-*z* HR **1.69 (1.33, 2.15, *P* < 0.001)**; the mortality pathway is direct-effect-dominated (Rcde = 1.823), not adiposity-mediated (mortality-pm = -0.003 *P* = 0.82)
- Cardiometabolic mortality (52 events) underpowered (Type-M 3-9×), treated as inconclusive

## Repository layout

```
scripts/                R analytic pipeline (run 00 → 33 + run_all + helpers)
  00_install_packages.R    one-time install
  01_download_data.R       download NHANES + mortality
  02_merge_data.R          merge across 7 cycles
  03_clean_data.R          phthalate + HOMA-IR + adiposity coding
  04_survey_design.R       svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_pooled)
  05_table1.R              Table 1 baseline
  06_bkmr_phth_homa_ir_checkpoint.R  BKMR 8-metabolite mixture
  07_cox_mortality.R       Cox + Schoenfeld PH check
  08_logistic_ir.R         logistic single-metabolite + Σ-DEHP
  09_rcs.R                 restricted cubic splines
  10_qgcomp.R              qgcomp mixture
  11_wqs.R                 WQS positive-direction
  12_cmaverse_4way.R       CMAverse 4-way decomposition
  13_iptw.R                IPTW (GBM + stabilised + 99% trim)
  14_evalue.R              VanderWeele-Ding E-value
  15_dag.R                 ggdag DAG
  16_subgroup.R            9-subgroup × interaction forest
  17_sensitivity.R         S1-S8 sensitivity analyses
  18_evalue_ci_bound.R     E-value CI bound
  19_tripod_ai.R           TRIPOD-AI 27-item compliance
  20_nri_idi.R             NRI + IDI
  21_dca.R                 decision curve analysis
  22_calibration.R         Hosmer-Lemeshow + Brier + calibration belt
  23_nomogram.R            nomogram + 5-tier risk stratification
  24_multiverse.R          144-spec multiverse
  25_bayesian_gcomp.R      Bayesian g-computation (rstanarm HMC)
  26_episensr.R            probabilistic bias analysis (Lash 2009)
  27_lag_analysis.R        negative-control exposure (Lipsitch 2010)
  28_ctat_methods_table.R  CTAT 4-quadrant methods summary
  29_fine_gray.R           Fine-Gray competing-risk for CM mortality
  30_reri_ap_si.R          RERI + AP + synergy index
  31_multiplicity_padjust.R  BH-FDR within 3 test families
  32_retrodesign.R         Type-S/M (Gelman-Carlin 2014)
  33_figures_rework.R      figure rendering
  run_all.R                end-to-end orchestrator
data/processed/         intermediate .RData files (start here to skip 01-03)
  nhanes_raw_merged.RData    raw merged across 7 cycles
  nhanes_final.RData         analytic sample N=2,239
  nhanes_design.RData        svydesign object
  bayesian_gcomp_009.RData   rstanarm fitted posterior
  episensr_009.RData         probabilistic bias-analysis MC draws
  tripod_ir_models.RData     4-model TRIPOD AUROC artefacts
output/
  tables/               primary + supplementary tables (CSV)
  figures/              CONSORT, CMAverse 4-way, multiverse, nomogram (PNG)
LICENSE                 MIT
.gitignore              excludes data/raw/ + .DS_Store + R IDE
```

## How to reproduce

Software: **R 4.6** + **pandoc** (for rendering, optional). All package install is in `scripts/00_install_packages.R`.

```r
# 1. Install packages (one-time, ~30-60 min including BKMR + CMAverse from GitHub)
Rscript scripts/00_install_packages.R

# 2. Download raw NHANES + mortality data (~15 min, ~700 MB to data/raw/)
Rscript scripts/01_download_data.R

# 3. Run the full pipeline (~6-10 h with BKMR; ~30 min without)
Rscript scripts/run_all.R
```

Or, to **skip the data download / merge / BKMR steps**, load the saved analytic sample directly:

```r
load("data/processed/nhanes_final.RData")   # → nhanes_final, N=2,239
load("data/processed/nhanes_design.RData")  # → nhanes_design, svydesign object
# then run scripts/05_table1.R onward
```

**Note on BKMR**: the BKMR posterior `output/tables/bkmr_phth_results.RData` (~311 MB) is **not** included in this repository (GitHub 100 MB per-file hard limit). To regenerate, run `scripts/06_bkmr_phth_homa_ir_checkpoint.R` (~6-10 h on 8-core machine; chain 1 + chain 2 × 10,000 iterations).

## Reproducibility notes

- The pipeline is deterministic. Section-stratified seeds (full table in Supplementary §S10):
  - `set.seed(20260523)` — Bayesian g-computation / cleaning / MICE / retrodesign
  - `set.seed(20260524)` — mixture and mediation bootstraps (qgcomp / WQS / CMAverse / IPTW)
  - `SEED_BASE + iter` — BKMR chain initialisation
- All survey-weighted analyses use the pooled fasting weight (`WTSAF2YR` scaled to 14 person-years across seven biennial cycles per NCHS Series 2 No. 190) with `survey::svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled, nest = TRUE)`.
- The CMAverse 4-way decomposition is unweighted by design (no native svydesign support); the weighted-vs-unweighted sensitivity (Supplementary §S25) shows ≤ 3.4% drift on the Q4-vs-Q1 contrast.
- BKMR partial-convergence (kernel-bandwidth *r*-parameter ESS 179-697 < 1,000 target with rhat 1.002-1.015 < 1.05) is recognised in the BKMR diagnostics literature [Bobb 2018 *Environ Health*]; posterior inclusion probabilities are stable under this partial state because PIP integrates over both posterior modes.

## Data

- **NHANES 2005–2018**: public domain. CDC/NCHS. https://wwwn.cdc.gov/nchs/nhanes/
- **NCHS Linked Mortality File 2019**: public-use release. https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

This repository does **not** include the raw `.XPT` NHANES files. `scripts/01_download_data.R` will download them into `data/raw/` (which is `.gitignore`-d).

## License

Code: MIT License (see `LICENSE`).
Data: NHANES is in the public domain (US government data).

## Citation

If you use this code, please cite:

> Li J, Sun X, Zhang J, Zhai L, Yu L. Urinary phthalate metabolites and adiposity-mediated insulin resistance: a multi-method triangulation analysis in NHANES 2005–2018. *Diabetes Care*. (Submitted 2026; volume / pages / DOI to be assigned upon acceptance.)

## Contact

**Corresponding author**: Ling Yu (yulingyxb@jlu.edu.cn, ORCID 0000-0001-7362-3581), Department of Pharmacy, The Second Hospital of Jilin University, Changchun, Jilin Province, China.

**Co-authors**:
- Jie Li (first author) — Department of Obstetrics and Gynecology, The Second Hospital of Jilin University
- Xiubo Sun, Jing Zhang, Lijie Zhai — Department of Pharmacy, The Second Hospital of Jilin University

## Funding

Jilin Provincial Higher Education Research Project (Grant No. JGJX2021D37), awarded to Ling Yu.

## Declaration of Generative AI and AI-assisted technologies in the writing process

AI-assisted writing tools were used in two scope-limited ways:

(i) code implementation for the post-protocol sensitivity layers (qgcomp, WQS, CMAverse, IPTW, Bayesian g-computation, episensr, multiverse) after the OSF v1.0 analytic plan was finalised by the authors;

(ii) sentence-level language polishing of the Methods and Results sections only.

The Background, Discussion, and Conclusions sections were authored without AI assistance. All study design choices, statistical model selection, scientific decisions, citation choices, numerical claims, and interpretations were independently made by the authors and verified against the underlying statistical outputs. Disclosure conforms to COPE 2023 guidance on authorship and AI tools.
