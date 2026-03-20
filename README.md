# Consciousness Dynamics: A Computable Framework Identifies α₂ as a Transdiagnostic Protective Target

Analysis code for the HCF-4α framework validation study using NHANES data (2017-2023).

## 📌 For Peer Review

This repository contains all analysis scripts for the manuscript **"Consciousness Dynamics: A Computable Framework Identifies α₂ as a Transdiagnostic Protective Target"**.

---

## 📁 Repository Structure

```
NHANES-HCF-4alpha/
├── 01_data_cleaning_L.R
├── 01_data_cleaning_P.R
├── 02_data_merge_L.R
├── 02_data_merge_P.R
├── 03_variable_basic_L.R
├── 03_variable_basic_P.R
├── 04_variable_biomarker_L.R
├── 04_variable_biomarker_P.R
├── 08_final_merge_L.R
├── 08_final_merge_P.R
├── 09_paper1_HCF_analysis_L.R
├── 09_paper1_HCF_analysis_P.R
├── 10_paper1_supplement_L.R
├── 10_paper1_supplement_P.R
├── 11_paper1_deep_analysis_L.R
├── 11_paper1_deep_analysis_P.R
├── 12_paper2_pathway_analysis_L.R
├── 12_paper2_pathway_analysis_P.R
├── 13_paper2_supplement_L.R
├── 13_paper2_supplement_P.R
├── 14_paper2_deep_analysis_L.R
├── 14_paper2_deep_analysis_P.R
├── 15_paper2_deep_analysis_continued_L.R
├── 15_paper2_deep_analysis_continued_P.R
├── 16_paper3_alpha_analysis_L.R
├── 16_paper3_alpha_analysis_P.R
├── 17_cross_cycle_validation_PL.R
├── 18_longitudinal_analysis_PL.R
├── 19_supplement_analysis_PL.R
├── 20_sensitivity_analysis_PL.R
├── 21_serious_research.R
├── 22_sensitivity_imputation.R
├── run_all.R
├── config.R
├── README.md
│
├── data/
│   ├── raw/                # ⚠️ User must place NHANES XPT files here
│   └── processed/           # ✅ Pre-computed RDS files (05-07 outputs)
│       ├── L_cycle/         # 2021-2023 cycle
│       │   ├── alpha_factors.rds
│       │   ├── HCF_typing.rds
│       │   ├── pathway_clustering_results.rds
│       │   ├── pathway_final_with_clusters.rds
│       │   ├── master.rds
│       │   ├──analysis_dataset_subset.rds	
│       │   └── final_analysis_dataset.rds
│       └── P_cycle/         # 2017-2020 cycle
│           ├── alpha_factors_P.rds
│           ├── HCF_typing_P.rds
│           ├── pathway_clustering_results_P.rds
│           ├── pathway_final_with_clusters_P.rds
│           ├── master_P.rds
│           ├── analysis_dataset_subset_P.rds
│           └── final_analysis_dataset_P.rds
│
└── outputs/                 # Generated after running (not in repo)
    ├── logs/
    ├── results/
    └── figures/
```

---

## ⚠️ Important Note About Scripts 05-07

**Scripts 05, 06, 07 are NOT included** in this repository due to commercial patent considerations.

However, **their outputs are provided** as RDS files in `data/processed/`:
- `alpha_factors_*.rds` (from script 05)
- `HCF_typing_*.rds` (from script 06)
- `pathway_*_*.rds` (from script 07)

This allows full reproduction of scripts 08-22 without needing the original code.

---

## 🚀 How to Reproduce All Results

### 1. Prerequisites

Install required R packages:

```r
packages <- c(
  "tidyverse", "survey", "lavaan", "semTools", "psych",
  "ggplot2", "pheatmap", "pROC", "gridExtra", "emmeans",
  "jtools", "interactions", "corrplot", "irr", "rpart",
  "MatchIt", "WeightIt", "cluster", "splines", "rms",
  "poLCA", "flexmix", "lmtest", "sandwich", "gtools",
  "reshape2", "cobalt", "cocor", "lme4", "tableone",
  "partykit", "forestplot", "kableExtra", "knitr"
)

install.packages(packages)
```

### 2. Prepare NHANES Raw Data

Download XPT files from:
- [NHANES 2021-2023 (L cycle)](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2021)
- [NHANES 2017-2020 (P cycle)](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017)

Place all 40 XPT files in `data/raw/`:
```
data/raw/
├── DEMO_L.xpt
├── DPQ_L.xpt
├── ... (20 L cycle files)
├── P_DEMO.xpt
├── P_DPQ.xpt
└── ... (20 P cycle files)
```

### 3. One-Click Run

```r
# In R, simply run:
source("run_all.R")
```

This will sequentially execute all scripts 01-22 and generate all results.

### 4. View Results

All outputs are saved in:
- `outputs/results/` - Tables and analysis results
- `outputs/figures/` - Generated figures (PDF/PNG)
- `outputs/logs/` - Detailed execution logs

---

## 📋 Script Overview

| Scripts | Description |
|---------|-------------|
| 01-04 | Data cleaning and preparation (requires XPT files) |
| 05-07 | ⚠️ **Not included** (commercial) - outputs provided as RDS |
| 08 | Final data merging |
| 09-11 | Paper 1: HCF typology analyses |
| 12-15 | Paper 2: Pathway clustering and alpha factors |
| 16 | Paper 3: Alpha factor validation |
| 17-22 | Cross-cycle validation, sensitivity analyses |

---

## 📊 Key Findings

- **α₂ identified as transdiagnostic protective factor**  
  - ORs: 0.39–0.82 for depression, CVD, and diabetes  
  - Protective threshold: α₂ > -0.5  
  - Number Needed to Treat (NNT): **21** in general population, **4** in high-risk subgroup  

- **Perseveration Paradox**  
  Highest depression (PHQ-9 = 11.3) with lowest inflammation (CRP = 3.3 mg/L)  

- **Compression Hypothesis**  
  Fatigue ↑ 1.10, Suicidal ideation ↓ 0.40 — population-level evidence for α₁ awakening  

- **Incremental Value**  
  Alpha factors explain **+9.4% variance**, **4× physiological measures** (BMI, CRP, heart rate)  

- **Cross-Cycle Replicability**  
  All core findings replicated across two independent NHANES cycles (2017–2020 and 2021–2023)

---

## 🔗 Quick Links

- 📦 [GitHub Repository](https://github.com/liuxiangyu17git/Consciousness-Dynamics-Validation)
- 📋 [OSF Preregistration](https://osf.io/ap58t/overview?view_only=1edfed21b8874a249f592cbf691bdfc3)

---

## 📄 License

MIT License — see [LICENSE](LICENSE) file for details.

---

## 📬 Contact

Xiangyu Liu  
Independent Scholar  
Email: liuxiangyu@liuxiangyu.com.cn  
ORCID: 0009-0004-5650-7780