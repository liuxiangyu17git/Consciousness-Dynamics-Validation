---
title: Consciousness Dynamics Validation
layout: default
nav_order: 1
---

# Consciousness Dynamics: A Computable Framework Identifies α₂ as a Transdiagnostic Protective Target

Analysis code for the HCF-4α framework validation study using NHANES data (2017-2023).

## 📌 For Peer Review

This repository contains all analysis scripts for the manuscript **"Consciousness Dynamics: A Computable Framework Identifies α₂ as a Transdiagnostic Protective Target"**.

---

## 📁 Repository Structure

## 📁 Repository Structure

### 📊 Data Preparation (Scripts 01-08)
Scripts for cleaning and merging raw NHANES data (each has `_L` and `_P` versions for two cycles).

### 📄 Paper 1: HCF Typology (Scripts 09-11)
- `09_paper1_HCF_analysis_L.R` / `_P.R`: Main HCF typology analyses
- `10_paper1_supplement_L.R` / `_P.R`: Supplementary analyses for Paper 1
- `11_paper1_deep_analysis_L.R` / `_P.R`: Deep mediation and D-layer analyses

### 📄 Paper 2: Pathway Clustering (Scripts 12-15)
- `12_paper2_pathway_analysis_L.R` / `_P.R`: K-means clustering and pathway identification
- `13_paper2_supplement_L.R` / `_P.R`: Supplementary materials for Paper 2
- `14_paper2_deep_analysis_L.R` / `_P.R`: Pathway × alpha interactions
- `15_paper2_deep_analysis_continued_L.R` / `_P.R`: Triple interactions, inflammation trajectories

### 📄 Paper 3: Alpha Factor Validation (Scripts 16)
- `16_paper3_alpha_analysis_L.R` / `_P.R`: Main alpha factor validation analyses

### 🔄 Cross-Cycle & Sensitivity (Scripts 17-21)
- `17_cross_cycle_validation_PL.R`: ICC, weighted kappa, cross-cycle comparisons
- `18_longitudinal_analysis_PL.R`: Age-period-cohort, COVID natural experiment
- `19_supplement_analysis_PL.R`: Supplemental analyses (FDR, NNT, AUC, etc.)
- `20_sensitivity_analysis_PL.R`: Sensitivity analyses (medication, thresholds, youth subgroups)
- `21_serious research.R`: EFA, correlation matrices, CART validation

### 📁 Data Files
- `analysis_dataset_subset_L.rds`: Cleaned dataset for 2021-2023 cycle (N = 7,015)
- `analysis_dataset_subset_P.rds`: Cleaned dataset for 2017-2020 cycle (N = 5,843)

## 📊 Key Findings

- **α₂** identified as transdiagnostic protective factor  
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

## 🧪 How to Reproduce the Results

### 1. Prerequisites

Install required R packages:

```r
packages <- c("survey", "lavaan", "lavaan.survey", "mice", 
              "emmeans", "ggplot2", "gtsummary", "dplyr")
install.packages(packages)
```

### 2. Run the analysis
Place the `.rds` data files in your working directory, then run scripts sequentially:

```r
# Start with scripts 09–11 for HCF analyses
source("09_paper1_HCF_analysis_L.R")
source("10_paper1_supplement_L.R")
source("11_paper1_deep_analysis_L.R")

# Then scripts 12–15 for pathway analyses
source("12_paper2_pathway_analysis_L.R")
# ... continue with subsequent scripts
```
All outputs (tables and figures) will be generated in the working directory.

🔗 Quick Links
📦 GitHub Repository

📋 OSF Preregistration

📄 View README (full documentation)

📄 License
This project is licensed under the MIT License — see the LICENSE file for details.

📬 Contact
Xiangyu Liu
Independent Scholar
Email: liuxiangyu@liuxiangyu.com.cn
ORCID: 0009-0004-5650-7780
