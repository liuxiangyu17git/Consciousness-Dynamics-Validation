# HCF-4α Framework Analysis Code

## For Peer Review Only

This repository contains analysis code for the manuscript 
*"H = R_carrier × R_consciousness: A Computable Theory of Psychophysiological Health"*.

---

## 📁 Repository Structure

| Folder | Contents |
|--------|----------|
| `/` | Analysis scripts 01-20 (excluding 05-06) |
| `/data` | Final analysis datasets (`analysis_dataset_subset_L.rds`, `analysis_dataset_subset_P.rds`) |

---

## 📊 Data Files

The `/data` folder contains:
- `analysis_dataset_subset_L.rds`: 2021-2023 cycle (N = 7,015)
- `analysis_dataset_subset_P.rds`: 2017-2020 cycle (N = 5,843)

These datasets include all variables needed to replicate scripts 09-20:
- Alpha factors (alpha1-alpha4, standardized)
- HCF types and pathway clusters
- Health outcomes and physiological measures
- Demographics and survey design variables

---

## 🧪 Running the Code

1. **Copy data files to working directory**
   ```r
   file.copy("data/analysis_dataset_subset_L.rds", "./")
