# Cross-Ancestry Transfer Learning for Endometriosis Polygenic Risk Prediction

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Code repository for the paper:

> **Cross-Ancestry Transfer Learning for Polygenic Risk Prediction of Endometriosis**  
> *[Author list — to be updated before submission]*

---

## Overview

Standard polygenic risk scores (PRS) are trained on large European-ancestry GWAS cohorts (e.g., UK Biobank) and lose predictive accuracy when applied to genetically diverse populations. This repository implements two transfer-learning strategies that improve cross-population generalisability:

| Method | Description |
|---|---|
| **SNP-level TL + Linear PRS** | Per-SNP exponential weight based on MAF similarity between source and target populations |
| **MetaGeno + TL-Loss** | Transformer-based deep genomic model with a cross-ancestry embedding alignment loss |

### Performance (Endometriosis)

| Dataset | Method | AUC | Precision | Recall | F1 |
|---|---|---|---|---|---|
| UKBB | Original GWAS + Linear PRS | 0.632 | 0.528 | 0.602 | 0.562 |
| UKBB | MetaGeno | 0.655 | 0.545 | 0.625 | 0.582 |
| IVF (IIF) | Original GWAS + Linear PRS | 0.625 | 0.714 | 0.455 | 0.556 |
| IVF (IIF) | **SNP-level TL + Linear PRS** | **0.659** | **0.875** | **0.636** | **0.737** |
| IVF (IIF) | **MetaGeno + TL-Loss** | **0.670** | **0.857** | **0.727** | **0.786** |

---

## Repository Structure

```
IIF_TL_GitHub/
├── src/
│   ├── 01_calculate_maf.py            # Compute SNP MAF from IIF genotypes
│   ├── 02_transfer_learning_weights.py # SNP-level TL weight computation
│   ├── 03_prs_ukbb.py                 # Baseline PRS evaluation on UKBB
│   ├── 04_prs_iif_baseline.py         # Cross-population baseline PRS on IIF
│   ├── 05_prs_iif_tl.py               # TL-adjusted PRS on IIF
│   ├── 06_normalize_scores.py         # Min-max normalisation of scores
│   └── 07_integrate_transformer.py    # Merge MetaGeno predictions
├── data/                              # ⚠ NOT included (see Data Access below)
│   ├── endo/
│   │   ├── gwas/                      # GWAS summary statistics (.xlsx)
│   │   ├── sample/                    # Genotype matrix (.csv)
│   │   ├── maf_iif.csv
│   │   └── maf_ukb.csv
│   └── uf/
│       ├── gwas/
│       └── sample/
├── results/                           # Output files written here
├── TL.ipynb                           # Full analysis notebook
├── requirements.txt
└── .gitignore
```

---


## Installation

```bash
git clone https://github.com/<your-username>/IIF_TL.git
cd IIF_TL
pip install -r requirements.txt
```

Python 3.9 or higher is required.

---

## Usage

Run the scripts in order from the repository root:

```bash
# Step 1 — Compute IIF cohort MAF
python src/01_calculate_maf.py

# Step 2 — Compute SNP-level TL weights
python src/02_transfer_learning_weights.py

# Step 3 — Evaluate baseline PRS on UKBB
python src/03_prs_ukbb.py

# Step 4 — Evaluate cross-population baseline PRS on IIF
python src/04_prs_iif_baseline.py

# Step 5 — Evaluate TL-adjusted PRS on IIF
python src/05_prs_iif_tl.py

# Step 6 — Normalise scores
python src/06_normalize_scores.py

# Step 7 — Integrate MetaGeno predictions
python src/07_integrate_transformer.py
```

Alternatively, run the full pipeline interactively in `TL.ipynb`.

---

## Data Access

Raw genotype data and clinical records cannot be shared publicly due to participant privacy. Researchers wishing to access the data should contact:

- **UK Biobank**: Apply via [https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)

GWAS summary statistics used in this study are publicly available from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

See `data/README.md` for the expected file formats.

---
