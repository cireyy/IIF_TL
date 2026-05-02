# Data Directory

**Raw data files are not included in this repository** due to participant privacy and data governance restrictions.

This directory shows the expected folder structure and file formats needed to run the analysis.

---

## Expected Structure

```
data/
├── endo/
│   ├── gwas/
│   │   ├── gwas_endo_ukb.xlsx      # UKBB endometriosis GWAS summary statistics
│   │   └── gwas_endo_iif.xlsx      # IIF endometriosis GWAS summary statistics
│   ├── sample/
│   │   └── endo_iif.csv            # IIF genotype matrix
│   ├── maf_ukb.csv                 # UKBB SNP minor allele frequencies
│   └── maf_iif.csv                 # IIF SNP minor allele frequencies
├── uf/
│   ├── gwas/
│   │   ├── gwas_uf_ukb.xlsx        # UKBB uterine fibroids GWAS (optional)
│   │   └── gwas_uf_iif.xlsx        # IIF uterine fibroids GWAS (optional)
│   └── sample/
│       └── uf_ukb.csv              # UKBB genotype matrix (uterine fibroids)
├── merged_clinical.xlsx            # IIF clinical metadata (Endo labels, Couple IDs)
└── 1kg_23_predictions_ENDO.csv     # MetaGeno predicted probabilities (1000 Genomes)
```

---

## File Format Specifications

### GWAS summary statistics (`gwas_endo_*.xlsx`)

| Column | Description |
|---|---|
| `rsID` | SNP identifier |
| `chr_name` | Chromosome |
| `chr_position` | Base-pair position (GRCh37) |
| `effect_allele` | Effect allele |
| `other_allele` | Non-effect allele |
| `effect_weight` | Log-odds / beta coefficient |
| `hm_chr` | Harmonised chromosome |
| `hm_pos` | Harmonised position |
| `hm_rsID` | Harmonised rsID |

### Genotype matrix (`endo_iif.csv`, `uf_ukb.csv`)

- Row per sample, column per SNP
- First column: `Sample ID`
- SNP columns named as `chr<CHR>_<POS>` (e.g. `chr1_5040823`)
- Genotype format: `"A,G"` (two alleles, comma-separated)

### MAF tables (`maf_iif.csv`, `maf_ukb.csv`)

| Column | Description |
|---|---|
| `chr_pos` | SNP identifier (`chr<CHR>_<POS>`) |
| `MAF` | Minor allele frequency (0–0.5) |

### Clinical metadata (`merged_clinical.xlsx`)

| Column | Description |
|---|---|
| `VPS` | Sample identifier (VPS number) |
| `Name` | Alternative sample identifier |
| `Couple ID` | Couple identifier; female samples end with `'F'` |
| `Endo` | Endometriosis status: `Y` = case, `N` = control, `U` = unknown |

---

## Data Access

- **IIF cohort**: Contact the corresponding author. Data sharing is subject to ethics approval and a Data Access Agreement.
- **UK Biobank**: Apply via [https://www.ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access).
- **GWAS summary statistics**: Publicly available from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) (search accession numbers provided in the paper).
