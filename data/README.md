# Data Directory

**Raw data files are not included in this repository** due to participant privacy and data governance restrictions.

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

- **UK Biobank**: Apply via [https://www.ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access).
- **GWAS summary statistics**: Publicly available from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) (search accession numbers provided in the paper).
