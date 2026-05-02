# baseline PRS on the UKBB cohort
# encodes genotypes as 0/1/2 effect allele counts and dots with GWAS weights

import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix,
)

N_CASES = 2379  # cases are at the top of the file


def encode_genotype(geno_str: str, effect_allele: str) -> int:
    """count effect alleles (0/1/2) from a genotype like 'A,G'"""
    try:
        if pd.isna(geno_str) or str(geno_str).strip() == "":
            return 0
        a1, a2 = str(geno_str).split(",")
        return int(a1 == effect_allele) + int(a2 == effect_allele)
    except Exception:
        return 0


def run_ukbb_prs(gwas_path: str, geno_path: str) -> None:
    # load
    gwas = pd.read_excel(gwas_path)
    gwas["chr_pos"] = "chr" + gwas["hm_chr"].astype(str) + "_" + gwas["hm_pos"].astype(str)

    geno = pd.read_csv(geno_path)

    # keep only SNPs present in both files
    common_snps = list(set(gwas["chr_pos"]).intersection(set(geno.columns)))
    gwas_filtered = gwas[gwas["chr_pos"].isin(common_snps)].set_index("chr_pos")
    geno_filtered = geno[common_snps].copy()

    effect_allele_map = gwas_filtered["effect_allele"].to_dict()

    # encode
    for col in geno_filtered.columns:
        ea = effect_allele_map[col]
        geno_filtered[col] = geno_filtered[col].astype(str).apply(
            lambda x, ea=ea: encode_genotype(x, ea)
        )

    # PRS
    weights = gwas_filtered.loc[geno_filtered.columns, "effect_weight"].astype(float).values
    prs_scores = geno_filtered.dot(weights)

    # labels
    labels = np.zeros(len(prs_scores))
    labels[:N_CASES] = 1

    # evaluate
    threshold = np.median(prs_scores)
    pred_binary = (prs_scores >= threshold).astype(int)

    auc       = roc_auc_score(labels, prs_scores)
    precision = precision_score(labels, pred_binary)
    recall    = recall_score(labels, pred_binary)
    f1        = f1_score(labels, pred_binary)
    tn, fp, fn, tp = confusion_matrix(labels, pred_binary).ravel()

    print("=== UKBB Endometriosis Baseline PRS ===")
    print(f"AUC:       {auc:.3f}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1:        {f1:.3f}")
    print(f"TP={tp}  TN={tn}  FP={fp}  FN={fn}")

    # enrichment: top vs bottom 10%
    top10_mask    = prs_scores >= np.percentile(prs_scores, 90)
    bottom10_mask = prs_scores <= np.percentile(prs_scores, 10)
    top10_rate    = labels[top10_mask].mean()
    bottom10_rate = labels[bottom10_mask].mean()
    fold_change   = top10_rate / bottom10_rate if bottom10_rate > 0 else float("nan")

    print(f"\nEnrichment — Top 10%: {top10_rate:.2%}  |  Bottom 10%: {bottom10_rate:.2%}")
    print(f"Fold-change (Top / Bottom): {fold_change:.2f}x")


if __name__ == "__main__":
    run_ukbb_prs(
        gwas_path="data/uf/gwas/gwas_uf_ukb.xlsx",
        geno_path="data/uf/sample/uf_ukb.csv",
    )
