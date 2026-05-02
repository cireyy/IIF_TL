# cross-population baseline: apply GWAS weights directly to IIF with no adjustment

import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix,
)

CLINICAL_FILE = "data/merged_clinical.xlsx"  # set to your clinical metadata file


def encode_genotype(geno_str: str, effect_allele: str) -> int:
    try:
        if pd.isna(geno_str) or str(geno_str).strip() == "":
            return 0
        a1, a2 = str(geno_str).split(",")
        return int(a1 == effect_allele) + int(a2 == effect_allele)
    except Exception:
        return 0


def run_iif_baseline_prs(
    gwas_path: str,
    geno_path: str,
    clinical_path: str,
    output_path: str = "results/IIF_Endometriosis_PRS_scores_baseline.csv",
) -> None:
    # load
    gwas = pd.read_excel(gwas_path)
    gwas["chr_pos"] = "chr" + gwas["hm_chr"].astype(str) + "_" + gwas["hm_pos"].astype(str)

    geno = pd.read_csv(geno_path)
    merged_df = pd.read_excel(clinical_path)

    # keep only female samples (Couple ID ends with 'F')
    merged_female = merged_df[merged_df["Couple ID"].astype(str).str.endswith("F")].copy()
    merged_female["VPS"]  = merged_female["VPS"].astype(str).str.strip().str.upper()
    merged_female["Name"] = merged_female["Name"].astype(str).str.strip().str.upper()

    sample_ids = geno["Sample ID"].astype(str).str.strip().str.upper()
    mask = sample_ids.isin(merged_female["VPS"]) | sample_ids.isin(merged_female["Name"])
    geno_f = geno[mask].copy()

    # build labels (Y=1, N=0, unknown=nan)
    label_map_vps  = merged_female.set_index("VPS")["Endo"].to_dict()
    label_map_name = merged_female.set_index("Name")["Endo"].to_dict()

    labels = []
    for sid in sample_ids[mask]:
        val = label_map_vps.get(sid, label_map_name.get(sid, np.nan))
        if str(val).strip().upper() == "Y":
            labels.append(1)
        elif str(val).strip().upper() == "N":
            labels.append(0)
        else:
            labels.append(np.nan)
    labels = np.array(labels)

    # encode genotypes
    common_snps = list(set(gwas["chr_pos"]).intersection(set(geno.columns)))
    gwas_filtered = gwas[gwas["chr_pos"].isin(common_snps)].set_index("chr_pos")
    geno_filtered = geno_f[common_snps].copy()

    effect_allele_map = gwas_filtered["effect_allele"].to_dict()
    for col in geno_filtered.columns:
        ea = effect_allele_map[col]
        geno_filtered[col] = geno_filtered[col].astype(str).apply(
            lambda x, ea=ea: encode_genotype(x, ea)
        )

    # PRS
    weights    = gwas_filtered.loc[geno_filtered.columns, "effect_weight"].astype(float).values
    prs_scores = geno_filtered.dot(weights)

    # drop samples with unknown label
    valid      = ~np.isnan(labels)
    labels_v   = labels[valid]
    prs_v      = prs_scores[valid]
    sample_ids_v = sample_ids[mask].values[valid]

    # evaluate
    threshold    = np.median(prs_v)
    pred_binary  = (prs_v >= threshold).astype(int)

    auc       = roc_auc_score(labels_v, prs_v)
    precision = precision_score(labels_v, pred_binary)
    recall    = recall_score(labels_v, pred_binary)
    f1        = f1_score(labels_v, pred_binary)
    tn, fp, fn, tp = confusion_matrix(labels_v, pred_binary).ravel()

    print("=== IIF Baseline PRS (no TL) ===")
    print(f"AUC:       {auc:.3f}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1:        {f1:.3f}")
    print(f"TP={tp}  TN={tn}  FP={fp}  FN={fn}")

    # save per-sample scores
    pd.DataFrame({
        "Sample_ID": sample_ids_v,
        "PRS_Score": prs_v,
        "Label": labels_v,
    }).to_csv(output_path, index=False)
    print(f"Saved -> {output_path}")


if __name__ == "__main__":
    run_iif_baseline_prs(
        gwas_path="data/endo/gwas/gwas_endo_iif.xlsx",
        geno_path="data/endo/sample/endo_iif.csv",
        clinical_path=CLINICAL_FILE,
    )
