# TL-adjusted PRS for the IIF cohort using weights from step 02
# uses beta_final = w_j * beta_src instead of raw GWAS effect sizes

import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix,
)

CLINICAL_FILE = "data/merged_clinical.xlsx"


def encode_genotype(geno_str: str, effect_allele: str) -> int:
    try:
        if pd.isna(geno_str) or str(geno_str).strip() == "":
            return 0
        a1, a2 = str(geno_str).split(",")
        return int(a1 == effect_allele) + int(a2 == effect_allele)
    except Exception:
        return 0


def run_iif_tl_prs(
    tl_weights_path: str,
    gwas_iif_path: str,
    geno_path: str,
    clinical_path: str,
    output_path: str = "results/IIF_Endo_TL_Score.csv",
) -> None:
    # load TL weights and attach effect_allele from IIF GWAS
    gwas_tl = pd.read_csv(tl_weights_path, index_col=0)

    gwas_iif = pd.read_excel(gwas_iif_path)
    gwas_iif["chr_pos"] = (
        "chr" + gwas_iif["hm_chr"].astype(str) + "_" + gwas_iif["hm_pos"].astype(str)
    )
    gwas_iif = gwas_iif.set_index("chr_pos")
    gwas_tl["effect_allele"] = gwas_iif["effect_allele"]

    # load genotype and clinical data
    geno       = pd.read_csv(geno_path)
    merged_df  = pd.read_excel(clinical_path)

    # filter females
    merged_female = merged_df[merged_df["Couple ID"].astype(str).str.endswith("F")].copy()
    merged_female["VPS"]  = merged_female["VPS"].astype(str).str.strip().str.upper()
    merged_female["Name"] = merged_female["Name"].astype(str).str.strip().str.upper()

    geno["Sample ID"] = geno["Sample ID"].astype(str).str.strip().str.upper()
    sample_ids = geno["Sample ID"]
    mask = sample_ids.isin(merged_female["VPS"]) | sample_ids.isin(merged_female["Name"])

    geno_f = geno[mask].copy()
    geno_f = geno_f.set_index("Sample ID")

    # SNP intersection and encoding
    common_snps = list(set(gwas_tl.index).intersection(set(geno_f.columns)))
    geno_filtered = geno_f[common_snps].copy()
    gwas_filtered = gwas_tl.loc[common_snps]
    effect_allele_map = gwas_filtered["effect_allele"].to_dict()

    for col in geno_filtered.columns:
        ea = effect_allele_map[col]
        geno_filtered[col] = geno_filtered[col].astype(str).apply(
            lambda x, ea=ea: encode_genotype(x, ea)
        )

    # build labels
    label_map: dict = {}
    for _, row in merged_female.iterrows():
        for key in [str(row["VPS"]).strip().upper(), str(row["Name"]).strip().upper()]:
            val = str(row["Endo"]).strip().upper()
            if val in ("Y", "N"):
                label_map[key] = val

    valid_samples, sample_labels = [], []
    for sid in geno_filtered.index:
        if sid in label_map:
            sample_labels.append(1 if label_map[sid] == "Y" else 0)
            valid_samples.append(sid)

    geno_final = geno_filtered.loc[valid_samples]
    labels     = np.array(sample_labels)

    # compute TL-adjusted PRS
    tl_weights = gwas_filtered.loc[geno_final.columns, "beta_final"].astype(float).values
    prs_scores = geno_final.dot(tl_weights).values

    # evaluate
    threshold   = np.median(prs_scores)
    pred_binary = (prs_scores >= threshold).astype(int)

    auc       = roc_auc_score(labels, prs_scores)
    precision = precision_score(labels, pred_binary)
    recall    = recall_score(labels, pred_binary)
    f1        = f1_score(labels, pred_binary)
    tn, fp, fn, tp = confusion_matrix(labels, pred_binary).ravel()

    print("=== IIF SNP-level TL + Linear PRS ===")
    print(f"AUC:       {auc:.3f}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1:        {f1:.3f}")
    print(f"TP={tp}  TN={tn}  FP={fp}  FN={fn}")

    # save
    pd.DataFrame({
        "Sample_ID": valid_samples,
        "TL_PRS_Score": prs_scores,
        "Label": labels,
    }).to_csv(output_path, index=False)
    print(f"Saved -> {output_path}")


if __name__ == "__main__":
    run_iif_tl_prs(
        tl_weights_path="results/gwas_transfer_learning_result.csv",
        gwas_iif_path="data/endo/gwas/gwas_endo_iif.xlsx",
        geno_path="data/endo/sample/endo_iif.csv",
        clinical_path=CLINICAL_FILE,
    )
