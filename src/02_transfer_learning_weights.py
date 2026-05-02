# compute per-SNP transfer learning weights based on MAF similarity
# between UKBB (source) and IIF (target)
# w_j = exp(-lambda * |MAF_src - MAF_tgt|), then multiply onto the UKBB effect size

import pandas as pd
import numpy as np

LAMBDA = 10


def compute_tl_weights(
    gwas_ukbb_path: str,
    gwas_iif_path: str,
    maf_ukbb_path: str,
    maf_iif_path: str,
    output_path: str,
    lambda_val: float = LAMBDA,
) -> pd.DataFrame:

    # load everything
    gwas_ukbb = pd.read_excel(gwas_ukbb_path)
    gwas_iif  = pd.read_excel(gwas_iif_path)
    maf_ukbb  = pd.read_csv(maf_ukbb_path).rename(columns={"MAF": "maf_ukbb"})
    maf_iif   = pd.read_csv(maf_iif_path).rename(columns={"MAF": "maf_iif"})

    # build chr_pos key
    for df in [gwas_ukbb, gwas_iif]:
        df["chr_pos"] = (
            "chr" + df["hm_chr"].astype(str) + "_" + df["hm_pos"].astype(int).astype(str)
        )

    # merge
    df = (
        gwas_ukbb
        .merge(gwas_iif, on="chr_pos", suffixes=("_ukbb", "_iif"))
        .merge(maf_ukbb, on="chr_pos")
        .merge(maf_iif,  on="chr_pos")
    )

    # weights and adjusted betas
    df["w_i"] = np.exp(-lambda_val * np.abs(df["maf_ukbb"] - df["maf_iif"]))
    df["beta_final"] = df["effect_weight_ukbb"] * df["w_i"]

    # save
    cols_out = [
        "chr_pos", "rsID_ukbb",
        "effect_weight_ukbb", "effect_weight_iif",
        "maf_ukbb", "maf_iif", "w_i", "beta_final",
    ]
    df_out = df[cols_out]
    df_out.to_csv(output_path, index=False)
    print(f"Saved TL weights ({len(df_out)} SNPs) -> {output_path}")
    return df_out


if __name__ == "__main__":
    compute_tl_weights(
        gwas_ukbb_path="data/endo/gwas/gwas_endo_ukb.xlsx",
        gwas_iif_path="data/endo/gwas/gwas_endo_iif.xlsx",
        maf_ukbb_path="data/endo/maf_ukb.csv",
        maf_iif_path="data/endo/maf_iif.csv",
        output_path="results/gwas_transfer_learning_result.csv",
    )
