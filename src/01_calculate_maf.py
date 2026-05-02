# calculate MAF for each SNP from the IIF genotype file
# genotypes come in as "A,G" format so we split and count alleles

import pandas as pd
from collections import Counter


def calculate_maf(geno_path: str, output_path: str) -> pd.DataFrame:
    df = pd.read_csv(geno_path)
    df_snps = df.drop(columns=["Sample ID"])

    maf_results = []
    for snp in df_snps.columns:
        genotypes = df_snps[snp].dropna()

        alleles = []
        for genotype in genotypes:
            alleles.extend(str(genotype).split(","))  # "A,G" -> ["A","G"]

        counts = Counter(alleles)
        total = sum(counts.values())

        if total == 0:
            maf = None
        else:
            freqs = [v / total for v in counts.values()]
            maf = min(freqs)

        maf_results.append({"chr_pos": snp, "MAF": maf})

    maf_df = pd.DataFrame(maf_results)
    maf_df.to_csv(output_path, index=False)
    print(f"Saved MAF table ({len(maf_df)} SNPs) -> {output_path}")
    return maf_df


if __name__ == "__main__":
    calculate_maf(
        geno_path="data/endo/sample/endo_iif.csv",
        output_path="results/calculated_maf.csv",
    )
