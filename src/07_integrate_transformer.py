# merge MetaGeno predicted probabilities into the TL-PRS table
# join on the numeric portion of sample IDs since the two files use different ID formats

import re
import pandas as pd


def extract_id_num(s: str) -> str | None:
    """pull out the numeric part of a sample ID"""
    nums = re.findall(r"\d+", str(s))
    return nums[0] if nums else None


def integrate_transformer(
    tl_score_path: str,
    transformer_pred_path: str,
    output_path: str,
    pred_col: str = "predicted_probability_normalized_normalized",
) -> pd.DataFrame:
    df_iif  = pd.read_csv(tl_score_path)
    df_pred = pd.read_csv(transformer_pred_path)

    # join on numeric ID
    df_iif["_id_num"]  = df_iif["Sample_ID"].apply(extract_id_num)
    df_pred["_id_num"] = df_pred.iloc[:, 0].apply(extract_id_num)

    # merge and clean up
    df_pred_slim = df_pred[["_id_num", pred_col]].rename(
        columns={pred_col: "MetaGeno_Prob"}
    )

    merged = df_iif.merge(df_pred_slim, on="_id_num", how="left")
    merged.drop(columns=["_id_num"], inplace=True)

    merged.to_csv(output_path, index=False)
    print(f"Merged {len(merged)} samples -> {output_path}")
    print(f"  MetaGeno predictions matched: {merged['MetaGeno_Prob'].notna().sum()}")
    return merged


if __name__ == "__main__":
    integrate_transformer(
        tl_score_path="results/IIF_Endo_TL_Score_Normalized.csv",
        transformer_pred_path="data/1kg_23_predictions_ENDO.csv",
        output_path="results/IIF_Endo_TL_Score_Transformer.csv",
    )
