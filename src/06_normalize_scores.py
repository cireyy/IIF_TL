# min-max normalize scores to [-1, 1] so PRS and Transformer predictions are on the same scale

import pandas as pd
from sklearn.preprocessing import MinMaxScaler


SCORE_COLS = ["TL_PRS_Score", "TL_Transformer"]


def normalize_scores(
    input_path: str,
    output_path: str,
    score_cols: list[str] | None = None,
) -> pd.DataFrame:
    if score_cols is None:
        score_cols = SCORE_COLS

    df = pd.read_csv(input_path)

    # Only normalize columns that actually exist in the file
    cols_to_scale = [c for c in score_cols if c in df.columns]
    if not cols_to_scale:
        raise ValueError(f"None of {score_cols} found in {input_path}")

    scaler = MinMaxScaler(feature_range=(-1, 1))
    df[cols_to_scale] = scaler.fit_transform(df[cols_to_scale])

    df.to_csv(output_path, index=False)
    print(f"Normalized columns {cols_to_scale} -> {output_path}")
    return df


if __name__ == "__main__":
    normalize_scores(
        input_path="results/IIF_Endo_TL_Score.csv",
        output_path="results/IIF_Endo_TL_Score_Normalized.csv",
    )
