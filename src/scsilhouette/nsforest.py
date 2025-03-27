from pathlib import Path
import pandas as pd

def load_nsforest_fscore(nsforest_file: Path) -> pd.DataFrame:
    """
    Load NS-Forest F-scores from a marker gene CSV.

    Expects at least two columns:
    - cluster / label (used to match)
    - F-score or f_score
    """
    df = pd.read_csv(nsforest_file)
    df.columns = df.columns.str.lower()

    # Standardize column names
    if "f-score" in df.columns:
        df.rename(columns={"f-score": "f_score"}, inplace=True)

    return df[["cluster", "f_score"]].dropna()

