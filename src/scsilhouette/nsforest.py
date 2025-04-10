from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import json


def load_nsforest_binary_genes(
    nsforest_file: Path, dump_json: bool = False, json_path: Path = None
) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    """
    Load NS-Forest binary gene definitions and F-score data.

    Args:
        nsforest_file: Path to the NS-Forest CSV file.
        dump_json: Whether to save the binary gene dictionary to JSON.
        json_path: Optional path to save JSON if dump_json is True.

    Returns:
        A tuple:
        - Dictionary mapping cluster to list of binary genes
        - DataFrame with columns: cluster, F-score, binary_genes
    """
    df = pd.read_csv(nsforest_file)

    if 'clusterName' not in df.columns or 'binary_genes' not in df.columns:
        raise ValueError("NS-Forest CSV must include 'clusterName' and 'binary_genes' columns.")

    df = df.rename(columns={
        'clusterName': 'cluster',
        'f_score': 'F-score'
    })

    gene_dict = {}
    for _, row in df.iterrows():
        cluster = str(row['cluster'])
        genes = str(row['binary_genes']).split(';') if pd.notna(row['binary_genes']) else []
        gene_dict[cluster] = [g.strip() for g in genes if g.strip()]

    if dump_json:
        if json_path is None:
            json_path = Path(nsforest_file).with_suffix(".binary_genes.json")
        with open(json_path, "w") as f:
            json.dump(gene_dict, f, indent=2)

    return gene_dict, df[['cluster', 'F-score', 'binary_genes']]


def validate_cluster_alignment(
    silhouette_df: pd.DataFrame, fscore_df: pd.DataFrame, label_key: str
) -> List[str]:
    """
    Identify unmatched cluster labels between silhouette and NS-Forest data.

    Args:
        silhouette_df: DataFrame with silhouette scores.
        fscore_df: DataFrame with NS-Forest F-scores.
        label_key: Column name for cluster labels in silhouette_df.

    Returns:
        List of unmatched labels.
    """
    silhouette_labels = set(silhouette_df[label_key].astype(str).unique())
    fscore_labels = set(fscore_df['cluster'].astype(str).unique())
    return sorted(list(silhouette_labels.symmetric_difference(fscore_labels)))

