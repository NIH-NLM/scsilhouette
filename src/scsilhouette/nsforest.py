from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import json
import ast


def load_nsforest_binary_genes(
    nsforest_file: Path, dump_json: bool = False, json_path: Path = None
) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
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


def extract_binary_genes(nsforest_path: str, output_path: str):
    df = pd.read_csv(nsforest_path)

    if "binary_genes" not in df.columns:
        raise ValueError("Missing 'binary_genes' column in NS-Forest file.")

    all_genes = set()
    for item in df["binary_genes"]:
        if pd.isna(item):
            continue
        try:
            genes = ast.literal_eval(item)
            if isinstance(genes, list):
                all_genes.update(genes)
        except Exception as e:
            raise ValueError(f"Failed to parse binary_genes entry: {item}") from e

    with open(output_path, "w") as f:
        for gene in sorted(all_genes):
            f.write(f"{gene}\n")


def load_binary_gene_list(gene_list_path: Path) -> List[str]:
    with open(gene_list_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

