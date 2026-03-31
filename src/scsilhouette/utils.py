# src/scsilhouette/utils.py

import pandas as pd
from mygene import MyGeneInfo
from typing import List, Tuple

def get_output_prefix(organ, first_author, journal, year, cluster_header, embedding="", dataset_version_id=""):
    """Build standardized output filename prefix with embedding and vid suffix for uniqueness."""
    journal_safe = journal.replace(" ", "_") if journal else "unknown"
    cluster_header_safe = cluster_header.replace(" ", "_")
    embedding_safe = embedding.replace(" ", "_") if embedding else "unknown"
    vid_suffix = f"{dataset_version_id[-6:]}" if dataset_version_id and len(dataset_version_id) >= 6 else ""
    return f"{organ}_{first_author}_{journal_safe}_{year}_{cluster_header_safe}_{embedding_safe}_{vid_suffix}"

def map_gene_symbols_to_ensembl(symbols: List[str]) -> Tuple[List[str], pd.DataFrame]:
    mg = MyGeneInfo()
    query_results = mg.querymany(symbols, scopes="symbol", fields="ensembl.gene", species="human")

    mapped = []
    rows = []
    for res in query_results:
        symbol = res["query"]
        ensg = None
        if "notfound" not in res and "ensembl" in res:
            ensembl = res["ensembl"]
            if isinstance(ensembl, dict):
                ensg = ensembl.get("gene")
            elif isinstance(ensembl, list):
                ensg = ensembl[0].get("gene")
        rows.append({"gene_symbol": symbol, "ensembl_id": ensg})
        if ensg:
            mapped.append(ensg)

    mapping_df = pd.DataFrame(rows)
    return mapped, mapping_df

