# src/scsilhouette/utils.py

import pandas as pd
from mygene import MyGeneInfo
from typing import List, Tuple

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

