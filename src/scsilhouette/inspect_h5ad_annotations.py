# file: inspect_h5ad_annotations.py

import scanpy as sc
import sys

def inspect_h5ad(file_path: str):
    adata = sc.read_h5ad(file_path)

    print(f"\nLoaded file: {file_path}")
    print(f"Shape: {adata.shape} (cells x genes)")

    print("\nCell Annotations (.obs):")
    print(adata.obs.columns.tolist())

    print("\nGene Annotations (.var):")
    print(adata.var.columns.tolist())

    print("\nUnstructured Metadata (.uns):")
    print(list(adata.uns.keys()))

    print("\nExample values from .obs:")
    for col in adata.obs.columns[:5]:  # limit output
        print(f"{col}: {adata.obs[col].unique()[:5]}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_h5ad_annotations.py <file.h5ad>")
    else:
        inspect_h5ad(sys.argv[1])
