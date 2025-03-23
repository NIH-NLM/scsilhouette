# src/scsilhouette/download.py

import os
import pooch

def download_h5ad(url: str, output_dir: str = "data") -> str:
    os.makedirs(output_dir, exist_ok=True)
    fname = os.path.join(output_dir, os.path.basename(url))

    print(f"Downloading file to: {fname}")
    file_path = pooch.retrieve(
        url=url,
        known_hash=None,
        fname=fname,
        path=output_dir,
        progressbar=True,
    )
    return file_path

