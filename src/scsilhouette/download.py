# src/scsilhouette/download.py

import os
import pooch
from urllib.parse import urlparse

def get_filename_from_url(url: str) -> str:
    return os.path.basename(urlparse(url).path)

def download_h5ad(url: str, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    filename = get_filename_from_url(url)

    # Let pooch manage the full path under output_dir
    file_path = pooch.retrieve(
        url=url,
        known_hash=None,
        fname=filename,
        path=output_dir,
        progressbar=True,
    )

    return file_path

