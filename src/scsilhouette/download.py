#!/usr/bin/env python

import sys
import os
import pooch

def download_h5ad(url: str, fname: str = None, out_dir: str = "data") -> str:
    if fname is None:
        fname = os.path.basename(url)

    file_path = pooch.retrieve(
        url=url,
        known_hash=None,  # Set SHA256 if known for reproducibility
        fname=fname,
        path=out_dir,
        progressbar=True
    )

    print(f"Downloaded to: {file_path}")
    return file_path

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python download_h5ad.py <url> [filename]")
        sys.exit(1)

    url = sys.argv[1]
    fname = sys.argv[2] if len(sys.argv) > 2 else None
    download_h5ad(url, fname)
