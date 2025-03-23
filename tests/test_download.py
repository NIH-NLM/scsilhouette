from scsilhouette.download import download_h5ad
from pathlib import Path
import anndata as ad

def test_download_h5ad(tmp_path):
    url = "https://datasets.cellxgene.cziscience.com/5daeaafe-c79e-4ee4-a9f0-ddf6649adc21.h5ad"
    filepath = download_h5ad(url, dst_dir=tmp_path)
    assert filepath.exists()
    adata = ad.read_h5ad(filepath)
    assert isinstance(adata, ad.AnnData)

