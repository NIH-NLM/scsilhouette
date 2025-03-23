import argparse
import os
import anndata as ad
from scsilhouette.compute import compute_silhouette_score
from scsilhouette.viz import (
    plot_score_distribution,
    plot_cluster_summary,
    plot_cluster_size_vs_score
)
from scsilhouette.download import download_h5ad

def main():
    parser = argparse.ArgumentParser(description="Compute silhouette scores from an .h5ad file")

    parser.add_argument("--input", "-i", help="Path to input .h5ad file")
    parser.add_argument("--from-url", help="Download .h5ad from URL before running")

    parser.add_argument("--label-field", "-l", required=True, help="Field in .obs for cluster labels")
    parser.add_argument("--embedding-key", "-e", default="X_pca", help="Key in .obsm for embedding (default: X_pca)")
    parser.add_argument("--metadata-fields", "-m", default="", help="Comma-separated .obs fields to include in output")
    parser.add_argument("--output-dir", "-o", default=".", help="Directory to write output files")
    parser.add_argument("--no-plot", action="store_true", help="Disable plot generation")

    args = parser.parse_args()

    if not args.input and not args.from_url:
        raise ValueError("You must provide --input or --from-url")

    if args.from_url:
        print(f"[↓] Downloading .h5ad from {args.from_url}")
        args.input = download_h5ad(args.from_url, args.output_dir)
        print(f"[✓] Downloaded to: {args.input}")

    print(f"[🔍] Loading .h5ad file: {args.input}")
    adata = ad.read_h5ad(args.input)

    metadata_fields = [f.strip() for f in args.metadata_fields.split(",") if f.strip()]

    result = compute_silhouette_score(
        adata,
        label_field=args.label_field,
        embedding_key=args.embedding_key,
        metadata_fields=metadata_fields
    )

    os.makedirs(args.output_dir, exist_ok=True)

    out_scores = os.path.join(args.output_dir, "silhouette_output.csv")
    out_summary = os.path.join(args.output_dir, "cluster_summary.csv")

    result["cell_scores"].to_csv(out_scores, index=False)
    result["cluster_summary"].to_csv(out_summary, index=False)

    print(f"[✓] Saved: {out_scores}")
    print(f"[✓] Saved: {out_summary}")
    print(f"[✓] Mean silhouette score: {result['mean_score']:.4f}")

    if not args.no_plot:
        plot_score_distribution(result["cell_scores"], args.output_dir)
        plot_cluster_summary(result["cluster_summary"], args.output_dir)
        plot_cluster_size_vs_score(result["cluster_summary"], args.output_dir)

if __name__ == "__main__":
    main()

