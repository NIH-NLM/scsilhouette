import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_score_distribution(df: pd.DataFrame, output_dir: str):
    plt.figure()
    sns.histplot(df["silhouette_score"], kde=True, bins=30)
    plt.title("Silhouette Score Distribution")
    plt.xlabel("Score")
    plt.ylabel("Count")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "silhouette_distribution.png"))
    plt.close()

def plot_cluster_summary(summary: pd.DataFrame, output_dir: str):
    plt.figure()
    sns.barplot(data=summary, x="label", y="mean_silhouette")
    plt.title("Mean Silhouette Score per Cluster")
    plt.xticks(rotation=90)
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "cluster_summary.png"))
    plt.close()

def plot_cluster_size_vs_score(summary: pd.DataFrame, output_dir: str):
    plt.figure()
    sns.scatterplot(
        data=summary,
        x="cluster_size",
        y="mean_silhouette",
        hue="label"
    )
    plt.title("Cluster Size vs Mean Silhouette Score")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "size_vs_score.png"))
    plt.close()

