Process Overview
================

The `scsilhouette` Python package computes silhouette scores on single-cell RNA-seq data
and performs visualizations and association analysis with external metrics like F-scores.

Usage Workflow
--------------
0. **Available Commands**
   ```
   scsilhouette --help
   ```

1. **compute-silhouette** given a dataset (e.g. from cellxgene), compute the silhouette scores and generate a cluster-summary and the silhouette scores fore each cluster.  Type *--help* for the latest options.
   ```
 scsilhouette compute-silhouette --help
                                                                                                           
 Usage: scsilhouette compute-silhouette [OPTIONS]                                                          
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --h5ad-path                                            PATH  [default: None] [required]              │
│ *  --label-key                                            TEXT  [default: None] [required]              │
│ *  --embedding-key                                        TEXT  [default: None] [required]              │
│    --use-binary-genes        --no-use-binary-genes              [default: no-use-binary-genes]          │
│    --gene-list-path                                       PATH  [default: None]                         │
│    --metric                                               TEXT  [default: euclidean]                    │
│    --save-scores             --no-save-scores                   [default: no-save-scores]               │
│    --save-cluster-summary    --no-save-cluster-summary          [default: no-save-cluster-summary]      │
│    --save-annotation         --no-save-annotation               [default: no-save-annotation]           │
│    --show-annotation         --no-show-annotation               [default: no-show-annotation]           │
│    --help                                                       Show this message and exit.             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

2. **viz-summary** given the silhouette-score (output from the compute-silhouette), and the NSForest output file, F-score result, generate a summary picture of the dataset for use to gauge the dataset quality.  Type
   ```
 scsilhouette viz-summary --help       
                                                                                                           
 Usage: scsilhouette viz-summary [OPTIONS]                                                                 
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --silhouette-score-path                 PATH  [default: None] [required]                             │
│ *  --label                                 TEXT  [default: None] [required]                             │
│ *  --score-col                             TEXT  [default: None] [required]                             │
│    --fscore-path                           PATH  [default: None]                                        │
│    --mapping-path                          PATH  [default: None]                                        │
│    --show                     --no-show          [default: no-show]                                     │
│    --sort-by                               TEXT  Sort by mean|median|std [default: median]              │
│    --help                                        Show this message and exit.                            │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

3. **viz-correlation** display the correlation between any two vectors of strings from the data obtained from the compute-silhouette score
 
   ```
 scsilhouette viz-correlation --help   
                                                                                                           
 Usage: scsilhouette viz-correlation [OPTIONS]                                                             
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --cluster-summary-path                 PATH  CSV with summary metrics (e.g., mean, fscore, count)    │
│                                                 [default: None]                                         │
│                                                 [required]                                              │
│ *  --x-metric                             TEXT  X-axis metric for correlation (e.g., fscore)            │
│                                                 [default: None]                                         │
│                                                 [required]                                              │
│ *  --y-metrics                            TEXT  Comma-separated list of Y-axis metrics (e.g.,           │
│                                                 mean,median,std)                                        │
│                                                 [default: None]                                         │
│                                                 [required]                                              │
│ *  --label                                TEXT  Label column used for filenames and groupings           │
│                                                 [default: None]                                         │
│                                                 [required]                                              │
│    --show                    --no-show          [default: no-show]                                      │
│    --fscore-path                          PATH  Optional path to fscore CSV [default: None]             │
│    --mapping-path                         PATH  Optional mapping file to match cluster labels           │
│                                                 [default: None]                                         │
│    --help                                       Show this message and exit.                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```
4. **nsforest-genes** this will convert gene symbol names to ENSG for consistency with the output generated expected from cellxgene h5ad files.

   ```
 scsilhouette nsforest-genes --help
                                                                                                           
 Usage: scsilhouette nsforest-genes [OPTIONS]                                                              
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --nsforest-path        PATH  [default: None] [required]                                              │
│    --help                       Show this message and exit.                                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

5. **viz-dotplot** will construct the dotplot for all the clusters using the chosen embedding and label

   ```
scsilhouette viz-dotplot --help   
                                                                                                           
 Usage: scsilhouette viz-dotplot [OPTIONS]                                                                 
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --h5ad-path            TEXT  Path to input .h5ad file [default: None] [required]                     │
│ *  --embedding-key        TEXT  Embedding key (X_umap, X_scanvi_emb) [default: None] [required]         │
│ *  --label-key            TEXT  Label key for color [default: None] [required]                          │
│    --help                       Show this message and exit.                                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

6. **viz-distribution** will display the distributions of the mean, median, std and cell counts per cluster

   ```
 scsilhouette viz-distribution --help
                                                                                                           
 Usage: scsilhouette viz-distribution [OPTIONS]                                                            
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --cluster-summary-path PATH  CSV file from compute-silhouette with mean, median, std, count          │
│                                 [default: None]                                                         │
│                                 [required]                                                              │
│ *  --label-key            TEXT  Column name used as x-axis (e.g., ann_finest_level) [default: None]     │
│                                 [required]                                                              │
│    --help                       Show this message and exit.                                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

7. **viz-dataset-summary** will perform summary analysis and provide descriptive statistics overall for the dataset for the chosen cluster using the cluster summary files generated by compute-silhouette

   ```
scsilhouette viz-dataset-summary --help
                                                                                                           
 Usage: scsilhouette viz-dataset-summary [OPTIONS]                                                         
                                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --cluster-summary-path                 TEXT  Path to cluster_summary file created by the compute     │
│                                                 silhouette score                                        │
│                                                 [default: None]                                         │
│                                                 [required]                                              │
│ *  --label                                TEXT  Label for the clusters [default: None] [required]       │
│    --show                    --no-show          Show the plot interactively [default: no-show]          │
│    --help                                       Show this message and exit.                             │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
   ```

