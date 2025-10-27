# dgex
Differential Gene Expression

## `dgex` Usage

The `dgex` module provides a `DifferentialExpression` class for performing differential gene expression analysis on AnnData objects.

### Quick Start

1.  **Import and Initialize:**
    Import the class and initialize it with your AnnData object. You can optionally specify a layer to use for expression counts via `use_key`.

    ```python
    from dgex import DifferentialExpression
    import anndata

    # Load your AnnData object (example)
    adata = anndata.read_h5ad("your_data.h5ad")

    # Initialize the DE analysis object
    de_analyzer = DifferentialExpression(adata)
    # Or, to use a specific layer (e.g., 'counts'):
    # de_analyzer = DifferentialExpression(adata, use_key="counts")
    ```

2.  **Run Differential Expression Analysis:**
    Call the instance of the `DifferentialExpression` class, specifying the `groupby` column in `adata.obs` that defines your groups, and the names of the two groups (`group1`, `group2`) you want to compare.

    You can also set parameters like:
    *   `pval_cutoff`: The adjusted p-value threshold for determining significance (default: 0.05).
    *   `min_frac_expr`: The minimum fraction of cells in at least one group where a gene must be detected to be included in the analysis (default: 0.05).
    *   `pseudocount`: A small value added to counts before calculating log fold changes to avoid division by zero (default: 1).

    ```python
    # Perform DE analysis between 'ClusterA' and 'ClusterB' based on 'louvain' clustering
    de_results_df = de_analyzer(
        groupby="louvain",
        group1="ClusterA",
        group2="ClusterB",
        pval_cutoff=0.05,
        min_frac_expr=0.1
    )

    # The returned de_results_df contains statistics like p-values, adjusted p-values,
    # log2 fold change, and mean expression for each group.
    print(de_results_df.head())
    ```

3.  **Accessing Filtered Results:**
    The `de_df` property provides a view of the results, filtered for genes that pass all FDR significance thresholds (both t-test and rank-sum test, if applicable), and sorted by `ranksum_adj_pvals` (or `adj_pval` as a fallback).

    ```python
    significant_de_genes_df = de_analyzer.de_df
    print(significant_de_genes_df.head())
    ```

### Output `DataFrame`

The resulting DataFrame typically includes:
*   `pval`: Raw p-value from the t-test.
*   `adj_pval`: Benjamini-Hochberg adjusted p-value from the t-test.
*   `FDR_BH_passing`: Boolean indicating if the gene passed the t-test FDR threshold.
*   `group1_mean`: Mean expression in group 1.
*   `group2_mean`: Mean expression in group 2.
*   `log_ratio`: Log2 fold change (group1_mean / group2_mean).
*   `ranksum_pvals`: Raw p-value from the Wilcoxon rank-sum test.
*   `ranksum_adj_pvals`: Benjamini-Hochberg adjusted p-value from the rank-sum test.
*   `ranksum_FDR_BH_passing`: Boolean indicating if the gene passed the rank-sum test FDR threshold.

The index of the `DataFrame` will be the gene identifiers.
