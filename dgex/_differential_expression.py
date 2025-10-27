# -- import packages: ---------------------------------------------------------
import logging
import numpy as np
import pandas as pd
import scipy.stats
import scipy.sparse
import statsmodels.sandbox.stats.multicomp
import statsmodels.stats.multitest

# -- configure logger: --------------------------------------------------------
logger = logging.getLogger(__name__)

# -- operational cls: ---------------------------------------------------------
class DifferentialExpression:

    """
    Source (in part): https://github.com/AllonKleinLab/klunctions/blob/77d2571e5083faed53b3950ffffc9dd07f60a9e5/sam/Analysis/scBasics/helper_functions.py#L1146-L1172
    """

    def __init__(self, adata, use_key=None) -> None:

        self._adata = adata
        self._obs_df = self._adata.obs.copy()
        self._var_df = self._adata.var.copy()
        self._var_idx = self._var_df.index
        self._use_key = use_key
        self._msg = "[DifferentialExpression] "

    def _group_data(
        self,
        groupby,
        group1,
        group2,
    ):

        self._grouped = self._obs_df.groupby(groupby, observed=False)

        self.group1 = self._grouped.get_group(group1)  # -> pd.DataFrame
        self.group2 = self._grouped.get_group(group2)  # -> pd.DataFrame
        self._group1_idx = self.group1.index
        self._group2_idx = self.group2.index

    def _filter_null_expression(self, idx):

        group_counts = (self._adata[idx].X > 0).sum(0) # .toarray()
        if scipy.sparse.issparse(self._adata.X) or isinstance(group_counts, np.matrix):
            group_counts = np.asarray(group_counts).flatten()
        return np.where(group_counts > 0)[0]

    def _filter_on_min_fraction_expressed(self, min_frac_expr):

        X_raw_group1 = self._adata[self._group1_idx].X.copy()
        X_raw_group2 = self._adata[self._group2_idx].X.copy()

        group1_expr = (X_raw_group1 > 0).sum(0) / len(self._group1_idx)
        group2_expr = (X_raw_group2 > 0).sum(0) / len(self._group2_idx)

        min_expr_filter_idx = (group1_expr > min_frac_expr) | (
            group2_expr > min_frac_expr
        )

        return self._var_idx[min_expr_filter_idx]

    def _filter_genes(self, min_frac_expr):

        self._group1_null_expr_mask = self._filter_null_expression(self._group1_idx)
        self._group2_null_expr_mask = self._filter_null_expression(self._group2_idx)

        self._filtered_var_idx = self._var_idx[
            np.unique(np.append(self._group1_null_expr_mask, self._group2_null_expr_mask))
        ]
        print(f"{self._msg}Filtering null-expression genes.\n\t   Genes remaining: {self.n_genes}")
        if min_frac_expr > 0:
            self._filtered_var_idx = self._filter_on_min_fraction_expressed(
                min_frac_expr
            )
            print(
                f"{self._msg}Filtering on minimum fraction of detection: {min_frac_expr*100}%.\n\t   Genes remaining: {self.n_genes}"
            )

        self._filtered_adata = self._adata[:, self._filtered_var_idx]

        self._group1_adata = self._filtered_adata[self._group1_idx]
        self._group2_adata = self._filtered_adata[self._group2_idx]

    def _register_counts(self, group_adata):

        if isinstance(self._use_key, type(None)):
            if scipy.sparse.issparse(group_adata.X):
                return group_adata.X.toarray()
            return group_adata.X

        if self._use_key in group_adata.layers:
            counts_matrix = group_adata.layers[self._use_key]
            if scipy.sparse.issparse(counts_matrix):
                return counts_matrix.toarray()
            return counts_matrix
        else:
            raise KeyError(f"Key '{self._use_key}' not found in AnnData layers.")

    @property
    def n_genes(self):
        return len(self._filtered_var_idx)

    def _single_gene_t_test(self, X1, X2):
        return scipy.stats.ttest_ind(X1, X2).pvalue

    def _adjust_pvals(self):
        """adjust p-values for multiple testing using the Benjamini-Hochberg procedure"""
        self._ADJUSTED = {}
        (
            self._ADJUSTED["reject"],
            self._ADJUSTED["corrected"],
            self._ADJUSTED["alphacSidak"],
            self._ADJUSTED["alphacBonf"],
        ) = statsmodels.stats.multitest.multipletests(self.p_values, method="fdr_bh")
        self.adjusted_pvals = self._ADJUSTED["corrected"]

    def _t_test(self):
        """run t-test loop"""
        self.p_values = np.empty(self.n_genes, dtype=float)

        for gene_idx in range(self.n_genes): # tqdm():
            X1 = self._group_1_counts[:, gene_idx]
            X2 = self._group_2_counts[:, gene_idx]
            self.p_values[gene_idx] = self._single_gene_t_test(X1, X2)

        self._adjust_pvals()

    def _return_de_genes(self, pval_cutoff=0.05):
        """initialize differential expression DataFrame with t-test results."""

        _de_df = pd.DataFrame(index=self._filtered_var_idx)
        _de_df["pval"] = self.p_values
        _de_df["adj_pval"] = self.adjusted_pvals
        _de_df["FDR_BH_passing"] = self.adjusted_pvals < pval_cutoff
        # _de_df = _de_df.reset_index().rename({"index": "cell_idx"}, axis=1) # Removed this line
        # _de_df = _de_df.sort_values("adj_pval").reset_index(drop=True) # Sorting will be done later

        return _de_df

    def _compute_group_means_ratio(self, pseudocount=0):

        self._group1_mean = self._group_1_counts.mean(axis=0) + pseudocount
        self._group2_mean = self._group_2_counts.mean(axis=0) + pseudocount
        self._log_ratio = np.log2(self._group1_mean / self._group2_mean)

        if hasattr(self, "_de_df"):
            self._de_df["group1_mean"] = self._group1_mean
            self._de_df["group2_mean"] = self._group2_mean
            self._de_df["log_ratio"] = self._log_ratio

    def _rank_sum_test(self):

        X1 = self._group_1_counts
        X2 = self._group_2_counts

        self._ranksum_pvals = np.empty(self.n_genes, dtype=float)
        for n, gene in enumerate(self._filtered_var_idx):
            pv = scipy.stats.ranksums(X1[:, n], X2[:, n])[1]
            self._ranksum_pvals[n] = pv

        ranksum_multiple_testing = statsmodels.sandbox.stats.multicomp.multipletests(
                    self._ranksum_pvals,
                    alpha=0.05,
                    method="fdr_bh",
                )
        self._adj_ranksum_pvals = ranksum_multiple_testing[1]
        self._adj_ranksum_fdr_passing = ranksum_multiple_testing[0]

        if hasattr(self, "_de_df"):            
            self._de_df['ranksum_pvals'] = self._ranksum_pvals
            self._de_df['ranksum_adj_pvals'] = self._adj_ranksum_pvals
            self._de_df['ranksum_FDR_BH_passing'] = self._adj_ranksum_fdr_passing

    @property
    def de_df(self) -> pd.DataFrame:
        """
        Return the differential expression dataframe, filtered for genes passing all FDR checks,
        and sorted by ranksum adjusted p-values.
        Returns an empty DataFrame if DE analysis has not been run.
        """
        if not hasattr(self, "_de_df") or self._de_df is None:
            logger.warning("Differential expression analysis has not been run. Returning empty DataFrame.")
            return pd.DataFrame()
        
        # Ensure all FDR columns exist before attempting to filter
        fdr_cols = [col for col in self._de_df.columns if "FDR" in col]
        if not fdr_cols:
            logger.warning("No FDR columns found in _de_df. Returning unfiltered DataFrame.")
            return self._de_df.copy()
            
        passing_all_fdr = self._de_df[fdr_cols].all(axis=1)
        filtered_df = self._de_df[passing_all_fdr].copy()
        
        sort_col = "ranksum_adj_pvals" if "ranksum_adj_pvals" in filtered_df.columns else "adj_pval"
        if sort_col in filtered_df.columns:
            return filtered_df.sort_values(sort_col)
        return filtered_df

    def __call__(
        self,
        groupby,
        group1,
        group2,
        pval_cutoff=0.05,
        min_frac_expr=0.05,
        pseudocount=1,
    ):

        self._group_data(groupby, group1, group2)
        self._filter_genes(min_frac_expr)

        self._group_1_counts = self._register_counts(self._group1_adata)
        self._group_2_counts = self._register_counts(self._group2_adata)

        self._t_test()
        self._de_df = self._return_de_genes(pval_cutoff)
        self._compute_group_means_ratio(pseudocount=pseudocount)

        self._rank_sum_test()

        n_t_test_fdr = str(self._de_df["FDR_BH_passing"].sum())
        n_ranksum_fdr = str(self._de_df["ranksum_FDR_BH_passing"].sum())

        logger.info(f"Identified {n_t_test_fdr} Benjamini-Hochberg correction-passing genes from t-test.")
        logger.info(f"Identified {n_ranksum_fdr} Benjamini-Hochberg correction-passing genes from Wilcoxon ranksum test.")

        # Sort the final _de_df
        sort_col = "ranksum_adj_pvals" if "ranksum_adj_pvals" in self._de_df.columns else "adj_pval"
        if sort_col in self._de_df.columns:
            self._de_df = self._de_df.sort_values(sort_col)
        else:
            # Fallback if neither column is present, though adj_pval should always be.
            self._de_df = self._de_df.sort_values("pval") 

        return self._de_df
