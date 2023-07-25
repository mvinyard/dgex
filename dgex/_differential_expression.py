
# -- import packages: ----------------------------------------------------------
from statsmodels.stats.multitest import multipletests
import statsmodels.sandbox.stats.multicomp
from licorice_font import font_format
from tqdm.notebook import tqdm
import scipy.stats
import numpy as np


# -- main operating class: -----------------------------------------------------
class DifferentialExpression:
    
    """
    Source (in part): https://github.com/AllonKleinLab/klunctions/blob/77d2571e5083faed53b3950ffffc9dd07f60a9e5/sam/Analysis/scBasics/helper_functions.py#L1146-L1172
    """
    
    def __init__(self, adata, use_key=None):

        self._adata = adata
        self._obs_df = self._adata.obs.copy()
        self._var_df = self._adata.var.copy()
        self._var_idx = self._var_df.index
        self._use_key = use_key
        msg = font_format("INFO", ['GREEN'])
        self._msg = f"- [{msg}] | "

    def _group_data(
        self,
        groupby,
        group1,
        group2,
    ):

        self._grouped = self._obs_df.groupby(groupby)

        self.group1 = self._grouped.get_group(group1)  # -> pd.DataFrame
        self.group2 = self._grouped.get_group(group2)  # -> pd.DataFrame
        self._group1_idx = self.group1.index
        self._group2_idx = self.group2.index

    def _filter_null_expression(self, idx):

        group_counts = (self._adata[idx].X > 0).sum(0).toarray()
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
            return group_adata.X

        if self._use_key in group_adata.layers:
            return group_adata.layers[self._use_key].toarray()

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
        ) = multipletests(self.p_values, method="fdr_bh")
        self.adjusted_pvals = self._ADJUSTED["corrected"]

    def _t_test(self):
        """run t-test loop"""
        self.p_values = []

        for gene_idx in tqdm(range(self.n_genes)):
            X1 = self._group_1_counts[:, gene_idx]
            X2 = self._group_2_counts[:, gene_idx]
            self.p_values.append(self._single_gene_t_test(X1, X2))

        self.p_values = np.array(self.p_values)
        self._adjust_pvals()
        
    def _return_de_genes(self, pval_cutoff=0.05):
        """identify differentially expressed genes with adjusted p-values"""

        de_adata = self._adata[:, self._filtered_var_idx]
        _de_df = de_adata.var["gene_ids"].to_frame()
        _de_df["pval"] = self.p_values
        _de_df["adj_pval"] = self.adjusted_pvals
        _de_df["FDR_BH_passing"] = self.adjusted_pvals < pval_cutoff
        _de_df = _de_df.reset_index().rename({"index": "cell_idx"}, axis=1)
        _de_df = _de_df.sort_values("adj_pval").reset_index(drop=True)

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

        self._ranksum_pvals = []
        for n, gene in enumerate(self._filtered_var_idx):
            pv = scipy.stats.ranksums(X1[:, n], X2[:, n])[1]
            self._ranksum_pvals.append(pv)
            
        self._ranksum_pvals = np.array(self._ranksum_pvals)
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
    def de_df():
        self._de_df.loc[self._de_df.filter(regex="FDR").all(1)].sort_values("ranksum_adj_pvals")

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
        
        n_t_test_fdr = font_format(str(self._de_df["FDR_BH_passing"].sum()), ['BOLD'])
        n_ranksum_fdr = font_format(str(self._de_df["ranksum_FDR_BH_passing"].sum()), ['BOLD'])
        
        print(f"{self._msg}Identified {n_t_test_fdr} Benjamini-Hochberg correction-passing genes from t-test.")
        print(f"{self._msg}Identified {n_ranksum_fdr} Benjamini-Hochberg correction-passing genes from Wilcoxon ranksum test.")
        
        return self._de_df
