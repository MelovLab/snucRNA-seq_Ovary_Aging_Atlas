"""
Utility functions for Human Ovary Atlas single-nucleus RNA-seq analysis.

Includes helper functions for scVI model training and diagnostics, marker gene
visualization, cell composition heatmaps, session information reporting, and
marker extraction from scVI and Scanpy differential expression results.
"""

# Imports
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
import sys
import importlib.metadata
import numpy as np
import os


# Run scVI on raw counts, compute latent embedding/clustering/UMAP, and return the annotated AnnData object plus trained model
def run_scvi(adata, batch_hv=None, batch_scvi='unique_sample_ID',
             cat_cov_scvi=['SenNet_ID', 'plex'], cont_cov_scvi=['pct_counts_mito'], exclude_mt_genes=True, **kwargs):

    # Prep adata for scVI and identify highly variable genes
    adata_scvi = sc.AnnData(X=adata.layers['raw_counts'].copy(), obs=adata.obs.copy(), var=adata.var.copy())
    sc.pp.highly_variable_genes(adata_scvi, flavor="seurat_v3", n_top_genes=10000, batch_key=batch_hv)
    selected_genes = list(set(adata_scvi.var.loc[adata_scvi.var['highly_variable']].index.tolist()))

    # Remove mitochondrial genes if specified
    if exclude_mt_genes:
        mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')]
        [selected_genes.remove(i) for i in mt_genes if i in selected_genes]

    # Subset adata to only include highly variable genes
    adata_scvi = adata_scvi[:, selected_genes].copy()

    # Setup anndata for scVI
    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key=batch_scvi,
                                  categorical_covariate_keys=cat_cov_scvi,
                                  continuous_covariate_keys=cont_cov_scvi)

    # Filter kwargs to only include valid scVI and VAE init parameters, setup vae
    scvi_varnames = set(scvi.model.SCVI.__init__.__code__.co_varnames) | set(scvi.module.VAE.__init__.__code__.co_varnames)
    scvi_kwargs = {k: v for k, v in kwargs.items() if k in scvi_varnames}
    vae = scvi.model.SCVI(adata_scvi, **scvi_kwargs)

    # Filter kwargs to only include valid training parameters, train vae
    train_kwargs = {k: v for k, v in kwargs.items() if k in vae.train.__code__.co_varnames}
    vae.train(**train_kwargs)

    # Get latent representation
    adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()

    # Run UMAP and Leiden clustering
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.leiden(adata_scvi)
    sc.tl.umap(adata_scvi, min_dist=0.3)

    # #Transfer metadata and UMAP/Leiden clustering to original adata
    adata_raw_scvi = adata.copy()
    adata_raw_scvi.obsm['X_scVI'] = adata_scvi.obsm['X_scVI'].copy()
    adata_raw_scvi.obsm['X_umap'] = adata_scvi.obsm['X_umap'].copy()
    adata_raw_scvi.obsp = adata_scvi.obsp.copy()
    adata_raw_scvi.uns = adata_scvi.uns.copy()
    adata_raw_scvi.obs['leiden'] = adata_scvi.obs['leiden'].copy()

    sc.pl.umap(adata_raw_scvi, color=["Age", batch_scvi], frameon=False, ncols=2)
    sc.pl.umap(adata_raw_scvi, color=['total_counts', 'pct_counts_mito'])
    results = {'data': adata_raw_scvi, 'vae': vae}

    return results


# Plot scVI training diagnostics, including reconstruction loss and ELBO across epochs.
def scvi_model_parameters(vae_model):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    axs[0].plot(vae_model.history["reconstruction_loss_train"], label="Training", color="blue")
    axs[0].plot(vae_model.history["reconstruction_loss_validation"], label="Validation", color="orange")
    axs[0].set_xlabel("Epochs")
    axs[0].set_ylabel("Reconstruction Loss")
    axs[0].legend()
    axs[0].set_title("Reconstruction Loss Over Epochs")

    axs[1].plot(vae_model.history["elbo_train"], label="Training", color="blue")
    axs[1].plot(vae_model.history["elbo_validation"], label="Validation", color="orange")
    axs[1].set_xlabel("Epochs")
    axs[1].set_ylabel("ELBO")
    axs[1].legend()
    axs[1].set_title("ELBO Over Epochs")

    plt.tight_layout()
    plt.show()


# Generate a Scanpy dotplot for a provided gene list, excluding genes not present in the dataset.
def make_dotplot_with_genelist(adata, gene_list, clustering, fighei=8):
    all_genes_filtered = [gene for gene in gene_list if gene in adata.var_names]
    missing_genes = [gene for gene in gene_list if gene not in adata.var_names]

    if missing_genes:
        print(f"Note: The following genes are missing from the dataset and will be excluded: {missing_genes}")

    sc.pl.dotplot(adata, layer='log1p_norm', var_names=all_genes_filtered, groupby=clustering, standard_scale="var", figsize=(20, fighei))


# Plot the proportional composition of one AnnData observation category within another as a heatmap.
def cell_proportions_heatmap(adata, obs1, obs2, figwid=16, fighei=8):
    contingency_table = pd.crosstab(adata.obs[obs1], adata.obs[obs2])
    proportion_table = contingency_table.div(contingency_table.sum(axis=1), axis=0) * 100

    plt.figure(figsize=(figwid, fighei))
    ax = sns.heatmap(proportion_table, cmap='coolwarm', annot=False, linewidths=0.5, linecolor='black',
                     xticklabels=True, yticklabels=True, cbar_kws={'label': 'Percentage'})
    ax.grid(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=10)
    plt.ylabel(obs1)
    plt.xlabel(obs2)
    plt.title(f'Proportion of Cells from Each {obs2} in Each {obs1}')
    plt.tight_layout()
    plt.show()


# Run scVI differential expression for group-vs-rest or pairwise group comparisons, with optional cell subsetting.
def scvi_de_wrapper(
    vae,
    adata,
    groupby_key: str,
    groups=None,
    *,
    # --- mode control ---
    mode: str = "standard",  # "standard" | "pairwise"
    group2=None,             # for pairwise OR can be list/str of clusters to merge
    # --- subsetting control ---
    subset_to=None,          # optional list/Series mask or list of obs_names to subset before DE
    # --- misc ---
    de_kwargs: dict | None = None,
):
    """
    Unified scVI DE wrapper.

    Parameters
    ----------
    vae : scvi.model.VAE-like
        Must have .adata and .differential_expression().
    adata : anndata.AnnData
        Source of obs annotations (e.g., leiden / cell types) to copy into vae.adata.obs.
    groupby_key : str
        obs column used for DE grouping.
    groups : None | str | list[str]
        Which group(s) to run as group1. If None, uses all categories of groupby_key.
    mode : "standard" | "pairwise"
        - standard: group1 vs rest (group2=None)
        - pairwise: group1 vs group2 (possibly merged)
    group2 : str | list[str] | None
        Used in pairwise mode (or ignored otherwise).
        If list, it will be merged into a single label.
    subset_to : None | list[str] | np.ndarray[bool] | pd.Series[bool]
        Subset cells before DE. If list[str], treated as obs_names.
        If boolean mask/Series, must align to vae.adata.obs_names length.
    de_kwargs : dict
        Extra kwargs passed to vae.differential_expression().
    """
    if de_kwargs is None:
        de_kwargs = {}

    # ---------- sync groupby_key from adata into vae.adata ----------
    if groupby_key not in adata.obs.columns:
        raise KeyError(f"'{groupby_key}' not found in adata.obs.")

    # Align by obs_names (safe even if different ordering)
    src = adata.obs.reindex(vae.adata.obs_names)[groupby_key]

    # If there are missing cells, fail loudly (usually indicates mismatched obs_names)
    if src.isna().any():
        n_missing = int(src.isna().sum())
        raise ValueError(
            f"After reindexing adata.obs to vae.adata.obs_names, "
            f"{n_missing} cells have NA for '{groupby_key}'. "
            f"Check that adata.obs_names and vae.adata.obs_names match."
        )

    # Write into vae
    vae.adata.obs[groupby_key] = src

    # Keep as category if source is category, otherwise keep as-is
    if pd.api.types.is_categorical_dtype(src):
        vae.adata.obs[groupby_key] = vae.adata.obs[groupby_key].astype("category")

    # ---------- helper: subset vae ----------
    def _subset_vae(_vae, obs_names):
        if hasattr(_vae, "copy_subset"):
            return _vae.copy_subset(obs_names)
        # fallback
        _vae2 = _vae
        _vae2.adata = _vae2.adata[obs_names].copy()
        return _vae2

    # ---------- optional global subset ----------
    vae_use = vae
    if subset_to is not None:
        if isinstance(subset_to, (list, tuple, pd.Index)):
            obs_names = pd.Index(subset_to)
        else:
            # assume boolean mask / Series aligned to adata.obs_names
            mask = np.asarray(subset_to)
            obs_names = vae.adata.obs_names[mask]
        vae_use = _subset_vae(vae, obs_names)

    # ---------- determine which groups to iterate ----------
    obs = vae_use.adata.obs
    if groups is None:
        if pd.api.types.is_categorical_dtype(obs[groupby_key]):
            groups_iter = list(obs[groupby_key].cat.categories)
        else:
            groups_iter = sorted(obs[groupby_key].astype(str).unique())
    elif isinstance(groups, str):
        groups_iter = [groups]
    else:
        groups_iter = list(groups)

    de_df_by_group = {}

    # ---------- helper: build merged labels for multi-group2 ----------
    def _make_merged_group_key(_vae, base_key, group1, group2_list, label2):
        tmp_key = f"{base_key}__{group1}__vs__{label2}"
        x = _vae.adata.obs[base_key].astype(str)

        _vae.adata.obs[tmp_key] = x.where(x.eq(str(group1)), other=label2)
        _vae.adata.obs[tmp_key] = _vae.adata.obs[tmp_key].where(
            _vae.adata.obs[tmp_key].isin([str(group1), label2]),
            other=np.nan
        )
        # keep only cells in group1 or merged group2 label
        keep = _vae.adata.obs[tmp_key].notna().values
        _vae2 = _subset_vae(_vae, _vae.adata.obs_names[keep])
        _vae2.adata.obs[tmp_key] = _vae2.adata.obs[tmp_key].astype("category")
        return _vae2, tmp_key

    # ---------- main loop ----------
    for g1 in groups_iter:
        g1 = str(g1)

        vae_local = vae_use
        gby = groupby_key
        g2 = None
        comparison_label = None

        if mode == "standard":
            # group1 vs rest (within current vae_local)
            g2 = None
            comparison_label = g1

        elif mode == "pairwise":
            if group2 is None:
                raise ValueError("mode='pairwise' requires group2 (str or list[str]).")

            if isinstance(group2, (list, tuple, set, pd.Index)):
                group2_list = [str(x) for x in group2]
                label2 = f"{'_'.join(group2_list)}"
                vae_local, gby = _make_merged_group_key(
                    vae_use, groupby_key, g1, group2_list, label2=label2
                )
                g2 = label2
                comparison_label = f"{g1}_vs_{label2}"
            else:
                g2 = str(group2)
                comparison_label = f"{g1}_vs_{g2}"

        else:
            raise ValueError("mode must be one of: 'standard', 'pairwise'.")

        print(
            f"[scVI-DE] mode={mode} | "
            f"groupby='{gby}' | "
            f"comparison: '{g1}' vs "
            f"{'ALL_OTHER' if g2 is None else f'({g2})'}"
            f"n_cells={vae_local.adata.n_obs}"
        )

        # ---------- run DE ----------
        de_df = vae_local.differential_expression(
            groupby=gby,
            group1=g1,
            group2=g2,
            **de_kwargs
        )

        # ---------- annotate  ----------
        de_df = de_df.copy()
        de_df["cluster"] = comparison_label
        de_df_by_group[comparison_label] = de_df

    # ---------- combine de_df ----------
    de_df_all = (
        pd.concat(de_df_by_group.values(), axis=0, ignore_index=False)
        if len(de_df_by_group) > 0
        else pd.DataFrame()
    )

    return {
        "de_df": de_df_all,
        "de_df_by_group": de_df_by_group,
    }


# Filter scVI differential expression results and extract top marker genes for each comparison.
def extract_top_scvi_markers(
    de_df_by_group: dict,
    *,
    # filtering thresholds
    min_lfc: float = 1.0,
    min_bayes: float = 2.0,
    min_nz1: float = 0.1,
    # ranking
    top_n: int = 8,
    sort_by: str = "lfc_mean",   # or "bayes_factor"
    ascending: bool = False,
    # output control
    columns: tuple = ("lfc_mean", "bayes_factor", "non_zeros_proportion1"),
    return_filtered_df: bool = False
):
    """
    Post-hoc filtering + top-marker extraction from scVI DE results.

    Behavior:
        - If len(de_df_by_group) > 1: extract top_n from the positive side only (group1-enriched).
        - If len(de_df_by_group) == 1: extract top_n from BOTH sides:
            group1-enriched (lfc_mean >  min_lfc)
            group2-enriched (lfc_mean < -min_lfc)
            returning 2*top_n genes total by default.
    """
    top_markers_gene_bayes_dict = {}
    top_marker_genes_list = []
    filtered_de_df_by_group = {}

    single_comparison = (len(de_df_by_group) == 1)

    for comparison, ddf in de_df_by_group.items():
        if ddf is None or ddf.empty:
            top_markers_gene_bayes_dict[comparison] = {}
            filtered_de_df_by_group[comparison] = ddf
            continue

        def _filter_common(xdf):
            return xdf.loc[
                (xdf["bayes_factor"] > min_bayes)
                & (xdf["non_zeros_proportion1"] > min_nz1)
            ]

        def _select_top(xdf, key_name):
            if xdf.empty:
                top_markers_gene_bayes_dict[key_name] = {}
                return
            top_df = (
                xdf.sort_values(sort_by, ascending=ascending)
                   .head(top_n)[list(columns)]
            )
            top_dict = top_df.to_dict("index")
            top_markers_gene_bayes_dict[key_name] = top_dict
            top_marker_genes_list.append(list(top_dict.keys()))

        if not single_comparison:
            # ---- original behavior: positive side only ----
            ddf_pos = ddf.loc[ddf["lfc_mean"] > min_lfc]
            ddf_pos = _filter_common(ddf_pos)
            filtered_de_df_by_group[comparison] = ddf_pos
            _select_top(ddf_pos, comparison)

        else:
            # ---- single comparison: take BOTH sides ----

            # group1-enriched (positive lfc)
            ddf_pos = ddf.loc[ddf["lfc_mean"] > min_lfc]
            ddf_pos = _filter_common(ddf_pos)

            # group2-enriched (negative lfc) -> flip sign for ranking/interpretation
            ddf_neg = ddf.loc[ddf["lfc_mean"] < -min_lfc].copy()
            ddf_neg["lfc_mean"] = -ddf_neg["lfc_mean"]
            ddf_neg = _filter_common(ddf_neg)

            # store filtered dfs
            filtered_de_df_by_group[f"{comparison}|group1"] = ddf_pos
            filtered_de_df_by_group[f"{comparison}|group2"] = ddf_neg

            # select top for each side
            _select_top(ddf_pos, f"{comparison}|group1")
            _select_top(ddf_neg, f"{comparison}|group2")

    top_marker_genes_list_flat = [g for sub in top_marker_genes_list for g in sub]

    out = {
        "top_markers_gene_bayes_dict": top_markers_gene_bayes_dict,
        "top_markers": top_marker_genes_list_flat,
    }

    if return_filtered_df:
        out["filtered_de_df_by_group"] = filtered_de_df_by_group

    return out


# Run Scanpy marker-gene testing and return filtered top markers plus the full marker table.
def get_top_scanpy_markers(
    adata,
    groupby,
    compartment,
    layer="log1p_norm",
    method="wilcoxon",
    key_added="rank_genes_groups",
    top_n=100,
    pval_adj_cutoff=0.05,
    fraction_expressing_cutoff=0.10,
    require_positive_lfc=True,
    use_raw=False,
    rank_by=("pval_adj", "logfoldchange", "score"),
    copy=True,
):
    """
    Run scanpy.tl.rank_genes_groups and return a cleaned top-marker table.

    Notes
    -----
    Default filtering:
        pval_adj <= 0.05
        fraction_expressing_in >= 0.10
        logfoldchange > 0, if require_positive_lfc=True
    """

    if groupby not in adata.obs.columns:
        raise ValueError(f"`groupby='{groupby}'` not found in adata.obs.")

    if layer is not None and layer not in adata.layers:
        raise ValueError(f"`layer='{layer}'` not found in adata.layers.")

    adata_out = adata.copy() if copy else adata
    adata_out.obs[groupby] = adata_out.obs[groupby].astype("category")
    groups = list(adata_out.obs[groupby].cat.categories)

    sc.tl.rank_genes_groups(
        adata_out,
        groupby=groupby,
        layer=layer,
        method=method,
        key_added=key_added,
        use_raw=use_raw,
        pts=True,
    )

    pts = adata_out.uns[key_added].get("pts")
    pts_rest = adata_out.uns[key_added].get("pts_rest")

    if pts is None or pts_rest is None:
        raise RuntimeError(
            "Scanpy did not return `pts`/`pts_rest`. "
            "Make sure rank_genes_groups was run with pts=True."
        )

    marker_tables = []

    for group in groups:
        df = sc.get.rank_genes_groups_df(
            adata_out,
            group=group,
            key=key_added,
        ).rename(columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "logfoldchange",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        })

        # Drop redundant Scanpy fraction columns if present.
        # These appear when pts=True and are equivalent to the clearer
        # fraction_expressing_in/out columns added below.
        df = df.drop(
            columns=[c for c in ["pct_nz_group", "pct_nz_reference"] if c in df.columns],
            errors="ignore"
        )

        df["compartment"] = compartment
        df["cell_state"] = str(group)
        df["fraction_expressing_in"] = df["gene"].map(pts[group])
        df["fraction_expressing_out"] = df["gene"].map(pts_rest[group])

        marker_tables.append(df)

    markers_all = pd.concat(marker_tables, ignore_index=True)

    keep = (
        (markers_all["pval_adj"] <= pval_adj_cutoff) &
        (markers_all["fraction_expressing_in"] >= fraction_expressing_cutoff)
    )

    if require_positive_lfc:
        keep = keep & (markers_all["logfoldchange"] > 0)

    markers_filt = markers_all.loc[keep].copy()

    ascending = []
    for col in rank_by:
        if col not in markers_filt.columns:
            raise ValueError(f"`rank_by` column '{col}' not found in marker table.")

        if col in {"pval", "pval_adj"}:
            ascending.append(True)
        else:
            ascending.append(False)

    markers_filt = markers_filt.sort_values(
        ["cell_state", *rank_by],
        ascending=[True, *ascending],
    )

    markers_filt["rank"] = (
        markers_filt.groupby("cell_state").cumcount() + 1
    )

    markers_top = markers_filt.loc[markers_filt["rank"] <= top_n].copy()

    preferred_cols = [
        "compartment",
        "cell_state",
        "rank",
        "gene",
        "score",
        "logfoldchange",
        "pval",
        "pval_adj",
        "fraction_expressing_in",
        "fraction_expressing_out",
    ]

    top_cols = [c for c in preferred_cols if c in markers_top.columns]
    top_cols += [c for c in markers_top.columns if c not in top_cols]
    markers_top = markers_top[top_cols]

    all_cols = [c for c in preferred_cols if c in markers_all.columns]
    all_cols += [c for c in markers_all.columns if c not in all_cols]
    markers_all = markers_all[all_cols]

    return markers_top, markers_all, adata_out


# Print Scanpy environment information and versions for loaded Python packages.
def get_session_info():
    # List of imported modules & versions of packages used by scanpy
    sc.settings.verbosity = 3
    sc.logging.print_header()

    # Print the version of each python module
    for module_name in sys.modules:
            try:
                version = importlib.metadata.version(module_name)
                print(f"{module_name} version: {version}")
            except importlib.metadata.PackageNotFoundError:
                # Skip modules that don't have a version available
                continue
            except Exception as e:
                # Handle any other exceptions that might occur
                print(f"Error retrieving version for {module_name}: {e}")