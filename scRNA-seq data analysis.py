# -*- coding: utf-8 -*-
## scRNA-seq Analysis of a Human Blood-Derived Immune Cell Dataset This notebook performs a full single-cell RNA-seq analysis workflow using Scanpy: quality control (QC), doublet detection (Scrublet), normalisation, dimensionality reduction, clustering, cell-type annotation (using decoupler + PanglaoDB), and visualization.

## The goal is to: 
#(1) produce a clean, reproducible analysis pipeline **and**
#(2) Infer the likely tissue/source of the sample based on the resulting cell-type composition.

#**Installations to set up Scanpy environment.**

!pip install scanpy
!pip install anndata
!pip3 install igraph
!pip install celltypist
!pip install decoupler
!pip install fa2-modified

#1. Load data & basic checks

#Load data and standardise gene/cell names

import scanpy as sc
import anndata as ad

#Load pre-processed count matrix (cells x genes)
b_marrow_adata = sc.read_h5ad("/content/bone_marrow.h5ad")
print(b_marrow_adata)

# Exploration of dataset (shape of the anndata)
b_marrow_adata.shape

# Exploration of dataset (first five genes in the anndata)
b_marrow_adata.var.head()

# Exploration of dataset(first five cell data)
b_marrow_adata.obs.head()

# Exploration of dataset
b_marrow_adata.to_df()

#2. Gene-level QC flags (MT, RIBO, HB)
#Why: identify mitochondrial, ribosomal and hemoglobin genes for per-cell QC metrics.

##Ensure unique gene and cell names (required by Scanpy)
b_marrow_adata.var_names_make_unique()
b_marrow_adata.obs_names_make_unique()

##Flag mitochondrial, ribosomal and haemoglobin genes

###Mitochondrial genes: MT- prefix (human-style)
b_marrow_adata.var['MT'] = b_marrow_adata.var_names.str.startswith("MT-")

###Ribosomal: RPS* and RPL*
b_marrow_adata.var['RIBO'] = b_marrow_adata.var_names.str.startswith(("RPS", "RPL"))

 ###Haemoglobin: selected known genes
hb_gene_list = ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBG1', 'HBG2']
# Convert var_names to uppercase for case-insensitive comparison with hb_gene_list
sym_upper = b_marrow_adata.var_names.str.upper()
# Set 'HB' to True for genes whose uppercase names are in hb_gene_list
b_marrow_adata.var['HB'] = sym_upper.isin(hb_gene_list)

# Quick sanity check: how many genes in each category?
b_marrow_adata.var[['MT', 'RIBO', 'HB']].sum()

##Calculating QC metrics
sc.pp.calculate_qc_metrics(
    b_marrow_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)
b_marrow_adata.obs.head()

b_marrow_adata.var.head()

#3. Cell-level QC metrics & filtering
#Why: remove low-quality cells and outliers (likely doublets / dying cells) before downstream analysis.


# Basic gene/cell filters
sc.pp.filter_cells(b_marrow_adata, min_genes=200)
sc.pp.filter_genes(b_marrow_adata, min_cells=3)

# Mito threshold (tune based on plots)
b_marrow_adata = b_marrow_adata[b_marrow_adata.obs['pct_counts_MT'] < 20, :].copy()

b_marrow_adata = b_marrow_adata[b_marrow_adata.obs['pct_counts_RIBO'] < 20, :].copy()

b_marrow_adata = b_marrow_adata[b_marrow_adata.obs['pct_counts_HB'] < 20, :].copy()

#Average number of genes with at least one detected identifier in each cell
sc.pl.violin(
    b_marrow_adata,
    ["n_genes_by_counts"],
    jitter=0.4,
    multi_panel=False,
)

#Total number of molecules (UMI) detected in cell.
sc.pl.violin(
    b_marrow_adata,
    ["total_counts"],
    jitter=0.4,
    multi_panel=False,
)

#Visualization of mt genes
sc.pl.scatter(b_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_MT")

#visualization of ribosomal genes
sc.pl.scatter(b_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_RIBO")

#Visualization of HB genes
sc.pl.scatter(b_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_HB")

b_marrow_adata.obs.columns

# 4. Doublet detection with Scrublet
#Why: identify droplets that likely contain more than one cell and could create artifactual clusters.

!pip install scrublet

import scrublet as scr
import scipy.sparse as sp

sc.pp.scrublet(b_marrow_adata)

# ---------- Scrublet doublet detection ----------

# Use raw counts
if "counts" in b_marrow_adata.layers:
    counts_matrix = b_marrow_adata.layers["counts"]
else:
    counts_matrix = b_marrow_adata.X

# Convert to dense scrublet
if sp.issparse(counts_matrix):
    counts_matrix = counts_matrix.toarray()

#Initialize Scrublet with a reasonable expected doublet rate (~5-10%)
scrub = scr.Scrublet(
    counts_matrix,
    expected_doublet_rate=0.06  # 6% prior; fine for coursework unless told otherwise
)

# Compute pre-cell doublet scores and binary calls
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Store outputs in pre-cell metadata
b_marrow_adata.obs["doublet_score"] = doublet_scores
b_marrow_adata.obs["predicted_doublet"] = predicted_doublets

# Sanity check
b_marrow_adata.obs[["doublet_score", "predicted_doublet"]].head()

# Compare QC metrics between singlets and predicted doublets
b_marrow_adata.obs['predicted_doublet'] = b_marrow_adata.obs['predicted_doublet'].astype('category')
sc.pl.violin(
    b_marrow_adata,
    ["n_genes_by_counts", "total_counts"],
    groupby="predicted_doublet",
    multi_panel=True
)

# ---------- Filter out predicted doublets ----------

before = b_marrow_adata.n_obs
b_marrow_adata = b_marrow_adata[b_marrow_adata.obs["predicted_doublet"] == False].copy()
after = b_marrow_adata.n_obs

print(f"Removed {before - after} predicted doublets; remaining {after} cells.")

# Compare QC metrics between singlets and predicted doublets
b_marrow_adata.obs['predicted_doublet'] = b_marrow_adata.obs['predicted_doublet'].astype('category')
sc.pl.violin(
    b_marrow_adata,
    ["n_genes_by_counts", "total_counts"],
    groupby="predicted_doublet",
    multi_panel=True
)

# 5. Normalization, log-transform and highly variable genes
#Why: stabilize variance and focus on informative genes for downstream dimensionality reduction.

#Save a copy of the bone marrow data
b_marrow_adata.layers["counts"] = b_marrow_adata.X.copy()

#Normalizing to median total counts
sc.pp.normalize_total(b_marrow_adata)
#Log of the data
sc.pp.log1p(b_marrow_adata)

#Selecting top 1000 most variable genes
sc.pp.highly_variable_genes(b_marrow_adata, n_top_genes=1000)

sc.pl.highly_variable_genes(b_marrow_adata)

#6. PCA, neighborhood graph, UMAP and Leiden clustering
#Why: embed cells into a low-dimensional manifold and identify transcriptionally similar groups.

sc.tl.pca(b_marrow_adata)
sc.pl.pca_variance_ratio(b_marrow_adata, n_pcs=10, log=False)

sc.pl.pca(
    b_marrow_adata,
    color=["pct_counts_MT"]
)

sc.pp.neighbors(b_marrow_adata)
sc.tl.umap(b_marrow_adata)

sc.pl.umap(
    b_marrow_adata,
    color=["pct_counts_RIBO"],
    size=8,
)

sc.tl.leiden(b_marrow_adata, flavor="igraph", n_iterations=2)

sc.pl.umap(
    b_marrow_adata,
    color=["leiden"],
    size=8,
)

sc.pl.umap(
    b_marrow_adata,
    color=["leiden"],
    wspace=0.5,
    size=3,
    ncols = 1
)

b_marrow_adata.obs.columns

sc.pl.umap(b_marrow_adata, color=["leiden", "doublet_score", "predicted_doublet"])

sc.tl.leiden(b_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_02", resolution=0.02)
sc.tl.leiden(b_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(b_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res2", resolution=2)

sc.pl.umap(
    b_marrow_adata,
    color=["leiden_res0_02"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1
)

sc.pl.umap(
    b_marrow_adata,
    color=["leiden_res0_5"],
    wspace=0.5,
    size=15,
    ncols = 1,
    legend_loc="on data"
)

sc.pl.umap(
    b_marrow_adata,
    color=["leiden_res2"],
    wspace=0.5,
    size=15,
    ncols = 1,
    legend_loc="on data"
)

#7. Cell-type scoring with decoupler (PanglaoDB)
#Why: use reference marker signatures to quantify how “neutrophil-like”, “B cell-like”, etc. each cell is.


import decoupler as dc

# Query Omnipath and get CellMarker
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Keep canonical cell type markers alone
markers = markers[markers["canonical_marker"]]

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

# Format because dc only accepts cell_type and genesymbol

markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]

markers.tail()

# This line was causing an error and is not a valid operation for a pandas DataFrame.
# If you intended to remove the 'markers' variable from memory, you could use 'del markers'.
# markers.free()

import pandas as pd
import decoupler as dc

# Query Omnipath and get CellMarker
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Keep canonical cell type markers alone
markers = markers[markers["canonical_marker"]]

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

# Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]

# Convert 'feature_name' column to a plain Index of strings before assigning to var_names
b_marrow_adata.var_names = pd.Index(b_marrow_adata.var['feature_name'].astype(str))
b_marrow_adata.var_names_make_unique()

dc.mt.ulm(data=b_marrow_adata,
          net=markers,
          tmin = 3)

# retrieve the score for each cell type
score = dc.pp.get_obsm(b_marrow_adata, key="score_ulm")

sc.pl.dotplot(
    score, # Use the 'score' AnnData object directly
    var_names=score.var_names.tolist(), # List of cell types (features)
    groupby='leiden_res0_5',
    cmap='RdBu_r', # A common colourmap for scores
    dendrogram=True,
    standard_scale='var', # Scale scores for better visualisation
    title='Cell Type Scores by Leiden Clusters'
)

#preview the data
b_marrow_adata.obsm["score_ulm"].head()

b_marrow_adata.obsm["score_ulm"].columns

sc.pl.umap(score, color=["NK cells", "leiden_res2"], cmap="RdBu_r")

sc.pl.umap(score, color=["NK cells", "leiden_res0_5"], cmap="RdBu_r")

#8. Cluster annotation and visualisation
#Why: assign human-readable cell-type labels to clusters and interpret the cellular composition.

# Get cell-type scores as DataFrame
score = dc.pp.get_obsm(b_marrow_adata, key="score_ulm")

#Rank signatures per cluster for your main clustering
rank_df = dc.tl.rankby_group(
    score,
    groupby="leiden_res2",  # use the same clustering
    reference="rest",
    method="t-test_overestim_var",
)

# Keep only positive stats (enriched vs rest)
rank_df = rank_df[rank_df["stat"] > 0]

#Take top signature per cluster
cluster_auto = (
    rank_df
    .sort_values(["group", "stat"], ascending=[True, False])
    .groupby("group")
    .head(1)
    .set_index("group")["name"]   # 'name' = signature (e.g. "Neutrophils")
    .to_dict()
)

cluster_auto

b_marrow_adata.obs["cell_type"] = (
    b_marrow_adata.obs["leiden_res2"].map(cluster_auto)
).astype("category")

sc.pl.umap(b_marrow_adata, color=["cell_type"], legend_loc="on data")

# Get cell-type scores as DataFrame
score = dc.pp.get_obsm(b_marrow_adata, key="score_ulm")

# Rank signatures per cluster for your main clustering
rank_df = dc.tl.rankby_group(
    score,
    groupby="leiden_res0_02",  # use the same clustering
    reference="rest",
    method="t-test_overestim_var",
)

# Keep only positive stats (enriched vs rest)
rank_df = rank_df[rank_df["stat"] > 0]

# Take top signature per cluster
cluster_auto = (
    rank_df
    .sort_values(["group", "stat"], ascending=[True, False])
    .groupby("group")
    .head(1)
    .set_index("group")["name"]   # 'name' = signature (e.g. "Neutrophils")
    .to_dict()
)

cluster_auto

b_marrow_adata.obs["cell_type_auto"] = (
    b_marrow_adata.obs["leiden_res0_02"].map(cluster_auto)
).astype("category")

sc.pl.umap(b_marrow_adata, color=["cell_type_auto"], legend_loc="on data")

# Get cell-type scores as DataFrame
score = dc.pp.get_obsm(b_marrow_adata, key="score_ulm")

# Rank signatures per cluster for your main clustering
rank_df = dc.tl.rankby_group(
    score,
    groupby="leiden_res0_5",  # use the same clustering
    reference="rest",
    method="t-test_overestim_var",
)

# Keep only positive stats (enriched vs rest)
rank_df = rank_df[rank_df["stat"] > 0]

# Take top signature per cluster
cluster_auto = (
    rank_df
    .sort_values(["group", "stat"], ascending=[True, False])
    .groupby("group")
    .head(1)
    .set_index("group")["name"]   # 'name' = signature (e.g. "Neutrophils")
    .to_dict()
)

cluster_auto

b_marrow_adata.obs["cell_type"] = (
    b_marrow_adata.obs["leiden_res0_5"].map(cluster_auto)
).astype("category")

sc.pl.umap(b_marrow_adata, color=["cell_type"], legend_loc="on data")

import seaborn as sns

sc.pl.violin(score, keys=["Neutrophils"], groupby="leiden_res0_5", rotation=90)

sc.pl.violin(score, keys=["NK cells"], groupby="leiden_res0_5", rotation=90)

sc.pl.violin(score, keys=["Plasma cells"], groupby="leiden_res0_5", rotation=90)

sc.pl.violin(score, keys=["B cells naive"], groupby="leiden_res0_5", rotation=90)

sc.pl.violin(score, keys=["Nuocytes"], groupby="leiden_res0_5", rotation=90)

#rank genes
b_marrow_adata_rank = dc.tl.rankby_group(score, groupby="leiden_res0_02", reference="rest", method="t-test_overestim_var")
b_marrow_adata_rank = b_marrow_adata_rank[b_marrow_adata_rank["stat"] > 0]

#rank genes
b_marrow_adata_rank = dc.tl.rankby_group(score, groupby="leiden_res0_5", reference="rest", method="t-test_overestim_var")
b_marrow_adata_rank = b_marrow_adata_rank[b_marrow_adata_rank["stat"] > 0]
b_marrow_adata_rank.head()

#rank genes
b_marrow_adata_rank = dc.tl.rankby_group(score, groupby="leiden_res2", reference="rest", method="t-test_overestim_var")
b_marrow_adata_rank = b_marrow_adata_rank[b_marrow_adata_rank["stat"] > 0]
b_marrow_adata_rank.head()

cluster_annotations = b_marrow_adata_rank[b_marrow_adata_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()

cluster_annotations

b_marrow_adata.obs['cell_type'] = b_marrow_adata.obs['leiden_res0_02'].map(cluster_annotations)

b_marrow_adata.obs['cell_type'] = b_marrow_adata.obs['leiden_res0_5'].map(cluster_annotations)

b_marrow_adata.obs['cell_type'] = b_marrow_adata.obs['leiden_res2'].map(cluster_annotations)

# Example of how to subset for multiple genes in the 'source' column
available_genes = set(b_marrow_adata.var_names)

neutro_markers = markers[markers['source'].isin(['Neutrophils'])]['target']
neutro_markers = neutro_markers[neutro_markers.isin(available_genes)]

macro_markers = markers[markers['source'].isin(['Macrophages'])]['target']
macro_markers = macro_markers[macro_markers.isin(available_genes)]

den_cells_markers = markers[markers['source'].isin(['Dendritic cells'])]['target']
den_cells_markers = den_cells_markers[den_cells_markers.isin(available_genes)]

kup_cells_markers = markers[markers['source'].isin(['Kupffer cells'])]['target']
kup_cells_markers = kup_cells_markers[kup_cells_markers.isin(available_genes)]

micro_markers = markers[markers['source'].isin(['Microglia'])]['target']
micro_markers = micro_markers[micro_markers.isin(available_genes)]

marker_genes_dict = {
    "Neutrophils":neutro_markers.head().tolist(),
    "Macrophages": macro_markers.head().tolist(),
    "Dendritic cells": den_cells_markers.head().tolist(),
     "Kupffer cells": kup_cells_markers.head().tolist(),
     "Microglia": micro_markers.head().tolist()
}

sc.tl.dendrogram

sc.tl.dendrogram(b_marrow_adata, groupby="leiden_res0_5")
sc.pl.matrixplot(
    b_marrow_adata,
    marker_genes_dict,
    "leiden_res0_5",
    dendrogram=True,
    cmap="Greens",
    use_raw=False
)

sc.tl.dendrogram(b_marrow_adata, groupby="leiden_res0_02")
sc.pl.matrixplot(
    b_marrow_adata,
    marker_genes_dict,
    "leiden_res0_02",
    dendrogram=True,
    cmap="Blues",
    use_raw=False
)

sc.tl.dendrogram(b_marrow_adata, groupby="leiden_res2")
sc.pl.matrixplot(
    b_marrow_adata,
    marker_genes_dict,
    "leiden_res2",
    dendrogram=True,
    cmap="Reds",
    use_raw=False
)

sc.pl.tracksplot(b_marrow_adata, marker_genes_dict, groupby="leiden_res2", dendrogram=False, use_raw=False)

b_marrow_adata.var.index = b_marrow_adata.var.index.astype(str)
b_marrow_adata.var_names_make_unique()
sc.pl.tracksplot(b_marrow_adata, marker_genes_dict, groupby="leiden_res0_5", dendrogram=False, use_raw=False)

b_marrow_adata.var.index = b_marrow_adata.var.index.astype(str)
b_marrow_adata.var_names_make_unique()
sc.pl.tracksplot(b_marrow_adata, marker_genes_dict, groupby="leiden_res0_02", dendrogram=False, use_raw=False)

(
    b_marrow_adata.obs["cell_type"]
    .value_counts(normalize=True)
    .sort_values(ascending=False)
)
