# scRNA-seq Analysis of a Human Blood-Derived Immune Cell Dataset

## 1. Overview

This project implements an end-to-end single-cell RNA-seq (scRNA-seq) analysis workflow in Python using Scanpy.  
The analysis starts from a count matrix stored in an `.h5ad` file (`bone_marrow.h5ad`) and performs:

- Quality control (QC) and doublet detection
- Normalisation and log-transformation
- Highly variable gene (HVG) selection
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering (Leiden)
- Cell-type scoring using decoupler and PanglaoDB marker signatures
- Manual cluster annotation and biological interpretation

Although the file is named `bone_marrow.h5ad`, the final cell-type composition is more consistent with a **blood-derived immune cell dataset (PBMC-like / peripheral blood leukocytes)** than classical bone marrow. The analysis therefore serves two purposes:

1. Provide a clean, reusable scRNA-seq analysis pipeline.
2. Infer and justify the underlying tissue/source from the resulting cellular landscape.

---

## 2. Environment and Installation

### 2.1. Requirements

Core Python packages:

- `scanpy`
- `anndata`
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `seaborn`
- `decoupler`
- `scrublet`
- `igraph`            (backend for Leiden clustering)
- `celltypist`        (optional, alternative annotation)
- `fa2-modified`      (optional, for graph layouts)



### 2.2. Setup

#### Create and activate a virtual environment (optional but recommended)
    python -m venv .venv
    source .venv/bin/activate    # Linux/macOS
    # .venv\Scripts\activate      # Windows

#### Install dependencies
pip install -r requirements.txt
Once installed, open and run the notebook:
scRNA_seq_data_analysis_Lavinia2.ipynb

----

## 3. Data

### 3.1. Input format
The workflow expects a single AnnData object in .h5ad format:

     bone_marrow.h5ad

Content:
    
    .X: raw counts (cells × genes), or raw counts in .layers["counts"]
    .obs: per-cell metadata (will be augmented by QC, doublet scores, clusters, etc.)
    .var: per-gene metadata, including a gene symbol column (e.g. gene_symbols or feature_name)

### 3.2. Important note on dataset naming
Despite the filename, the analysis shows that the sample behaves more like a peripheral blood / PBMC immune dataset than canonical bone marrow. This conclusion is based on:

- Observed cell-type composition
- Lack of strong stem/progenitor and erythroid compartments
- Agreement with known PBMC scRNA-seq profiles

This is documented in the Findings section below.

---

## 4. Analysis Workflow
This section summarises the complete scRNA-seq analysis pipeline implemented in the notebook.

### 4.1. Data loading and identifier standardisation
Load the count matrix as an AnnData object from bone_marrow.h5ad using scanpy.read_h5ad.

Ensure unique gene and cell identifiers:
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

Why it matters:
Non-unique cell or gene identifiers break indexing and downstream operations. Standardising names is required for a robust, reproducible pipeline.


### 4.2. Gene-level QC flags and per-cell QC metrics
Select a gene symbol column from .var (e.g. gene_symbols or feature_name) and convert to upper case.

Define gene categories:
- Mitochondrial genes: symbols starting with MT-
- Ribosomal genes: symbols starting with RPS or RPL
- Hemoglobin genes: predefined list, e.g. HBA1, HBA2, HBB, HBD, HBG1, HBG2

Compute per-cell QC metrics using:

    sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['MT', 'RIBO', 'HB'],
    inplace=True
        )

This produces, among others:

    n_genes_by_counts
    total_counts
    pct_counts_MT
    pct_counts_RIBO
    pct_counts_HB (if haemoglobin genes are present)

Why it matters:
These metrics identify low-quality cells (few genes, low counts, high mitochondrial percentage) and potential technical biases that should be removed before clustering.

### 4.3. Doublet detection with Scrublet

    Extract raw counts from .X or .layers["counts"].

Run Scrublet:

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

Store results in per-cell metadata:

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

Why it matters:
Doublets (two cells in one droplet) produce hybrid expression profiles and can form artificial “mixed” clusters, especially at higher resolution. Flagging or removing them improves cluster purity and interpretability.


### 4.4. Cell-level QC filtering
Define thresholds based on QC plots (violins, scatterplots):

Examples (to be tuned to the dataset):

- n_genes_by_counts: minimum and maximum range (e.g. 200–5000)
- total_counts: minimum and maximum range (e.g. 500–50,000)
- pct_counts_MT: maximum allowed mitochondrial percentage (e.g. ≤ 20%)

Construct a Boolean mask:

    cell_mask = (
      (adata.obs['n_genes_by_counts'] >= MIN_GENES) &
      (adata.obs['n_genes_by_counts'] <= MAX_GENES) &
     (adata.obs['total_counts'] >= MIN_COUNTS) &
     (adata.obs['total_counts'] <= MAX_COUNTS) &
      (adata.obs['pct_counts_MT'] <= MAX_MT)
    )

Subset the object:

    adata = adata[cell_mask].copy()

Optionally remove doublets:

    adata = adata[adata.obs['predicted_doublet'] == False].copy()

Why it matters:
Filtering ensures the downstream analysis is performed on high-quality singlet cells, reducing noise and artifactual structure.

### 4.5. Normalisation, log-transform, and highly variable genes (HVGs)

Save raw counts for reference:

    adata.layers['counts'] = adata.X.copy()

Normalise and log-transform:

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

Identify highly variable genes:

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
    adata = adata[:, adata.var['highly_variable']].copy()

Why it matters:
Normalisation and log-transformation stabilise the variance structure, while restricting to HVGs focuses the analysis on informative genes and improves downstream dimensionality reduction and clustering.

### 4.6. PCA, neighbourhood graph, UMAP, and Leiden clustering

Run PCA on the HVG expression matrix:

    sc.tl.pca(adata, svd_solver='arpack')

Construct a k-nearest neighbour graph:

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)

Compute UMAP embedding:

    sc.tl.umap(adata, random_state=0)

Run Leiden clustering at one or more resolutions (e.g. 0.02, 0.5, 2.0):

    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_res0_5')

Why it matters:
PCA and neighbours capture the main structure of variation. UMAP provides a 2D visualisation of this manifold, and Leiden clustering identifies transcriptionally similar cell groups (candidate cell types or states).

### 4.7. Cell-type scoring with decoupler and PanglaoDB

Retrieve PanglaoDB marker signatures via decoupler.

Convert the marker table into a network format with columns:
- source: cell type
- target: gene symbol

Run the ULM method:

    dc.mt.ulm(data=adata, net=markers, tmin=3)

This stores per-cell activity scores for each cell-type signature in adata.obsm ['score_ulm'].

Aggregate scores at the cluster level using dotplots or matrixplots to see which signatures dominate each Leiden cluster.

Why it matters:
Cell-type scoring provides a reference-based, quantitative measure of how strongly each cell or cluster matches known cell-type programs (e.g. neutrophils, monocytes, B cells, T/NK cells), instead of relying on a few manually chosen markers.

### 4.8. Differential expression and marker-based validation
Use sc.tl.rank_genes_groups to find marker genes per cluster.

Compare top-ranked genes with canonical lineage markers:
- Neutrophils (e.g. S100A8, S100A9, LST1)
- Monocytes/macrophages (myeloid markers)
- B cells (e.g. MS4A1, CD79A)
- T/NK cells (e.g. CD3D, TRAC, NKG7, GNLY)

Why it matters:
Differential expression validates cell-type assignments and refines cluster identities beyond what signature scoring alone can provide.

### 4.9. Cluster annotation
Choose a main clustering resolution (e.g. leiden_res0_5).

Combine:
- decoupler-derived cell-type scores
- cluster-level marker genes
- known immunology of blood/bone marrow

Define a mapping:

    cluster_to_celltype = {
       "0": "Neutrophils",
      "1": "Monocytes / Macrophages",
      "2": "B cells",
      "3": "T / NK cells",
      "4": "Dendritic cells",
       "5": "Other/rare population",
    }
Apply it:

    adata.obs['cell_type'] = (
    adata.obs['leiden_res0_5'].map(cluster_to_celltype)
    ).astype('category')

Why it matters:
Cluster annotation translates numeric cluster IDs into interpretable biological populations, enabling comparison to external datasets and literature.

### 4.10. Interpretation and dataset-type inference

Summarise the final composition of annotated cell types.

Compare the observed mixture to reference profiles of:
- Healthy bone marrow (HSPCs, erythroid progenitors, megakaryocytes, etc.)
- PBMC / peripheral blood (mature T, B, NK, monocytes, with minimal HSPCs/erythroid cells)
- Infer whether the dataset behaves like bone marrow, PBMC, or another immune compartment, and justify this with the observed distribution.

Why it matters:
Validating the dataset’s tissue-of-origin against known biology protects against mislabeled data and strengthens the interpretation of results.

---

## 5. Findings and Interpretation

### 5.1. What cell types did you identify

From the final cell_type annotations and the proportion table, the dataset contains the following main cell types (collapsed across clusters):

- Nuocytes (type 2 innate lymphoid-like cells)
- NK cells
- Gamma delta T cells (γδ T cells)
- Neutrophils
- Naive B cells
- Plasma cells
- Monocytes
- Platelets / platelet-like cells
- Thymocytes (immature T-cell-like population)

Those are the dominant identities the decoupler/PanglaoDB + clustering pipeline resolved.


### 5.2. Biological role of each cell type (concise)

- Nuocytes (ILC2-like cells): Innate lymphoid cells that rapidly produce type 2 cytokines and help shape downstream Th2 responses.

- NK cells: Innate cytotoxic lymphocytes that kill virally infected or stressed cells without prior antigen priming.

- Gamma delta T cells (γδ T cells): Unconventional T cells with γδ T-cell receptors that recognise non-classical ligands, bridge innate and adaptive immunity, respond rapidly to stress/infection, and can produce inflammatory cytokines or exert cytotoxicity.

- Neutrophils: Short-lived granulocytes generated in bone marrow and released into the blood. First-line phagocytes that rapidly migrate to infected or damaged sites, engulf microbes.

- Naive B cells: Antigen-inexperienced B lymphocytes that circulate through the blood and secondary lymphoid organs. Upon encountering antigen and T-cell help, they undergo activation, class switching and differentiation into memory B cells or plasma cells.

- Plasma cells: Terminally differentiated B cells specialised for high-rate antibody secretion. They maintain humoral immunity, reside in bone marrow and secondary lymphoid tissues.

- Monocytes: Circulating myeloid cells that patrol the blood and can migrate into tissues to differentiate into macrophages or dendritic cells. They phagocytose pathogens and debris and secrete inflammatory mediators.

- Platelets / platelet-like cells : Anucleate fragments derived from megakaryocytes. They mediate hemostasis and thrombosis, but also interact with leukocytes and endothelium and modulate inflammatory responses.

- Thymocytes (immature T-cell–like cells): T-lineage cells in early or intermediate stages of differentiation (normally located in the thymus).

### 5.3. Is the tissue source really bone marrow?

Conclusion:
On biological grounds, the dataset is much more consistent with a peripheral blood / PBMC-like immune sample than with true bone marrow.

#### 5.3.1. What we actually see

- From the proportions you reported:
  - Nuocytes: 32.6 %
  -  NK cells: 23.9 %
   - γδ T cells: 14.3 %
   - Neutrophils: 8.3 %
   - Naive B cells: 7.3 %
   - Plasma cells: 5.5 %
    - Monocytes: 5.4 %
  -  Platelets: 1.8 %
   - Thymocytes: 0.8 %

- Key features:
  - Dominance of mature lymphoid and innate lymphoid populations.
  - Presence of mature myeloid cells (neutrophils, monocytes, DC-like signatures).
   - Modest plasma cell and platelet compartments.
  -  No clear, large CD34+ HSPC cluster and no obvious multi-stage erythroid progenitor/erythroblast series with strong haemoglobin expression.

#### 5.3.2. What true bone marrow should look like

- Healthy bone marrow single-cell profiles typically include:
  - A substantial hematopoietic stem and progenitor cell (HSPC) compartment (CD34+, early progenitors).
  - Multiple erythroid progenitor and maturation clusters (proerythroblasts → erythroblasts → reticulocytes) with high globin expression.
  - Megakaryocyte progenitors, in addition to platelets.
  - A mixture of developing and mature lymphoid and myeloid cells.

In other words, bone marrow is characterised not just by mature immune cells but by visible differentiation hierarchies, especially HSPCs and erythroid lineages.

#### 5.3.3. Comparison

The dataset is dominated by mature immune effectors (innate lymphoid, NK, γδ T, neutrophils, monocytes, B, plasma, platelets) and lacks a robust HSPC or erythroid progenitor architecture.

This pattern matches blood / PBMC-type samples, where one expects mainly mature T, B, NK, monocytes, neutrophils, small plasma-cell fractions and platelets, with very few progenitors.

Therefore, despite the file name, the most biologically coherent interpretation is:
  - The sample represents a blood-derived immune cell compartment (PBMC-like / peripheral leukocytes) rather than canonical bone marrow.


### 5.4. Based on the relative abundance of cell types, is the patient healthy or infected?

- Given the proportions:
  -  Neutrophils: 8.3 %
  - Monocytes: 5.4 %
  - Combined myeloid (neutrophils + monocytes): ~13.7 %
  -  NK cells: 23.9 %
  - Nuocytes + NK + γδ T (innate/innate-like lymphoid): ~71 %
  - Naive B + plasma cells: ~12.8 %
  -  Thymocyte-like: 0.8 %


#### 5.4.1. What an acute, infected / highly inflamed profile would look like

In acute systemic infection or strong inflammatory states, typically, we would see, in the blood:
- Neutrophilia and/or inflammatory monocyte expansion, neutrophils and monocytes becoming a very large fraction of total leukocytes.
- Often relative lymphopenia (drop in total lymphocyte fraction).
- NK and CD8 T cells may be highly activated and sometimes expanded, but usually in the context of a strong myeloid response and/or global IFN signature.


#### 5.4.2. What the proportions suggest

- Neutrophils and monocytes are not expanded: together they are ~14 %, not the dominant population. There is no myeloid takeover that would flag severe acute bacterial sepsis–type physiology.
- Lymphoid / innate lymphoid populations dominate (nuocytes, NK, γδ T), and B cells + plasma cells are well represented. There is no lymphopenia; if anything, lymphoid compartments are abundant.
- A ~5–6 % plasma cell fraction indicates that the immune system has seen an antigen and is producing antibodies, but this is compatible with past infection, vaccination or low-grade/chronic exposure and does not by itself prove acute infection.

There isn't an overwhelming interferon-stimulated gene signature or extremely activated NK/CD8 states across most cells, which would strengthen an “acutely infected” interpretation.


#### 5.4.3. Defensible call

Based on proportions alone:
- The absence of neutrophil/monocyte expansion and the preservation (even dominance) of lymphoid populations argue against a classic acute systemic infection profile.
- The dataset looks more like a peripheral immune compartment in a non-acutely infected or mildly activated state: robust lymphoid presence, modest myeloid fractions, some plasma cells and platelets.


### Conclusion
Using only relative cell-type abundances, the dataset is more consistent with a non-acute/near-baseline immune state than with a frank acute infection. There may be evidence of immunological experience or mild activation (plasma cells, strong innate-like lymphoid representation), but there is no clear neutrophil/monocyte expansion or lymphocyte collapse that would justify calling this an overtly infected patient from these data alone.


### References
    *(Karolyn A. Oetjen, J. Philip McCoy, Christopher S. Hourigan. Human bone marrow assessment by single-cell RNA sequencing, mass cytometry, and flow cytometry.  Published December 6, 2018. JCI Insight. 2018;3(23):e124928. https://doi.org/10.1172/jci.insight.124928.)*
     *( Xiaowei Xie, Mengyao Liu, Yawen Zhang, Bingrui Wang, Caiying Zhu, Chenchen Wang, Qing Li, Yingying Huo, Jiaojiao Guo, Changlu Xu, Linping Hu, Aiming Pang, Shihui Ma, Lina Wang, Wenbin Cao, Shulian Chen, Qiuling Li, Sudong Zhang, Xueying Zhao, Wen Zhou, Hongbo Luo, Guoguang Zheng, Erlie Jiang, Sizhou Feng, Lixiang Chen, Lihong Shi, Hui Cheng, Sha Hao, Ping Zhu, Tao Cheng, Single-cell transcriptomic landscape of human blood cells, National Science Review, Volume 8, Issue 3, March 2021, nwaa180, https://doi.org/10.1093/nsr/nwaa180)*
    *(Pellin, D., Loperfido, M., Baricordi, C. et al. A comprehensive single cell transcriptional landscape of human hematopoietic progenitors. Nat Commun 10, 2395 (2019). https://doi.org/10.1038/s41467-019-10291-0)*
    *(Oelen, R., de Vries, D.H., Brugge, H. et al. Single-cell RNA-sequencing of peripheral blood mononuclear cells reveals widespread, context-specific gene expression regulation upon pathogenic exposure. Nat Commun 13, 3267 (2022). https://doi.org/10.1038/s41467-022-30893-5)*


# Summary
A full scRNA-seq analysis pipeline was implemented, including QC, doublet detection, normalisation, dimensionality reduction, clustering, and cell-type annotation.
The pipeline identified major immune cell populations and produced a coherent map of the sample’s immune composition.
Based on the cell-type landscape, the dataset is most consistent with a PBMC-like / peripheral blood leukocyte sample, not canonical bone marrow, illustrating the importance of validating metadata labels against biological evidence.
"""
