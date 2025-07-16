# Predicting cell type identity from Platybrowser

We are attempting to assign cell type identities in an electron microscopy atlas, Platybrowser. For that, we have registered datasets of RNA stainings, and a scRNAseq atlas with known cell types. Using shared key markers, we will assign probabilities to each nuclei of being part of a certain tissue or a certain granular cell type. This is a "spatial transcriptomic"-like approach.

## Motivation

### Source data and processing

We are using a sub-atlas of 30850 cells, which have been curated to remove developmental precursors and have clearly defined identities. They were picked by Phil Oel/Detlev Arendt. We have ~25,000 protein coding genes, but ~60,000 XLOC features in our transcriptome. This is single nuclei data, which may not always be in full complete agreement with RNA staining methods used for Platybrowser, which don't differentiate cytoplasmic vs nuclear transcripts.

Originally, we have a Seurat object with various extra assays (`SCT`, `integrated`, used for batch correction) and dimensional reductions. I keep only the raw counts (`RNA` assay), and the same PCA and UMAP embeddings/loadings in the paper to ease the conversion to Anndata. I then load it in Python, rename `RNA` to the more canonical `counts` and reperform minor steps that were lost in translation (`log1p` and `normalization`). There is no new dimensional reductions / processing, or search for variable features because this is already done in R for the main paper.

The key metadata are:
* `group_full` (plus `_slug`/`_color` for plotting color keys) related to large tissue level divisions of cell types. There are ~20 groups like this, but some like the neuron types are hard to separate. Some special groups like enteric neurons have transcriptomic similarities to both neurons and the midgut.
* `clade`, `subclade`, `subsubclade` related to how the object was iteratively clustered in clade and subclade
* `curated_high_res_annotations` which is the cell type level information, like `cladeXsubYsubsubZ`. There are ~250 cell types, some quite tiny.

### Analysis

* `s1-alyona-python-side.ipynb` Syncing Platybrowser gene symbols like Arx and `XLOC` gene identifiers across single cell and Platybrowser
* `s2-alyona-R-side.ipynb` Seurat to Anndata conversion
* `s3-alyona-scanpy.ipynb` Rough Anndata processing and comparisons with Platybrowser

## Contents

- `./env_conda/`: 
    * `./env_conda/6dpf_paper.yaml`: R code to read in Seurat. The conversion to `anndata` requires https://github.com/settylab/convert2anndata whiich has to be manually installed, but is very fast.
    * `./env_conda/env-python.yaml`: Scanpy etc... for the Python side of scRNAseq

    Rscript -e "convert2anndata::cli_convert()"  --input /scratch/cros/curated_6dpf_cleaner.RDS --output /scratch/cros/seurat_curated_6dpf.h5ad


- `./prepped_data/`: data where the XLOC are Platybrowser names
    * `./prepped_data/genes_metadata.tsv`  Gene metadata for reference, 1st thing generated
    * `./prepped_data/platybrowser_counts.tsv`  Platybrowser expression matrix, renamed. Values are in $[0,1]$.

- `./source_data/`: no proper XLOC or genes we don't have spatial transcriptomics data for
    * `./source_data/platybrowser/*` Relevant source contents from /g/arendt/EM_6dpf_segmentation/platy-browser-data/data/1.0.1/tables/sbem-6dpf-1-whole-segmented-cells/
    * `./source_data/6dpf neuron types - Prospr_db_copy.tsv`  Old Prospr/Platybrowser information
    * `./source_data/markers_wilcox_202407_CuratedTotal6dpfAtlas_31knuclei_all_pos_pct2_0.04.csv` Yuyao from our lab use the Seurat FindMarkers function to identify good marker genes. The format is for each cell type a line with a marker gene. These are "positive" markers (logFC > 0, upregulated). You want high pct1 (% cells expressing the marker in the cluster) and low pct2 (% cells expressing outside the cluster). This provides a way of binarizing expression by extracting `cell_type, marker gene` and turning it into a matrix.
    * `./source_data/seurat_curated_6dpf.h5ad` Full Seurat object with RNA assay and PCA/UMAP, very little processing done. **In gitignore**
