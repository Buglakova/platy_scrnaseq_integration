import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np

def preprocess_scrnaseq(scrnaseq_atlas_path, gene_metadata_path, cluster_column, gene_list, out_dir, count_threshold=2, visualize=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    
    adata = sc.read_h5ad(scrnaseq_atlas_path)
    logging.info(adata)

    # Normalize
    adata.layers['counts'] = adata.layers.pop('RNA')
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)    

    if visualize:
        sc.pl.umap(
        adata,
        color=cluster_column,
        size=2,
        )
        plt.savefig(out_dir / "umap.png", dpi=300, bbox_inches="tight")

    logging.info("Number of cells per cluster/n" + str(adata.obs[cluster_column].value_counts()))

    # Assign metadata
    genes_metadata = pd.read_csv(gene_metadata_path, sep="\t")
    genes_metadata = genes_metadata[genes_metadata["Platybrowser name"].isin(gene_list)]
    logging.info(genes_metadata)
    genes_lookup = dict(zip(genes_metadata['XLOC'], genes_metadata['XLOC_full']))
    pb_names_lookup = dict(zip(genes_metadata['XLOC'], genes_metadata['Platybrowser name']))

    assert(not set(genes_lookup.keys()).difference(set(adata.var_names.to_list())))
    subset_adata = adata[:, adata.var_names.isin(genes_lookup.keys())].copy()
    subset_adata.var['XLOC'] = subset_adata.var_names
    subset_adata.var['pb_name'] = [pb_names_lookup[xloc] for xloc in subset_adata.var_names]
    subset_adata.var_names = subset_adata.var_names.to_series().map(genes_lookup).tolist()

    pd.set_option('display.max_rows', None)

    logging.info(subset_adata.var_names)

    if visualize:
        plt.figure(figsize=(6, 5))
        for name in subset_adata.var.index:
            expr = subset_adata[:, name].X.toarray().T[0]
            expr_sorted = np.sort(expr)
            plt.plot(expr_sorted, label=name)
            # plt.hist(expr, 100, cumulative=True, label=name, histtype="step", orientation="horisontal")
        plt.legend()
        plt.title("Cumulative distribution of expression levels")
        plt.ylabel("Normalized expression level")
        plt.xlabel("Cell number")
        plt.savefig(out_dir / "norm_count_cumulative.png", dpi=300, bbox_inches="tight")
                 

    logging.info(subset_adata)
    logging.info(subset_adata.var)

    subset_adata.write(out_dir / "scrnaseq_norm_expr.h5ad")

    subset_adata_bin = subset_adata.copy()
    subset_adata_bin.X = (subset_adata_bin.X > count_threshold).astype(np.int16)
    subset_adata_bin.write(out_dir / "scrnaseq_bin_expr.h5ad")
    
    return subset_adata_bin



def get_cluster_profiles(adata, cluster_column, out_dir):

    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    adata_aggregated = sc.get.aggregate(adata, by=cluster_column, func="count_nonzero")
    logging.info(adata_aggregated)
    logging.info(adata_aggregated.layers["count_nonzero"])

    cluster_sizes = adata.obs[cluster_column].value_counts()
    cluster_names = adata.obs[cluster_column].unique()
    logging.info(cluster_sizes)
    profiles = {}
    for cluster_name in cluster_names:
        profiles[cluster_name] = adata_aggregated[cluster_name, :].layers["count_nonzero"] / cluster_sizes[cluster_name]
        profiles[cluster_name] = np.array(profiles[cluster_name]).flatten()



    logging.info(profiles)

    profiles["zero"] = np.zeros_like(profiles[cluster_name])


    return profiles