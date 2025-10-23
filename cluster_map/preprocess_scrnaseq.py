import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import logging
from pathlib import Path
import numpy as np
import pickle

def preprocess_scrnaseq(scrnaseq_atlas_path, gene_metadata_path, cluster_column, gene_list, out_dir, count_threshold=2, visualize=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    
    adata = sc.read_h5ad(scrnaseq_atlas_path)
    logging.info(adata)

    # Normalize
    adata.layers['counts'] = adata.layers.pop('RNA_counts')
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if visualize:
        sc.pl.umap(
        adata,
        color=cluster_column,
        size=2,
        show=False
        )
        plt.savefig(out_dir / "umap.png", dpi=300, bbox_inches="tight")
        plt.close("all")

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

    thresholds = []
    for name in subset_adata.var.index:
        expr = subset_adata[:, name].X.toarray().T[0]
        expr_sorted = np.sort(expr)
        # Calculate the median of non-zero expression values
        thr = np.median(expr_sorted[expr_sorted > 0])
        thresholds.append(thr)
        
    # Now, if visualize is True, create the plots using the calculated thresholds.
    if visualize:
        logging.info("Generating cumulative expression plots...")
        for i, name in enumerate(subset_adata.var.index):
            print(f"Plot cumulative expression {name}")
            plt.figure(figsize=(6, 5))
            expr = subset_adata[:, name].X.toarray().T[0]
            expr_sorted = np.sort(expr)
            
            # Use the pre-calculated threshold for the plot
            thr = thresholds[i] 
            
            plt.plot(expr_sorted, label=name)
            plt.hlines(y=thr, xmin=0, xmax=len(expr_sorted), color="red")
            plt.legend()
            plt.title("Cumulative distribution of expression levels")
            plt.ylabel("Normalized expression level")
            plt.xlabel("Cell number")
            plt.savefig(out_dir / f"{subset_adata.var.loc[name, 'pb_name']}_norm_count_cumulative.png", dpi=300, bbox_inches="tight")
            plt.close("all")
                 
    

    logging.info(subset_adata)
    logging.info(subset_adata.var)

    subset_adata.write(out_dir / "scrnaseq_norm_expr.h5ad")

    subset_adata_bin = subset_adata.copy()
    subset_adata_bin.X = (subset_adata_bin.X > np.array(thresholds)).astype(np.int16)
    subset_adata_bin.write(out_dir / "scrnaseq_bin_expr.h5ad")
    
    return subset_adata_bin


def get_cluster_profiles(adata, cluster_column, out_dir, base_value: float = 0.001):

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
        nonzero_counts = adata_aggregated[cluster_name, :].layers["count_nonzero"]
        
        # This ensures P is never exactly 0 or 1.
        profiles[cluster_name] = (nonzero_counts + base_value) / (cluster_sizes[cluster_name] + 2 * base_value)
        profiles[cluster_name] = np.array(profiles[cluster_name]).flatten()
        

    # Create empty profile for the nuclei where nothing was measured
    profiles["zero"] = np.zeros(adata.shape[1]) + 0.0000001

    # Create profile that expresses everything for the autofluorescent areas
    profiles["autofluorescence"] = np.ones(adata.shape[1]) * 0.99

    # INSERTION POINT: SAVE THE PROFILES DICTIONARY
    output_file = out_dir / "profiles_dictionary.pkl"
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(profiles, f)
        logging.info(f"Successfully saved cluster profiles dictionary to: {output_file}")
    except Exception as e:
        logging.error(f"Failed to save cluster profiles via pickle: {e}")

    return profiles


""" def visualize_profiles(gene_list, profiles, out_dir):
    n_clusters = len(profiles)
    fig, axs = plt.subplots(n_clusters, figsize=(12, n_clusters * 2 + 3))
    # fig.suptitle('Expression profiles of clusters')
    fig.tight_layout()
    
    for i, k in enumerate(profiles):
        axs[i].set_title(k)
        axs[i].bar(gene_list, profiles[k])
        axs[i].set_ylim(0, 1)
        axs[i].tick_params("x", rotation=90, labelsize=4)

    
    plt.tight_layout(pad=1.5, h_pad=2)


    plt.tight_layout(pad=1.5, h_pad=2)

    plt.savefig(out_dir / "expression_profiles.png", dpi=300)
    plt.close("all") """

def visualize_profiles(gene_list, profiles, out_dir):
    """
    Generates and saves the gene expression profile for each cluster as a separate PNG file.

    Args:
        gene_list (list): List of gene names (for the x-axis).
        profiles (dict): Dictionary where keys are cluster names and values are 
                         the expression profile arrays (0-1).
        out_dir (Path): The directory path (e.g., log_dir / "profiles") to save the plots.
    """
    logging.info(f"Generating and saving individual expression profile plots for {len(profiles)} clusters...")
    
    # Ensure the output directory exists
    out_dir.mkdir(parents=True, exist_ok=True)

    # Iterate through each cluster profile
    for cluster_name, profile_data in profiles.items():
        
        # 1. Create a new figure for each cluster
        # Use a reasonable size that accommodates all genes on the x-axis
        # Size suggestion: 0.15 inches per gene + margin, but we'll use a fixed wide size.
        fig, ax = plt.subplots(figsize=(15, 6)) 

        # 2. Plot the profile data
        ax.bar(gene_list, profile_data, color='darkblue', alpha=0.9)
        
        # 3. Set limits and labels
        ax.set_ylim(0, 1)
        ax.set_title(f'Cluster Profile: {cluster_name}', fontsize=14)
        ax.set_ylabel('Proportion of Cells Expressing Gene')
        ax.set_xlabel('Gene')
        
        # Rotate x-axis ticks for readability (important for many genes)
        ax.tick_params("x", rotation=90, labelsize=6)
        
        # 4. Save the figure
        # Sanitize the filename to handle spaces/slashes in cluster names
        safe_name = cluster_name.replace(' ', '_').replace('/', '_').replace('\\', '_')
        plot_path = out_dir / f"profile_{safe_name}.png"
        
        # Use tight_layout only for this specific figure
        plt.tight_layout() 
        plt.savefig(plot_path, dpi=300)
        plt.close(fig) # Close the figure to free memory

    logging.info(f"Successfully saved {len(profiles)} individual cluster profile plots to {out_dir}")


    # for idx, name in enumerate(gene_list):
    #     plt.figure()
    #     plt.bar(cluster_ids, nonzero_cluster[idx, :], color="yellowgreen")
    #     plt.xlabel("Cluster ID")
    #     plt.ylabel("Proportion of cells with nonzero expression")
    #     plt.title(name, fontsize=20)
