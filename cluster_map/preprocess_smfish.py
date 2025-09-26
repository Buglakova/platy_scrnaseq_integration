import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import anndata as ad
import h5py
import logging
from pathlib import Path

from cluster_map.plot_utils import plot_cumulative_expression


def percentile_norm(expr, pmin, pmax, eps=1e-10):
    """Percentile normalization of an expression vector. Values are clipped to [0, 1] range
    Args:
        expr: expression vector
        pmin: in percents
        pmax: in percents
        eps: _description_. Defaults to 1e-10.

    Returns:
        _description_
    """

    pmin = np.percentile(expr, pmin)
    pmax = np.percentile(expr, pmax)
    return np.clip(((expr.astype(np.float32) - pmin) / (pmax - pmin + eps)), a_min=0, a_max=1)


def preprocess_registration_output(aggregated_table_path, out_dir, count_threshold=20, visualize=True):

    expression_matrix_path = aggregated_table_path
    df = pd.read_csv(expression_matrix_path, sep="\t")
    df = df.set_index("label_id")
    columns = [col for col in df.columns if col.endswith("spots")]
    df_expr = df[columns]
    df_expr = df_expr.rename(columns={col: col.split("_")[0].lower() for col in df_expr.columns})

    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    logging.info(f"Quality control and preprocessing of the smFISH data")

    pd.set_option('display.max_columns', None)
    logging.info("\n" + repr(df_expr.describe().T))

    if visualize:
        plt.figure(figsize=(10, 1 * len(df_expr.columns)))
        sns.stripplot(df_expr, jitter=0.4, alpha=0.6, orient="h")
        plt.title("Distribution of raw counts")
        plt.savefig(out_dir / "raw_count_stripplot.png", dpi=300)
        plt.close("all")

        plt.figure(figsize=(10, 1 * len(df_expr.columns)))
        sns.violinplot(df_expr, alpha=0.6, orient="h", cut=0)
        plt.title("Distribution of raw counts")
        plt.savefig(out_dir / "raw_count_violinplot.png", dpi=300)
        plt.close("all")


    # Threshold by raw number of counts
    df_thr_expr = df_expr
    df_thr_expr[df_expr <= count_threshold] = 0

    # Quantile norm
    df_norm_expr = pd.DataFrame()
    df_norm_expr.index = df_thr_expr.index
    for col in df_thr_expr.columns:
        df_norm_expr[col] = percentile_norm(df_thr_expr[col], 0, 99)


    if visualize:
        plt.figure(figsize=(10, 10))
        for col in df_norm_expr.columns:
            print(col)
            plot_cumulative_expression(df_norm_expr[col], col)
            plt.tight_layout()
            plt.savefig(out_dir / "norm_count_cumulative.png", dpi=300)
            plt.close("all")

    df_bin_expr = (df_norm_expr > 0).astype(np.int16)

    logging.info("Normalized expression df stats\n" + repr(df_norm_expr.describe().T))
    logging.info("Binarized expression df stats\n" + repr(df_bin_expr.describe().T))

    logging.info("Number of cells with a given number of nonzero genes/n" + str(df_bin_expr.sum(axis=1).value_counts()))

    if visualize:
        plt.figure(figsize=(1 * len(df_expr.columns), 5))
        plt.bar(df_bin_expr.columns, df_bin_expr.sum(axis=0))
        plt.title("Nonzero cells per gene in binarized data")
        plt.xlabel("genes")
        plt.ylabel("Number of cells")
        plt.savefig(out_dir / "bin_nonzero_per_gene.png", dpi=300)
        plt.close("all")

        plt.figure(figsize=(0.5 * len(df_expr.columns), 5))
        plt.bar(df_bin_expr.sum(axis=1).value_counts().index, df_bin_expr.sum(axis=1).value_counts())
        plt.xlabel("Number of nonzero genes")
        plt.ylabel("Number of cells")
        plt.savefig(out_dir / "bin_nonzero_per_cell.png", dpi=300)
        plt.close("all")

    df_bin_expr.to_csv(out_dir / "smfish_bin_expr.tsv", sep="\t", index=True, index_label="label_id")
    df_norm_expr.to_csv(out_dir / "smfish_norm_expr.tsv", sep="\t", index=True, index_label="label_id")


    return df_bin_expr


    