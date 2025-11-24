from cluster_map.preprocess_spatial_data import preprocess_registration_output
from cluster_map.preprocess_atlas import (
    preprocess_scrnaseq,
    get_cluster_profiles,
    visualize_profiles,
)
from cluster_map.bayesian_prob import assign_clusters

from pathlib import Path

# from platy_reg.n5_utils import read_volume, write_volume, get_attrs
from skimage.transform import rescale
import numpy as np
import pandas as pd
import argparse
import logging
import sys
import matplotlib.pyplot as plt


def calculate_spatial_specificity_weights(
    smfish_df: pd.DataFrame, min_weight: float = 0.25, max_weight: float = 1.75
) -> np.ndarray:
    """
    Calculates a weight for each gene based on how many nuclei express it.
    """
    logging.info("Calculating spatial gene specificity weights...")

    num_genes = smfish_df.shape[1]
    num_nuclei = smfish_df.shape[0]
    nuclei_counts = smfish_df.sum(axis=0).to_numpy()
    ubiquity_score = nuclei_counts / num_nuclei
    specificity_score = 1 - ubiquity_score
    spatial_weights = min_weight + specificity_score * (max_weight - min_weight)

    logging.info(
        f"Spatial weights calculated, ranging from {np.min(spatial_weights):.2f} to {np.max(spatial_weights):.2f}"
    )

    return spatial_weights


def calculate_gene_specificity_weights(
    profiles: dict,
    gene_list: list,
    expression_threshold: float = 0.1,
    min_weight: float = 0.25,
    max_weight: float = 1.75,
) -> np.ndarray:
    """
    Calculates a weight for each gene based on its expression specificity across clusters.
    ... (rest of the docstring) ...
    """
    logging.info("Calculating gene specificity weights...")

    num_genes = len(gene_list)
    num_clusters = len(profiles)
    profile_matrix = np.array(
        list(profiles.values())
    )  # Shape: (num_clusters, num_genes)
    is_expressed_matrix = profile_matrix > expression_threshold
    gene_expression_counts = is_expressed_matrix.sum(axis=0)

    # Calculate "ubiquity score" (0 = perfectly specific, 1 = in all clusters)
    # We use max(1, ...) to avoid division by zero if there's only one cluster
    ubiquity_score = (gene_expression_counts - 1) / max(1, num_clusters - 1)
    ubiquity_score = np.clip(ubiquity_score, 0, 1)  # Ensure score is between 0 and 1

    # Invert to get "specificity score" (1 = perfectly specific, 0 = ubiquitous)
    specificity_score = 1 - ubiquity_score

    # Map specificity score to the desired weight range
    gene_weights = min_weight + specificity_score * (max_weight - min_weight)

    logging.info(
        f"Gene weights calculated, ranging from {np.min(gene_weights):.2f} to {np.max(gene_weights):.2f}"
    )

    return gene_weights, gene_expression_counts


def main():
    parser = argparse.ArgumentParser(
        description="""Assign scRNAseq atlas cluster id to nuclei in EM."""
    )
    # Replaced single positional arg with two optional args
    parser.add_argument(
        "--scrnaseq_atlas", type=str, required=True, help="Path of the scRNAseq atlas"
    )
    parser.add_argument(
        "--gene_metadata",
        type=str,
        required=True,
        help="Correspondence between gene names and XLOCs",
    )
    parser.add_argument(
        "--cluster_column",
        type=str,
        required=True,
        help="Column in obs to use as cluster ids",
    )

    # Optional file arguments for both techniques
    parser.add_argument(
        "--smfish_data", type=str, default=None, help="Path to smFISH expression table."
    )
    parser.add_argument(
        "--prospr_data",
        type=str,
        default=None,
        help="Path to ProSPeR expression table.",
    )

    # Log directory and path remain mandatory
    parser.add_argument(
        "log_dir", type=str, help="Directory to store diagnostic plots etc"
    )
    parser.add_argument("log_path", type=str, help="Path to log file")

    # Separate threshold arguments
    parser.add_argument(
        "--smfish_thresh",
        type=float,
        default=20,
        help="Binarization threshold for smFISH data.",
    )
    parser.add_argument(
        "--prospr_thresh",
        type=float,
        default=0.1,
        help="Binarization threshold for ProSPR data.",
    )

    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Add weighting to probability assignment.",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="If passed, visualizes all plots, including cluster profiles",
    )

    parser.add_argument(
        "--run_diagnostics",
        action="store_true",
        help="If passed, generate diagnostic files in a 'diagnostics' subdirectory.",
    )
    parser.add_argument(
        "--profile_constant",
        type=float,
        default=0.001,
        help="Base value for cluster profiles to avoid probabilities of 0 or 1. Default: 0.001",
    )

    args = parser.parse_args()

    if not args.smfish_data and not args.prospr_data:
        logging.error(
            "At least one input file (--smfish_data or --prospr_data) must be provided."
        )
        sys.exit(1)

    log_dir = Path(args.log_dir)
    log_dir.mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(args.log_path, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    binarized_dfs = []

    diagnostics_path = log_dir / "diagnostics"

    # Process smFISH data (if provided)
    if args.smfish_data:
        logging.info(
            f"Preprocessing confocal smFISH data with threshold {args.smfish_thresh}"
        )
        df_smfish = preprocess_registration_output(
            args.smfish_data,
            log_dir / "qc_smfish_confocal",
            visualize=args.visualize,
            count_threshold=args.smfish_thresh,
            technique="smfish",
        )
        # SAVE THE PROCESSED SMFISH DF
        if args.run_diagnostics:
            diagnostics_path.mkdir(exist_ok=True)  # Ensure the folder exists
            smfish_diagnostic_path = (
                diagnostics_path / "smfish_confocal_processed_diagnostic.tsv"
            )
            df_smfish.to_csv(
                smfish_diagnostic_path, sep="\t", index=True, index_label="label_id"
            )
            logging.info(
                f"Saved processed smFISH diagnostic to: {smfish_diagnostic_path}"
            )
        binarized_dfs.append(df_smfish)

    # Process prospr data (if provided)
    if args.prospr_data:
        logging.info(f"Preprocessing ProSPR data with threshold {args.prospr_thresh}")
        df_prospr = preprocess_registration_output(
            args.prospr_data,
            log_dir / "qc_smfish_prospr",
            visualize=args.visualize,
            count_threshold=args.prospr_thresh,
            technique="prospr",
        )
        binarized_dfs.append(df_prospr)

    # Combining the binarized tables

    if len(binarized_dfs) == 1:
        smfish_df_combined = binarized_dfs[0]
    else:
        # Concatenate and sum: sum > 0 means the gene is expressed in at least one method
        smfish_df_combined = pd.concat(binarized_dfs, axis=0).groupby(level=0).max()
        logging.info("Combined smFISH and ProSPeR data using OR logic.")

    if args.run_diagnostics:
        diagnostics_path.mkdir(exist_ok=True)
        output_combined_path = diagnostics_path / "smfish_prospr_combined_binarized.tsv"
        smfish_df_combined.to_csv(
            output_combined_path, sep="\t", index=True, index_label="label_id"
        )
        logging.info(f"Saved combined binarized diagnostic to: {output_combined_path}")

    # preprocess scrnaseq

    logging.info(f"Quality control and preprocessing of the scRNAseq data...")
    scrnaseq_ad = preprocess_scrnaseq(
        args.scrnaseq_atlas,
        args.gene_metadata,
        args.cluster_column,
        smfish_df_combined.columns,  # Use combined columns for filtering
        log_dir / "qc_scrnaseq",
        visualize=args.visualize,
    )

    # Take only columns of smfish df which were found in scrnaseq atlas metadata
    gene_list = scrnaseq_ad.var.pb_name
    # gene_list = scrnaseq_ad.var.pb_name.tolist() # Get the final scRNA-seq gene list

    if args.run_diagnostics:
        diagnostics_path.mkdir(exist_ok=True)
        # Save the definitive scRNA-seq gene list
        pd.Series(gene_list).to_csv(
            diagnostics_path / "final_scrnaseq_gene_list.tsv",
            sep="\t",
            index=False,
            header=False,
        )
        logging.info(
            f"Final scRNA-seq gene list saved to: {diagnostics_path / 'final_scrnaseq_gene_list.tsv'}"
        )

        # Save the smFISH combined columns before final filtering
        pd.Series(smfish_df_combined.columns).to_csv(
            diagnostics_path / "smfish_combined_columns_unfiltered.tsv",
            sep="\t",
            index=False,
            header=False,
        )
        logging.info(
            f"Unfiltered combined smFISH columns saved to: {diagnostics_path / 'smfish_combined_columns_unfiltered.tsv'}"
        )

    smfish_df_combined = smfish_df_combined[gene_list]
    logging.info(f"Final gene list size: {len(gene_list)}")

    # Save the FINAL, filtered smFISH data for diagnostic use
    if args.run_diagnostics:
        diagnostics_path.mkdir(exist_ok=True)
        final_smfish_path = diagnostics_path / "smfish_prospr_combined_filtered.tsv"
        smfish_df_combined.to_csv(
            final_smfish_path, sep="\t", index=True, index_label="label_id"
        )
        logging.info(
            f"Saved final filtered expression data for diagnostics to: {final_smfish_path}"
        )

    # Prepare the expression profiles
    profiles = get_cluster_profiles(
        scrnaseq_ad,
        args.cluster_column,
        log_dir / "profiles",
        base_value=args.profile_constant,
    )
    total_profile = np.array(
        scrnaseq_ad.X.sum(axis=0) / scrnaseq_ad.X.shape[0]
    ).flatten()
    if args.visualize:
        visualize_profiles(
            gene_list, {**profiles, "total": total_profile}, log_dir / "profiles"
        )
    else:
        logging.info(
            "Skipping cluster profile plot generation. Use the --visualize flag to enable."
        )

    if args.weighted:
        logging.info(
            "Applying TWO-WAY gene specificity weighting as requested (--weighted)."
        )

        sc_weights, sc_cluster_counts = calculate_gene_specificity_weights(
            profiles, gene_list
        )
        spatial_weights = calculate_spatial_specificity_weights(smfish_df_combined)

        final_gene_weights = sc_weights * spatial_weights

        if args.run_diagnostics:
            weights_df = pd.DataFrame(
                {
                    "gene_name": gene_list,
                    "sc_cluster_count": sc_cluster_counts,
                    "sc_weight": sc_weights,
                    "spatial_weight": spatial_weights,
                    "final_weight": final_gene_weights,
                }
            )
            weights_df = weights_df.sort_values(by="final_weight", ascending=False)

            diagnostics_path.mkdir(exist_ok=True)
            weights_output_path = diagnostics_path / "gene_specificity_weights.tsv"
            weights_df.to_csv(weights_output_path, sep="\t", index=False)
            logging.info(
                f"Saved detailed two-way gene weights to: {weights_output_path}"
            )

    else:
        logging.info(
            "Using standard (unweighted) probabilistic assignment. To use weights, add the --weighted flag."
        )
        final_gene_weights = np.ones(len(gene_list))

    # Given the binarized smFISH expression and expression profiles, assign cluster probability to each nucleus in EM
    assign_clusters(
        smfish_df_combined,
        gene_list,
        profiles,
        total_profile,
        log_dir / "results",
        args.cluster_column,
        gene_weights=final_gene_weights,
        run_diagnostics=args.run_diagnostics,
        diagnostics_dir=diagnostics_path,
    )


if __name__ == "__main__":
    main()
