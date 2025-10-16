from cluster_map.preprocess_smfish import preprocess_registration_output
from cluster_map.preprocess_scrnaseq import preprocess_scrnaseq, get_cluster_profiles, visualize_profiles
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



    
def main():
    parser = argparse.ArgumentParser(
        description="""Assign scRNAseq atlas cluster id to nuclei in EM."""
    )
    # Replaced single positional arg with two optional args
    parser.add_argument("--scrnaseq_atlas", type=str, required=True, help="Path of the scRNAseq atlas")
    parser.add_argument("--gene_metadata", type=str, required=True, help="Correspondence between gene names and XLOCs")
    parser.add_argument("--cluster_column", type=str, required=True, help="Column in obs to use as cluster ids")
    
    # Optional file arguments for both techniques
    parser.add_argument("--smfish_data", type=str, default=None, help="Path to smFISH expression table.")
    parser.add_argument("--prospr_data", type=str, default=None, help="Path to ProSPeR expression table.")
    
    # Log directory and path remain mandatory
    parser.add_argument("log_dir", type=str, help="Directory to store diagnostic plots etc")
    parser.add_argument("log_path", type=str, help="Path to log file")
    
    # Separate threshold arguments
    parser.add_argument("--smfish_thresh", type=float, default=20, help="Binarization threshold for smFISH data.")
    parser.add_argument("--prospr_thresh", type=float, default=0.1, help="Binarization threshold for ProSPR data.")

    args = parser.parse_args()
    
    if not args.smfish_data and not args.prospr_data:
        logging.error("At least one input file (--smfish_data or --prospr_data) must be provided.")
        sys.exit(1)
    
    log_dir = Path(args.log_dir)
    log_dir.mkdir(exist_ok=True)
    
    logging.basicConfig(level=logging.INFO,
             format="%(asctime)s [%(levelname)s] %(message)s",
             handlers=[logging.FileHandler(args.log_path, mode="w"),
                       logging.StreamHandler(sys.stdout)],
             datefmt='%Y-%m-%d %H:%M:%S')
    
        
    binarized_dfs = []
    
    # Process smFISH data (if provided)
    if args.smfish_data:
        logging.info(f"Preprocessing confocal smFISH data with threshold {args.smfish_thresh}")
        df_smfish = preprocess_registration_output(
            args.smfish_data, 
            log_dir / "qc_smfish_confocal", 
            visualize=True, 
            count_threshold=args.smfish_thresh, 
            technique="smfish"
        )
        # 🚨 DIAGNOSTIC INSERTION: SAVE THE PROCESSED SMFISH DF 🚨
        smfish_diagnostic_path = log_dir / "smfish_confocal_processed_diagnostic.tsv"
        df_smfish.to_csv(
            smfish_diagnostic_path, 
            sep='\t', 
            index=True, 
            index_label='label_id'
        )
        binarized_dfs.append(df_smfish)

    # Process prospr data (if provided)
    if args.prospr_data:
        logging.info(f"Preprocessing ProSPR data with threshold {args.prospr_thresh}")
        df_prospr = preprocess_registration_output(
            args.prospr_data, 
            log_dir / "qc_smfish_prospr", 
            visualize=True, 
            count_threshold=args.prospr_thresh, 
            technique="prospr"
        )
        binarized_dfs.append(df_prospr)
        
    # Combining the binarized tables
    
    if len(binarized_dfs) == 1:
        smfish_df_combined = binarized_dfs[0]
    else:
        # Concatenate and sum: sum > 0 means the gene is expressed in at least one method
        smfish_df_combined = pd.concat(binarized_dfs, axis=0).groupby(level=0).max()
        logging.info("Combined smFISH and ProSPeR data using OR logic.")

    output_combined_path = log_dir / "smfish_prospr_combined_binarized.tsv"
    smfish_df_combined.to_csv(
        output_combined_path, 
        sep='\t', 
        index=True, 
        index_label='label_id'
    )
    logging.info(f"Combined binarized expression saved to: {output_combined_path}")

    # preprocess scrnaseq

    logging.info(f"Quality control and preprocessing of the scRNAseq data...")
    scrnaseq_ad = preprocess_scrnaseq(
        args.scrnaseq_atlas, 
        args.gene_metadata, 
        args.cluster_column, 
        smfish_df_combined.columns, # Use combined columns for filtering
        log_dir / "qc_scrnaseq", 
        visualize=True
    )
    
    # Take only columns of smfish df which were found in scrnaseq atlas metadata
    # gene_list = scrnaseq_ad.var.pb_name                                               ############## CHANGED FOR TESTING!!!!!
    gene_list = scrnaseq_ad.var.pb_name.tolist() # Get the final scRNA-seq gene list
    
    """ # --- DIAGNOSTIC STEP 1: Save the definitive scRNA-seq gene list ---
    pd.Series(gene_list).to_csv(log_dir / "final_scrnaseq_gene_list.tsv", sep='\t', index=False, header=False)
    logging.info(f"Final scRNA-seq gene list saved to: {log_dir / 'final_scrnaseq_gene_list.tsv'}")
    
    # --- DIAGNOSTIC STEP 2: Save the smFISH combined columns before final filtering ---
    pd.Series(smfish_df_combined.columns).to_csv(log_dir / "smfish_combined_columns_unfiltered.tsv", sep='\t', index=False, header=False)
    logging.info(f"Unfiltered combined smFISH columns saved to: {log_dir / 'smfish_combined_columns_unfiltered.tsv'}") """

    smfish_df_combined = smfish_df_combined[gene_list]

    logging.info(f"Final gene list size: {len(gene_list)}")

    # Prepare the expression profiles
    profiles = get_cluster_profiles(scrnaseq_ad, args.cluster_column, log_dir / "profiles")
    total_profile = np.array(scrnaseq_ad.X.sum(axis=0) / scrnaseq_ad.X.shape[0]).flatten()
    visualize_profiles(gene_list, {**profiles, "total": total_profile}, log_dir / "profiles")
    

    # Given the binarized smFISH expression and expression profiles, assign cluster probability to each nucleus in EM
    assign_clusters(smfish_df_combined, gene_list, profiles, total_profile, log_dir / "cluster_prob")

    

if __name__=="__main__":
    main()