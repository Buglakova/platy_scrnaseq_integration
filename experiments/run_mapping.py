from cluster_map.preprocess_smfish import preprocess_registration_output
from cluster_map.preprocess_scrnaseq import preprocess_scrnaseq, get_cluster_profiles

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
        description="""Assign scRNAseq atlas cluster id to nuclei in EM.
        """
    )
    parser.add_argument("scrnaseq_atlas", type=str, help="Path of the scRNAseq atlas")
    parser.add_argument("gene_metadata", type=str, help="Correspondence between gene names and XLOCs")
    parser.add_argument("cluster_column", type=str, help="Column in obs to use as cluster ids")
    parser.add_argument("smfish_aggregated_expression", type=str, help="Path to the smFISH expression table: each row is a nucleus in EM")
    # parser.add_argument("segm_n5", type=str, help="Path of the EM segmentation for visualizations")
    # parser.add_argument("segm_n5_key", type=str, help="Key in the EM segmentation n5")
    parser.add_argument("log_dir", type=str, help="Directory to store diagnostic plots etc")
    parser.add_argument("log_path", type=str, help="Path to log file")

    args = parser.parse_args()
    
    log_dir = Path(args.log_dir)
    log_dir.mkdir(exist_ok=True)
    
    logging.basicConfig(level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[logging.FileHandler(args.log_path, mode="w"),
                        logging.StreamHandler(sys.stdout)],
            datefmt='%Y-%m-%d %H:%M:%S')  

    # segm = read_volume(args.segm_n5, args.segm_n5_key)
    # logging.info(f"EM segmentation shape {segm.shape}")
    # logging.info(f"Calculate nuclei props")
    # nuclei_props = pd.DataFrame(regionprops_table(segm, properties=("label", "area", 'centroid')))

    logging.info(f"Quality control and preprocessing of the smFISH data")
    smfish_df = preprocess_registration_output(args.smfish_aggregated_expression, log_dir / "qc_smfish", visualize=False)

    logging.info(f"Quality control and preprocessing of the scRNAseq data")
    scrnaseq_ad = preprocess_scrnaseq(args.scrnaseq_atlas, args.gene_metadata, args.cluster_column, smfish_df.columns, log_dir / "qc_scrnaseq", visualize=False)
    
    # Take only columns of smfish df which were found in scrnaseq atlas metadata

    gene_list = scrnaseq_ad.var.pb_name
    smfish_df = smfish_df[gene_list]

    logging.info("Final gene list")
    logging.info(gene_list)

    # Prepare the expression profiles

    profiles = get_cluster_profiles(scrnaseq_ad, args.cluster_column, log_dir / "profiles")




    

if __name__=="__main__":
    main()