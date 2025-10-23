import pandas as pd
import numpy as np
import logging
import sys
import argparse
import pickle
from pathlib import Path 

# --- Assuming get_individual_gene_probabilities is defined/imported from cluster_map.bayesian_prob ---
# For independence, we assume the user has access to this function,
# but it is not defined in this file to keep focus.
from cluster_map.bayesian_prob import get_individual_gene_probabilities


# Setup logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)],
                    datefmt='%Y-%m-%d %H:%M:%S')

# ----------------------------------------------------------------------
# Helper Functions for Data Loading
# ----------------------------------------------------------------------

def load_smfish_expr(path: Path, nucleus_id: float) -> pd.DataFrame:
    """Loads the binarized smFISH expression data and ensures the index is float."""
    try:
        df = pd.read_csv(path, sep='\t', index_col='label_id')
        df.index = df.index.astype(float)

        logging.critical(f"Loaded SMFISH Data has {len(df.columns)} columns.")
        logging.critical(f"--- 215 Gene List Start ---")
        for i, col in enumerate(df.columns):
            if i < 215: # Print the first few columns
                print(f"Column {i+1}: {col}")
            else:
                print("... (Truncated)")
                break
        logging.critical(f"--- 215 Gene List End ---")


        # Ensure the nucleus exists after potential data issues
        if nucleus_id not in df.index:
            raise KeyError(f"Nucleus ID {nucleus_id} not found in smFISH data after loading.")
        return df
    except Exception as e:
        logging.error(f"Failed to load smFISH expression file ({path.name}): {e}")
        sys.exit(1)

def load_assignments(path: Path) -> pd.DataFrame:
    """Loads cluster probability results to find the assigned cluster."""
    try:
        df = pd.read_csv(path, sep=r'\s+', index_col='label_id')
        df.index = df.index.astype(float)
        # We need the assignment column for the lookup
        if 'most probable cluster' not in df.columns:
            # Handle the case where the header split 'most probable cluster'
            df = df.rename(columns={'cluster': 'most probable cluster'}) 
        return df
    except Exception as e:
        logging.error(f"Failed to load cluster assignment file ({path.name}): {e}")
        sys.exit(1)

def load_profiles(profiles_dir: Path) -> dict:
    """Loads the pickled cluster profiles dictionary."""
    profile_path = profiles_dir / "profiles_dictionary.pkl"
    try:
        with open(profile_path, 'rb') as f:
            profiles = pickle.load(f)
        return profiles
    except Exception as e:
        logging.error(f"Failed to load cluster profiles file: {e}")
        sys.exit(1)

# ----------------------------------------------------------------------
# TARGETED DIAGNOSTIC FUNCTION (Provided by user, integrated here)
# ----------------------------------------------------------------------

def run_targeted_gene_likelihood_diagnostic(
    nucleus_id: float, 
    smfish_df: pd.DataFrame, 
    profiles: dict, 
    gene_list: list, 
    assigned_cluster: str, 
    target_cluster_list: list, 
    log_dir: Path
):
    """
    Calculates the gene likelihood for a single nucleus against specified target clusters 
    and saves the detailed output for manual review.
    """
    
    # 1. Select the cell's expression vector
    # Access is guaranteed safe by the outer validation loop
    cell_expr = smfish_df.loc[nucleus_id].to_numpy()
    
    df_diagnostics = []
    
    # Ensure the list includes the assigned cluster for context
    clusters_to_test = list(set(target_cluster_list + [assigned_cluster]))

    # 2. Calculate the likelihood vector for this cell against *each* target cluster
    for cluster_name in clusters_to_test:
        if cluster_name not in profiles:
            logging.warning(f"Profile for target cluster '{cluster_name}' not found. Skipping.")
            continue
            
        profile = profiles[cluster_name]
        gene_likelihoods = get_individual_gene_probabilities(cell_expr, profile)

        df_row = pd.DataFrame({
            'label_id': nucleus_id,
            'assigned_cluster': assigned_cluster,
            'tested_cluster': cluster_name,
            'gene_name': gene_list, 
            'likelihood_term': gene_likelihoods
        })
        df_diagnostics.append(df_row)

    # 3. Save the targeted file
    if df_diagnostics:
        df_final = pd.concat(df_diagnostics, ignore_index=True)
        
        # Ensure nucleus_id is safely converted to integer for filename
        safe_nucleus_id = int(nucleus_id) if nucleus_id == int(nucleus_id) else str(nucleus_id).replace('.', '_')
        
        output_path = log_dir / "targeted_diagnostics" / f"cell_{safe_nucleus_id}_vs_targets.tsv"
        output_path.parent.mkdir(parents=True, exist_ok=True) # Directory creation fix
        
        df_final.to_csv(output_path, sep='\t', index=False)
        logging.info(f"Targeted diagnostic saved for cell {nucleus_id} to {output_path}")

# ----------------------------------------------------------------------
# Main Execution Block
# ----------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnostic plot for single-nucleus cluster assignment.")
    parser.add_argument("nucleus_id", type=float, help="The ID of the nucleus to plot (e.g., 7264.0).")
    parser.add_argument("base_dir", type=str, help="Base directory where mapping outputs are stored (e.g., data/experiments).")
    parser.add_argument("--compare_to", nargs='+', default=[], help="List of specific cluster names to include for comparison (e.g., Midgut 'Nervous System').")

    args = parser.parse_args()

    BASE_DIR = Path(args.base_dir)
    TARGET_ID = args.nucleus_id
    
    # --- Define File Paths ---
    SMFISH_EXPR_PATH = BASE_DIR / "smfish_prospr_combined_filtered.tsv"
    ASSIGNMENT_PATH = BASE_DIR / "cluster_prob" / "cluster_probability.tsv"
    PROFILES_DIR = BASE_DIR / "profiles" 
    
    # --- Load Necessary Data ---
    df_smfish = load_smfish_expr(SMFISH_EXPR_PATH, TARGET_ID)
    df_assignments = load_assignments(ASSIGNMENT_PATH)
    profiles = load_profiles(PROFILES_DIR)
    
    # --- Derive Required Variables ---
    
    # 1. Get Assigned Cluster ID
    if TARGET_ID not in df_assignments.index:
        logging.error(f"Target Nucleus ID {TARGET_ID} not found in the assignment file.")
        sys.exit(1)
        
    ASSIGNED_CLUSTER = df_assignments.loc[TARGET_ID, 'most probable cluster']
    
    # 2. Get Gene List (from the first profile array)
    # The gene list should align with the order of the profiles.
    # We must ensure the smfish data is correctly filtered before passing its columns.
    
    # Use the keys from the smfish_df to get the gene list
    GENE_LIST = df_smfish.columns.tolist() 

    # --- Run Targeted Diagnostic ---
    run_targeted_gene_likelihood_diagnostic(
        TARGET_ID,
        df_smfish,
        profiles,
        GENE_LIST,
        ASSIGNED_CLUSTER,
        args.compare_to,
        BASE_DIR
    )