import pandas as pd
import logging
from pathlib import Path
import sys

# Setup basic logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)],
                    datefmt='%Y-%m-%d %H:%M:%S')

def integrate_prospr_data(genes_path, mapping_path, output_path):
    """
    Merges ProSPeR gene data with nuclei mapping, filters out cells without nuclei,
    renames the nucleus_id to label_id, and reorders the columns.
    """
    
    # Define column names
    PROSPR_CELL_COL = 'label_id'      # Cell ID column in the genes file
    MAPPING_CELL_COL = 'label_id'     # Cell ID column in the mapping file
    NUCLEI_COL = 'nucleus_id'       # Nucleus ID column in the mapping file
    FINAL_LABEL_COL = 'label_id'      # Target name for the nucleus ID column
    
    # Define a temporary name for the conflicting column during merge
    MAPPING_CELL_COL_TEMP = 'cell_id_for_merge'

    logging.info(f"Loading ProSPeR gene data from {genes_path}...")
    try:
        # Load the gene expression data
        df_genes = pd.read_csv(genes_path, sep='\t')
        
        # Temporarily rename the gene file's 'label_id' to avoid conflict during merge
        df_genes = df_genes.rename(columns={PROSPR_CELL_COL: 'original_cell_id'})
        
    except Exception as e:
        logging.error(f"Error loading genes file: {e}")
        return

    logging.info(f"Loading cell-to-nuclei mapping from {mapping_path}...")
    try:
        # Load the mapping file
        df_mapping = pd.read_csv(mapping_path, sep='\t')
        
        # 🚨 FIX: Rename the mapping's cell ID column to a temporary name BEFORE MERGE 🚨
        df_mapping = df_mapping.rename(columns={MAPPING_CELL_COL: MAPPING_CELL_COL_TEMP})
        
    except Exception as e:
        logging.error(f"Error loading mapping file: {e}")
        return

    # --- 1. PREP & MERGE ---
    logging.info("Preparing columns for merge and merging ProSPeR data with nuclei mapping...")
    
    # Ensure merge columns have the same type
    df_genes['original_cell_id'] = df_genes['original_cell_id'].astype(str)
    df_mapping[MAPPING_CELL_COL_TEMP] = df_mapping[MAPPING_CELL_COL_TEMP].astype(str)
    
    # Merge using the temporary cell ID column
    df_merged = df_genes.merge(
        df_mapping[[MAPPING_CELL_COL_TEMP, NUCLEI_COL]],
        left_on='original_cell_id',
        right_on=MAPPING_CELL_COL_TEMP,
        how='inner' 
    )
    
    # --- 2. FILTER & RENAME ---
    
    # Drop the temporary and redundant cell columns FIRST
    df_filtered = df_merged.drop(columns=['original_cell_id', MAPPING_CELL_COL_TEMP]) 
    
    logging.info("Filtering out cells not assigned to a nucleus (nucleus_id = 0)...")
    
    # Filter: Keep only rows where 'nucleus_id' > 0.
    df_filtered = df_filtered[pd.to_numeric(df_filtered[NUCLEI_COL], errors='coerce') > 0].copy()

    # 🚨 FIX: Now it is SAFE to rename 'nucleus_id' to 'label_id' 🚨
    df_filtered = df_filtered.rename(columns={NUCLEI_COL: FINAL_LABEL_COL})
    
    
    # --- 3. REORDER & SAVE ---
    logging.info("Reordering columns to place the new 'label_id' (nucleus ID) first...")
    
    # Get all column names
    cols = df_filtered.columns.tolist()
    
    # Move 'label_id' to the beginning of the list
    cols.remove(FINAL_LABEL_COL)
    cols.insert(0, FINAL_LABEL_COL)
    
    # Apply the new column order
    df_final = df_filtered[cols]

    logging.info(f"Original cells: {len(df_genes)}")
    logging.info(f"Cells assigned to a nucleus: {len(df_final)}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(output_path, sep='\t', index=False)
    
    logging.info(f"Successfully created and saved merged file to {output_path}")


# --- Execute Script ---
if __name__ == "__main__":
    # Define the file paths (ensure these match the files in your directory)
    GENES_FILE = Path("other_data/smfish/prospr_genes.tsv")
    MAPPING_FILE = Path("other_data/smfish/cells_to_nuclei.tsv")
    OUTPUT_FILE = Path("other_data/smfish/prospr_genes_nuclei.tsv")
    
    integrate_prospr_data(GENES_FILE, MAPPING_FILE, OUTPUT_FILE)