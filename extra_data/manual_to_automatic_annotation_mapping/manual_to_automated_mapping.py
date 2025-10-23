import pandas as pd
import re

def create_clade_mapping_file():
    """
    Creates the family name to clade mapping file from the provided string.
    This is a helper function to ensure the script is fully self-contained.
    """
    mapping_string = "brain_5HT_IN_Phox2_Vsxclade11sub20subsub5brain_ACh_MN_mnx_phox2clade11sub37brain_ACh_MN_mnx_phox2_Lhx15clade11sub51brain_ACh_SSN_bsx_Dlxclade11sub48brain_DA_IN_Emx_Six12clade10sub7subsub4brain_GABA_IN_Dbx_Ptf1aclade11sub14subsub9brain_Glu_IN_TAL_lhx3clade11sub21brain_Glu_LLE-PRC2_FEZF_AP2clade11sub45subsub0brain_Glu_mechSN_FEZF_PouIV_Vsx_Nkx2-1clade10sub9brain_ls_IN_TAL_GATA_PA258clade11sub34brain_ls_pyg_ACh_SN_Phox2_HANDclade11sub22subsub4brain_or_ls_pyg_Glu_mechSN_POUVI_Lhx3clade10sub3subsub5-9brain_SSN_bsx_COEclade11sub3subsub9brain_SSN_six4_PDF_G0-Rclade11sub60brain_TYR/DA/ACh_NSC_foxA_nk2.1clade11sub16subsub7commisural_MNs_Lbx_Evxclade11sub6subsub12commisural_MNs_Lbx_Evxclade11sub11fg_GABA_MN_lmx1_sim_mnxclade11sub66fg_GABA_SN_Dbx_Ptf1aclade11sub28subsub3hg_Glu_EN_nkx22_lmx1_mnxnocladesub16ls_ACh_MN_Mnx_Lhx3_Pitxclade11sub29subsub5ls_GABA_cSN_Dbx_Ptf1aclade11sub5subsub10ls_HIS_VSN_foxQ2_phox2clade11sub6subsub11ls_pyg_Glu_cIN_Evxclade11sub13subsub5ls_pyg_Glu_cIN_Evx_Pou4clade11sub13subsub6ls1_5HT_MN_Pitx_GATA123_lhx15clade11sub20subsub0-1-2-4-6ls2-3_5HT_MNclade11sub20subsub3pyg_Glu_SN_POU4_BarH1_Islclade11sub6subsub6"
    
    # Use a regular expression to robustly parse the mapping string
    pattern = re.compile(r'(.+?)(clade.+?|noclade.+?)(?=(brain|commisural|fg|hg|ls|pyg|$))')
    matches = pattern.findall(mapping_string)
    
    with open('family_name_to_clade.tsv', 'w') as f:
        f.write("curated_cell_type\tmanual_clade\n")
        for name, clade, _ in matches:
            # Handle specific naming inconsistencies if necessary
            if "brain_TYR/DA/ACh_NSC_foxA_nk2.1" in name:
                name = name.replace("brain_TYR/DA/ACh_NSC_foxA_nk2.1", "brain_TYR-DA-ACh_NSC_foxA_nk2-1")
            
            f.write(f"{name}\t{clade}\n")
    print("Helper file 'family_name_to_clade.tsv' created.")


def main():
    """
    Merges manual curation with nucleus mapping, clade information, and pipeline results.
    """
    try:
        # Step 1: Load all three required files
        manual_curation_df = pd.read_csv('other_experiments/cell-types-manual-curation.tsv', sep='\t')
        cells_to_nuclei_df = pd.read_csv('other_experiments/cells_to_nuclei.tsv', sep='\t')
        clade_mapping_df = pd.read_csv('family_name_to_clade.tsv', sep='\t')
        cluster_prob_df = pd.read_csv('data/experiments/small_groups/cluster_prob/cluster_probability.tsv', sep='\t')

    except FileNotFoundError as e:
        print(f"❌ Error: Could not find the file '{e.filename}'.")
        print("Please make sure all necessary TSV files are in the same directory as this script.")
        return

    # Step 2: Merge the dataframes
    # Link manual cell curations to nucleus IDs
    merged_df = pd.merge(manual_curation_df, cells_to_nuclei_df, on='label_id', how='left')

    # Add the corresponding clade names from the mapping file
    final_df = pd.merge(merged_df, clade_mapping_df, on='curated_cell_type', how='left')
    
    # Rename columns for clarity and merge with the pipeline's final assignments
    final_df = final_df.rename(columns={'label_id': 'cell_id'})
    cluster_prob_df = cluster_prob_df.rename(columns={'label_id': 'nucleus_id'})
    
    # Join the manual data with the pipeline's output
    comparison_df = pd.merge(
        final_df,
        cluster_prob_df[['nucleus_id', 'most probable cluster']],
        on='nucleus_id',
        how='left'
    )

    # Step 3: Clean the final DataFrame
    # Drop any rows where a mapping was incomplete (e.g., no nucleus or no clade)
    comparison_df.dropna(subset=['nucleus_id', 'manual_clade', 'most probable cluster'], inplace=True)
    comparison_df = comparison_df[comparison_df['nucleus_id'] != 0]

    # Convert ID columns to clean integers
    comparison_df = comparison_df.astype({'nucleus_id': int, 'cell_id': int})

    # Step 4: Determine which assignments match
    comparison_df['match_status'] = (comparison_df['manual_clade'] == comparison_df['most probable cluster'])
    
    # Reorder columns for the final output file
    final_output = comparison_df[['nucleus_id', 'cell_id', 'curated_cell_type', 'manual_clade', 'most probable cluster', 'match_status']]

    # Step 5: Save the complete comparison to a file
    output_filename = 'comparison_results.tsv'
    final_output.to_csv(output_filename, sep='\t', index=False)

    print(f"\n✅ Success! The final comparison has been saved as '{output_filename}'.")
    
    # --- Print a summary of the results ---
    matches = final_output['match_status'].sum()
    mismatches = len(final_output) - matches
    
    print("\n--- Analysis Summary ---")
    print(f"Total Nuclei Compared: {len(final_output)}")
    print(f"Correct Matches: {matches}")
    print(f"Mismatches: {mismatches}")

if __name__ == "__main__":
    # Create the mapping file first, then run the main process
    create_clade_mapping_file()
    main()