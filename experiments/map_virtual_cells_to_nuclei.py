from platy_reg.n5_utils import read_volume
import matplotlib.pyplot as plt
from skimage.transform import resize
from pathlib import Path
from skimage.measure import regionprops_table
import pandas as pd
import numpy as np

def main():

    cell_segm = read_volume("/scratch/buglakova/data/prospr/sbem-6dpf-1-whole-segmented-cells.n5", "setup0/timepoint0/s4")
    prospr_image_dir = Path("/scratch/buglakova/data/prospr/prospr")
    prospr_image_paths = list(prospr_image_dir.glob("*"))
    print(len(prospr_image_paths))

    

    prospr_img = read_volume(prospr_image_paths[0], "setup0/timepoint0/s0")
    cell_segm = resize(cell_segm, prospr_img.shape, order=0)

    cell_ids = np.unique(cell_segm)[1:]
    merged_df_max = pd.DataFrame(index=cell_ids)
    merged_df_max.index.set_names("label_id", inplace=True)

    merged_df_mean = pd.DataFrame(index=cell_ids)
    merged_df_mean.index.set_names("label_id", inplace=True)
    for p in prospr_image_paths:
        if p not in merged_df_max.columns:
            prospr_img = read_volume(p, "setup0/timepoint0/s0")
            gene_name = p.stem.split("-")[-1].split(".")[0]
            print(gene_name)
            props = regionprops_table(cell_segm, intensity_image=prospr_img, properties=["label", "intensity_max", "intensity_mean"])
            data = pd.DataFrame(props)
            data = data.set_index("label")

            merged_df_max[gene_name] = data.loc[cell_ids, "intensity_max"]
            merged_df_max.to_csv("prospr_genes_max_per_cell.csv")

            merged_df_mean[gene_name] = data.loc[cell_ids, "intensity_mean"]
            merged_df_mean.to_csv("prospr_genes_mean_per_cell.csv")


if __name__=="__main__":
    main()