# platy_scrnaseq_integration
Map scRNAseq clusters to imaging data using smFISH and ProSPR signal with [scLocator](https://github.com/kreshuklab/scLocator).




## How to run

- Install [scLocator](https://github.com/kreshuklab/scLocator)
- Activate environment
```
conda activate sclocator_env
```
- Put the full scRNAseq `h5ad` (or symlink to it) to `data/scrnaseq_atlas/atlas_6dpf.h5ad`. Otherwise, change the corresponding path in the config. 
- Use configs to run mapping

for broad types:
```
sclocator --config experiments/broad_types/config_broad_types.yaml
```

for high-resolution annotations:
```
sclocator --config experiments/curated_high_res_annotations/config_highres_annotations.yaml
```

**Final output**: 

`data/experiments/large_groups/results/cluster_probability.tsv` - this table can be viewed in Platybrowser if you open EM nuclei segmentation and add table from file system.


