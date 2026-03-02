# platy_scrnaseq_integration
Map scRNAseq clusters to imaging data using smFISH and ProSPR signal with [scLocator](https://github.com/kreshuklab/scLocator).




## How to run

- Install [scLocator](https://github.com/kreshuklab/scLocator)
- Activate environment
```
conda activate sclocator_env
```
- Use configs to run mapping
```
sclocator --config /g/kreshuk/buglakova/projects/scLocator/tests/config_test_small_groups.yaml
```

**Final output**: 

`data/experiments/large_groups/results/cluster_probability.tsv` - this table can be viewed in Platybrowser if you open EM nuclei segmentation and add table from file system.


