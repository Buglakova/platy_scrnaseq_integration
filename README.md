# platy_scrnaseq_integration
Map scRNAseq clusters to imaging data using smFISH signal

## Algorithm

**Given** 
- scRNAseq atlas consisting of $G$ genes and $N_{cells}$ cells split into $N$ clusters
- For gene $k$ number of cells in which this gene is expressed in cluster $i$ is $nonzero_{k,i} \in [0; 1]$

**Input**
- Cell in EM with expression vector $\vec{e}$, for example, $(0, 0, 1, 1, 0)$, so binary measurements "expressed/ not expressed" for a subset of genes of the original atlas. $\vec{e}_g$ is the binary number corresponding to the gene $g$.

**Output**

Probability that given the expression measurements the cell belongs to cluster $i$: $P(cell \in cluster_i | \vec{e}), i \in \{1 .. N\}$

**Method**

Using Bayes theorem
$$
P(cell \in cluster_i | \vec{e}) = \frac{P(\vec{e}|cell \in cluster_i)\cdot P(cell \in cluster_i)}{P(\vec{e})}
$$

All parts can be easily estimated from the scRNAseq atlas:

$$P(\vec{e}|cell \in cluster_i) = \prod\limits_{g=0}^{G}[\vec{e}_g \cdot nonzero_{g,i} + (1 -\vec{e}_g) \cdot (1 - nonzero_{g,i})]$$

$$P(cell \in cluster_i) = \frac{1}{N_{clusters}}$$
$$P(\vec{e}) = \sum\limits_{k=0}^{N_{clusters}}P(\vec{e}|cell \in cluster_k)\cdot P(cell \in cluster_k)$$


Result is a table of size $N_{cells} \times N_{clusters}$ where each element is a probability that given cell belongs to the given cluster.


To accomodate for the missing signal and autofluorescent areas, two artificial clusters are added:
- **zero**: expression profile with value $0.01$ for each gene, describes cells in which no expression was measured
- **autofluorescence$**:  expression profile with value $0.99$ for each gene, describes cells in which expression appears to be measured for all genes, indicating areas with some artifacts



## How to run

Create environment using `environment.yml`

Run:
```python experiments/run_mapping.py data/scrnaseq_atlas/scrnaseq_atlas.h5ad data/scrnaseq_atlas/gene_metadata.tsv group_full data/smfish/expression_matrix_aggregated.tsv data/experiments/large_groups/ data/experiments/large_groups/cluster_mapping.log --smfish_thresh 20 --scrnaseq_thresh 2 ```

Here:

`data/scrnaseq_atlas/scrnaseq_atlas.h5ad` - scRNAseq atlas object

`data/scrnaseq_atlas/gene_metadata.tsv` - metadata table from Cyril to convert gene names to XLOC

`group_full` - column in the scRNAseq atalas `obs` to use as cluster names

`data/smfish/expression_matrix_aggregated.tsv` - output of the registration pipeline

`data/experiments/large_groups/` - results (plots and tables) are saved in this directory

`data/experiments/large_groups/cluster_mapping.log` - log with some QC information

Optional arguments `--smfish_thresh` and `--scrnaseq_thresh` define the thresholds for binarizing each data type

**Final output**: `data/experiments/large_groups/cluster_prob/cluster_probability.tsv` - this table can be viewed in Platybrowser if you open EM nuclei segmentation and add table from file system.


## Tasks

- Add missing genes from smFISH to the `gene_metadata.tsv`
- Try with different binariation thresholds and clusters
- Implement using Prospr as additional evidence/to increase number of genes

