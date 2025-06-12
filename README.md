# platy_scrnaseq_integration
Map scRNAseq clusters to imaging data using smFISH signal

**Given** 
- scRNAseq atlas consisting of $G$ genes and $N_{cells}$ cells split into $N$ clusters
- Number of cells in cluster $i$: $|cluster_i|$
- For gene $k$ number of cells in which this gene is expressed is $nonzero_k \in [0; 1]$
- For gene $k$ number of cells in which this gene is expressed in cluster $i$ is $nonzero_{k,i} \in [0; 1]$

**Input**
- Cell with expression vector $\vec{e}$, for example, $(0, 0, 1, 1, 0)$, so binary measurements "expressed/ not expressed" for a subset of genes of the original atlas

**Output**

Probability that given the expression measurements the cell belongs to cluster $i$: $P(cell \in cluster_i | \vec{e}), i \in \{1 .. N\}$

**Method**

Using Bayes theorem
$$
P(cell \in cluster_i | \vec{e}) = \frac{P(\vec{e}|cell \in cluster_i)\cdot P(cell \in cluster_i)}{P(\vec{e})}
$$

All parts can be easily estimated from the scRNAseq atlas:

$$P(\vec{e}|cell \in cluster_i) = \prod\limits_{k}nonzero_{k,i}$$
$$P(cell \in cluster_i) = \frac{|cluster_i|}{N_{cells}}$$
$$P(\vec{e}) = \prod\limits_{k}nonzero_{k}$$