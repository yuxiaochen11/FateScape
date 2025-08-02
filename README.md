# FateScape

<div style="float: right; margin-left: 15px;">
  <img src="man/logo.png" alt="FateScape logo" width="150">
</div>

FateScape leverages multimodal single-cell sequencing to infer cell division trees using a maximum likelihood framework. By combining CRISPR/Cas9-based lineage tracing technology with single-cell transcriptome sequencing data, FateScape can reconstruct cell division history and infer the plasticity and inheritance of cell states.


## Introduction

This repository provides a comprehensive suite of R functions for reconstructing cell division trees. The tools in this repository enable users to:

-   Refine lineage barcodes by addressing stochastic dropout issues.
-   Construct sub-cell division trees along different state lineages.
-   Optimize sub-cell division trees using a maximum likelihood framework combined with subtree exchanges.
-   Integrate sub-cell division trees into a complete cell division tree.
-   Infer the plasticity and heritability of cell states.

These functions are designed to help researchers analyze single-cell lineage data, explore cell lineage evolution, and gain insights into cell division patterns.

##  Installation

To install and use these functions, follow these steps:

``` bash
devtools::install_github(‘yuxiaochen11/FateScape')
```

##  Quick start

If you have a object with paired single-cell gene expression and lineage barcodes data, you can start right away with inferring the regulatory network:

``` r
# Load required libraries
library(FateScape)
library(PATH)

# Load scRNA-seq data
load("RNA_object.RData")

# Load lineage barcode data
load("barcode_object.RData")

# psedo-trajectory inference
umap_res <- Embeddings(RNA_object, reduction = "umap")
labels_new <- as.character(RNA_object$cell_states)

sce <- slingshot(data = umap_res,clusterLabels = labels_new, dist.method = "simple",start.clus = "Ectoderm")
state_lineages <- slingLineages(sce)

# Barcode imputation
barcodes <- imputation_dropout_alter(as.matrix(barcode_object), ...)

#Construction of the  sub-cell division trees 
Trees_initial <- initial_tree_construction(state_lineages, barcodes)
bestsubtree <- subtree_refinement(Trees_initial, ...)

#Removing duplicated leaf nodes across sub-cell division trees
bestsubtree <- drop_duplicated_tips(bestsubtree,...)

#Decomposition and reassembly of sub-cell division trees
subtrees_rootbar <- get_subtree_root_barcodes(bestsubtree, ...)
decomposed_subtrees <- decompose_subtrees(bestsubtree, ...)

# Integration of the sub-cell division trees
cell_division_tree <- merge_subcell_trees_ward(subtrees_rootbar, decomposed_subtrees)
```
## More
More info about FateScape can be found on our [website](https://yuxiaochen11.github.io/FateScape/). There you can find an API reference and a number of tutorial vignettes that give an introduction on how to use FateScape most effectively.


## Citation

If you use **FateScape** in your research, please cite the accompanying paper:

> Yu, Xiaochen (2025). _FateScape: Reconstructing cell division history and phenotypic dynamics from single-cell barcode and transcriptomic data_. GitHub. https://github.com/yuxiaochen11/FateScape

Or use the following BibTeX entry:

```bibtex
@software{FateScape,
  author  = {Yu, Xiaochen},
  title   = {FateScape: Reconstructing cell division history and phenotypic dynamics from single-cell barcode and transcriptomic data},
  year    = {2025},
  url     = {https://github.com/yuxiaochen11/FateScape},
}
```
