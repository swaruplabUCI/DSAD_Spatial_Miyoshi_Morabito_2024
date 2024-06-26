# ADDS_paper

Miyoshi & Morabito et al. 2024 (Nature Genetics, in press)

This repository contains the code used for data processing and analysis in our 
manuscript titled "Spatial and single-nucleus trancsriptomic analysis of 
genetic and sporadic forms of Alzheimer's disease".

TODO: Add a brief intro 

## TODO:

* Clean up all un-used code 
* Provide basic comments in all scripts

## Data generated in this study

The raw and processed ST (10X Genomics Visium) and snRNA-seq
(Parse Biosciences) datasets have been deposited on the NCBI 
Gene Expression Omnibus (GEO) at accession number XXXXXX.

## Processing sequencing data and quantifying gene expression

## Spatial and single-nucleus clustering analysis (Fig. 1)

### snRNA-seq clustering analysis 

New snRNA-seq data
* [Clustering analysis for DSAD and Control data generated in this study](snRNA/clustering/clustering_scvi.ipynb)
* [Format dataset as a Seurat Object](snRNA/clustering/DSAD_create_seurat.Rmd)

Integration 
* [Integration with 3 published studies using scANVI](snRNA/clustering/integration_scanvi.ipynb)
* [Subclustering analysis of major cell lineages](snRNA/clustering/subclustering_scvi.ipynb)
* [Format integrated dataset as a Seurat Object](snRNA/clustering/integrated_create_seurat.Rmd)


* [Differential cell state abundance analysis with MiloR](snRNA/clustering/diff_cell_abundance_MiloR.Rmd)


### ST clustering analysis 

Human dataset
* [Load human spaceranger outputs and format as a Seurat Object](spatial/clustering/human_create_seurat.Rmd)
* [Format human ST dataset as an anndata object, pt1](spatial/clustering/seurat_to_anndata_human.Rmd)
* [Format human ST dataset as an anndata object, pt2](spatial/clustering/seurat_to_anndata_human.ipynb)
* [BayesSpace clustering analysis in the human dataset](spatial/clustering/clustering_human.Rmd)

Mouse dataset
* [Load mouse spaceranger outputs and format as a Seurat Object](spatial/clustering/mouse_create_seurat.Rmd)
* [Format mouse ST dataset as an anndata object, pt1](spatial/clustering/seurat_to_anndata_mouse.Rmd)
* [Format mouse ST dataset as an anndata object, pt2](spatial/clustering/seurat_to_anndata_mouse.ipynb)
* [BayesSpace clustering analysis in the mouse dataset](spatial/clustering/clustering_mouse.Rmd)

Additional plotting 
* [Additional plotting scripts](spatial/clustering/plotting_for_paper.Rmd)

## Differential expression analysis (Fig. 2)

Find the spatial condition DEGs from Emily's code and move it 
Also need to find or re-make the 5xFAD vs WT DEG scripts...

All snRNA-seq DEG scripts are there

## hdWGCNA co-expression network analysis (Fig. 3)

## Spatial and single-nucleus genetic enrichment analysis (Fig. 3, XXX)

TODO: get the .ipynb files from HPC3!!!

## Female vs Male differential expression (Fig. 4)

## Imaging mass cytometry (IMC) analysis (Fig. 5)

## Predicting spatial coordinates for snRNA-seq data (Extended Data Fig. 5)

## Cell-cell communication (CCC) network analysis (Fig. 6)

## Amyloid-associated gene expression signatures (Fig. 7)

