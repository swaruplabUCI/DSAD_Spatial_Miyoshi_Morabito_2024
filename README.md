# ADDS_paper

Miyoshi & Morabito et al. 2024 (Nature Genetics, in press)

This repository contains the code used for data processing and analysis in our 
manuscript titled "Spatial and single-nucleus trancsriptomic analysis of 
genetic and sporadic forms of Alzheimer's disease".

TODO: Add a brief intro 

## TODO:

* Clean up all un-used code 
* Provide basic comments in all scripts
* In each README section, write a short para or few sentences about the analysis.

## Data generated in this study

The raw and processed ST (10X Genomics Visium) and snRNA-seq
(Parse Biosciences) datasets have been deposited on the NCBI 
Gene Expression Omnibus (GEO) at accession number XXXXXX.

## Processing sequencing data and quantifying gene expression

* [Parse Biosciences snRNA-seq preprocessing scripts](snRNA/preprocessing/)
* [Human 10X Visium ST preprocessing scripts](spatial/preprocessing/human/)
* [Mouse 10X Visium ST preprocessing scripts](spatial/preprocessing/mouse/)

## Spatial and single-nucleus clustering analysis (Fig. 1)

### snRNA-seq clustering analysis 

New snRNA-seq data
* [Clustering analysis for DSAD and Control data generated in this study](snRNA/clustering/clustering_scvi.ipynb)
* [Format dataset as a Seurat Object](snRNA/clustering/DSAD_create_seurat.Rmd)

Integration 
* [Integration with 3 published studies using scANVI](snRNA/clustering/integration_scanvi.ipynb)
* [Subclustering analysis of major cell lineages](snRNA/clustering/subclustering_scvi.ipynb)
* [Format integrated dataset as a Seurat Object](snRNA/clustering/integrated_create_seurat.Rmd)


Differential cell state analysis
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

## Differential expression analysis (Figs. 2 and 4)

NoteL:Need to find from Emily's codebase:

* Human ST marker genes 
* Human ST condition 
* Mouse ST marker genes (?)

All snRNA-seq DEG scripts are there

### snRNA-seq differential expression 

TODO: Split the Sex analysis by snRNA vs ST ?

* [Marker genes](snRNA/differential_expression/markers/)
* [Disease vs Control DEG scripts](snRNA/differential_expression/condition/)
* [Disease vs Control plotting and downstream analysis](snRNA/differential_expression/condition_DEG_analysis_snRNA.Rmd)
* [Female vs Male plotting and downstream analysis](spatial/differential_expression/human/sex_DEGs_human_ST.Rmd)

### ST differential expression

Human dataset 
* [Cluster marker genes](spatial/differential_expression/human/cluster_markers_human.Rmd)
* [Disease vs Control DEG scripts](spatial/differential_expression/human/condition/)
* [Disease vs Control plotting and downstream analysis](spatial/differential_expression/condition_DEGs_ST.Rmd)
* [Female vs Male DEG scripts](spatial/differential_expression/human/sex/)
* [Female vs Male plotting and downstream analysis](spatial/differential_expression/human/sex_DEGs_human_ST.Rmd)

Mouse dataset 
* [Cluster marker genes](spatial/differential_expression/mouse/cluster_markers_mouse.Rmd)
* [5xFAD vs WT DEG scripts](spatial/differential_expression/mouse/condition)
* [5xFAD vs WT plotting and downstream analysis](spatial/differential_expression/condition_DEGs_ST.Rmd)
* [Female vs Male DEG scripts](spatial/differential_expression/mouse/sex/)
* [Female vs Male plotting and downstream analysis](spatial/differential_expression/mouse/sex_DEGs_mouse_ST.Rmd)


## hdWGCNA co-expression network analysis (Fig. 3)

* [hdWGCNA in human ST](spatial/hdWGCNA/hdWGCNA_human.Rmd)
* [hdWGCNA in mouse ST](spatial/hdWGCNA/hdWGCNA_mouse.Rmd)

## Spatial and single-nucleus genetic enrichment analysis (Fig. 3, XXX)

* [Run scDRS in human snRNA datasets](snRNA/genetic_enrichment/run_scDRS/)
* [Plotting and downstream analysis of scDRS results in snRNA](snRNA/genetic_enrichment/scDRS_results_plotting.Rmd)
* [Run scDRS in human ST](spatial/genetic_enrichment/human_ST_scDRS.sub)
* [Plotting and downstream analysis of scDRS results in human ST](spatial/genetic_enrichment/scDRS_results_plotting_human.Rmd)
* [Run scDRS in mouse ST](spatial/genetic_enrichment/mouse_ST_scDRS.sub)
* [Plotting and downstream analysis of scDRS results in mouse ST](spatial/genetic_enrichment/scDRS_results_plotting_mouse.Rmd)

## Imaging mass cytometry (IMC) analysis (Fig. 5)

* [IMC clustering analysis and downstream plotting](IMC/hyperion_clustering.Rmd)

## Predicting spatial coordinates for snRNA-seq data (Extended Data Fig. 5)

## Cell-cell communication (CCC) network analysis (Fig. 6)

## Amyloid-associated gene expression signatures (Fig. 7)

