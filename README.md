# Spatial and single-nucleus transcriptomic analysis of genetic and sporadic forms of Alzheimer’s disease

[Miyoshi & Morabito et al. 2024 (Nature Genetics)](https://www.nature.com/articles/s41588-024-01961-x)

*The pathogenesis of Alzheimer’s disease (AD) depends on environmental and heritable factors, with its molecular etiology still unclear. Here we present a spatial transcriptomic (ST) and single-nucleus transcriptomic survey of late-onset sporadic AD and AD in Down syndrome (DSAD). Studying DSAD provides an opportunity to enhance our understanding of the AD transcriptome, potentially bridging the gap between genetic mouse models and sporadic AD. We identified transcriptomic changes that may underlie cortical layer-preferential pathology accumulation. Spatial co-expression network analyses revealed transient and regionally restricted disease processes, including a glial inflammatory program dysregulated in upper cortical layers and implicated in AD genetic risk and amyloid-associated processes. Cell–cell communication analysis further contextualized this gene program in dysregulated signaling networks. Finally, we generated ST data from an amyloid AD mouse model to identify cross-species amyloid-proximal transcriptomic changes with conformational context.*

**This repository contains the code used for data processing and analysis in our manuscript, and is generally organized in sync with the presentation of the data in the corresponding paper.** 

## Data generated in this study

The raw and processed ST (10X Genomics Visium) and snRNA-seq(Parse Biosciences) datasets have been deposited on the NCBI Gene Expression Omnibus (GEO) at accession number [GSE233208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208). Please contact the corresponding author of the paper (Vivek Swarup) with any queries related to the dataset.

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

### snRNA-seq differential expression 

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


## hdWGCNA co-expression network analysis (Fig. 3, Extended Data Fig. 9)

* [hdWGCNA in human ST](spatial/hdWGCNA/hdWGCNA_human.Rmd)
* [hdWGCNA in mouse ST](spatial/hdWGCNA/hdWGCNA_mouse.Rmd)

## Spatial and single-nucleus genetic enrichment analysis (Fig. 3, Extended Data Fig. 4)

* [Run scDRS in human snRNA datasets](snRNA/genetic_enrichment/run_scDRS/)
* [Plotting and downstream analysis of scDRS results in snRNA](snRNA/genetic_enrichment/scDRS_results_plotting.Rmd)
* [Run scDRS in human ST](spatial/genetic_enrichment/human_ST_scDRS.sub)
* [Plotting and downstream analysis of scDRS results in human ST](spatial/genetic_enrichment/scDRS_results_plotting_human.Rmd)
* [Run scDRS in mouse ST](spatial/genetic_enrichment/mouse_ST_scDRS.sub)
* [Plotting and downstream analysis of scDRS results in mouse ST](spatial/genetic_enrichment/scDRS_results_plotting_mouse.Rmd)

## Imaging mass cytometry (IMC) analysis (Fig. 5)

* [IMC clustering analysis and downstream plotting](IMC/hyperion_clustering.Rmd)

## Predicting spatial coordinates for snRNA-seq data (Extended Data Fig. 5)

* [Functions to run CellTrek in parallel for each snRNA-ST pair](snRNA/predict_spatial_coords/celltrek_parallel.R)
* [Running CellTrek in parallel for each snRNA-ST pair](snRNA/predict_spatial_coords/cellrtrek_parallel_run.sub)
* [Downstream plotting script pt. 1](snRNA/predict_spatial_coords/cellrtrek_downsrtream_plotting.Rmd)
* [Downstream plotting script pt. 2](snRNA/predict_spatial_coords/spatial_mapping_celltrek.Rmd)


## Cell-cell communication (CCC) network analysis (Fig. 6)

* [Cell-cell-communication network analysis with CellChat](snRNA/cell_cell_communication/spatial_cellchat.Rmd)
* [Cell-cell-communication network analysis with LIANA and comparing with CellChat](snRNA/cell_cell_communication/LIANA_ccc.Rmd)

## Amyloid-associated gene expression signatures (Fig. 7)

* [Integration of amyloid imaging with ST datasets](spatial/amyloid/splitseq_human_amyloid.Rmd)
* [Human amyloid analysis](spatial/amyloid/voyager_gene_signatures_human.Rmd)
* [Mouse amyloid analysis](spatial/amyloid/voyager_gene_signatures_5x.Rmd)

## TODO:
* Clean up code and analyis that we did which did not end up included in the paper.
