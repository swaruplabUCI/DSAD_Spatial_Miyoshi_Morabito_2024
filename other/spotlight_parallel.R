


library(optparse)

option_list = list(
  make_option(
    c('--seurat-vis'), type='character', default=NULL,
    help='Seurat object .rds file for Visium data', metavar='character'
  ),
  make_option(
    c('--seurat-sc'), type='character', default=NULL,
    help='Seurat object .rds file for single-cell data', metavar='character'
  ),
  make_option(
    c('--markers'), type='character', default=NULL,
    help='DEG results as a delimited file (.csv by default)', metavar='character'
  ),
  make_option(
    c('--outdir'), type='character', default='./',
    help='Directory to place output files'
  ),
  make_option(
    c('--index'), type='numeric', default=NULL,
    help='SLURM task array number goes here, selects which Visium sample to process on this task'
  ),
  make_option(
    c('--cluster-name'), type='character', default=NULL,
    help='name of the cluster / cell group column in the single-cell Seurat object corresponding to the marker DEGs'
  ),
  make_option(
    c('--sample-col'), type='character', default="Sample",
    help='name of the Sample column in the visium seurat object'
  ),
  make_option(
    c('--topn'), type='numeric', default=100,
    help='How many marker genes to use for each group?'
  ),
  make_option(
    c('--logfc'), type='numeric', default=NULL,
    help='Log FC cutoff for marker genes to use.'
  ),
  make_option(
    c('--n-hvg'), type='numeric', default=3000,
    help='Number of HVGs to use for SPOTLight.'
  ),
  make_option(
    c('--n-cells'), type='numeric', default=100,
    help='Number of cells per group to use.'
  ),
  make_option(
    c('--seed'), type='numeric', default=12345,
    help='Seed for random numbers.'
  )
)


# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

################################################################################
# setup variables based on optparse
################################################################################

seurat_sc_file = as.character(opt["seurat-sc"])
seurat_vis_file = as.character(opt["seurat-vis"])
marker_file = as.character(opt$markers)
outdir = as.character(opt$outdir)
cluster_name = as.character(opt["cluster-name"])
top_n_markers = as.numeric(opt$topn)
logfc_cutoff = as.numeric(opt$logfc)
index = as.numeric(opt$index)
sample_col = as.character(opt['sample-col'])
rng_seed = as.numeric(opt$seed)
n_cells = as.numeric(opt['n-cells'])
n_hvg = as.numeric(opt['n-hvg'])

################################################################################
# Load required libraries
################################################################################

library(Seurat)
library(tidyverse)
library(Matrix)
library(SPOTlight)

# set random seed for reproducibility
set.seed(rng_seed)

################################################################################
# Load Seurat objects
################################################################################

print(seurat_vis_file)
print(seurat_sc_file)

seurat_vis <- readRDS(seurat_vis_file)
seurat_obj <- readRDS(seurat_sc_file)

print('Loaded Seurat objects')

################################################################################
# Load the markers, and filter based on criteria:
################################################################################

if(is.null(top_n_markers) & is.null(logfc_cutoff)){
  stop("Need to make a selection for --logfc or --topn")
}

# handle case when logFC is not specified
if(is.null(logfc_cutoff)){
  logfc_cutoff <- 0
}

# handle case when top markers is null
if(is.null(top_n_markers)){
  top_n_markers <- nrow(seurat_obj)
}

# load the file
markers <- read.csv(marker_file)
print('Loaded Marker genes')

# check if there is a column called group, and rename it to cluster
if("group" %in% colnames(markers)){
  markers <- markers %>% dplyr::rename(cluster = group)
}

# apply filters:
markers <- markers %>%
  subset(p_val_adj <= 0.05 & avg_log2FC >= logfc_cutoff) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(top_n_markers, wt=avg_log2FC)

print('Filtered marker genes')

################################################################################
# Set up visium seurat obj for selected sample
################################################################################

samples <- unique(seurat_vis@meta.data[[sample_col]])
cur_sample <- samples[index]
print(paste0("Subsetting sample ", cur_sample))

sp_subset <- seurat_vis[,seurat_vis@meta.data[[sample_col]] == cur_sample]
cur_image <- names(sp_subset@images)[sapply(names(sp_subset@images), function(x){nrow(sp_subset@images[[x]]@coordinates) > 0})]
sp_subset@images <- list(cur_image = sp_subset@images[[cur_image]])


################################################################################
# Run spotlight
################################################################################

print("Running SPOTLight")

spotlight_ls <- spotlight_deconvolution(
  se_sc = seurat_obj,
  counts_spatial = sp_subset@assays$Spatial@counts,
  clust_vr = cluster_name, # Variable in sc_seu containing the cell-type annotation
  cluster_markers = markers, # Dataframe with the marker genes
  cl_n = n_cells, # number of cells per cell type to use
  hvg = n_hvg, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
)


# get outputs
nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

decon_mtrx[, colnames(decon_mtrx) != "res_ss"] %>% dim

#This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(sp_subset)



################################################################################
# Run spotlight
################################################################################

# save whole spotlight output:
save(spotlight_ls,  file=paste0(outdir, '/', cur_sample, '_spotlight.rda'))

# save the decon table:
write.csv(as.data.frame(decon_mtrx), file=paste0(outdir, '/', cur_sample, '_spotlight.csv'), quote=FALSE)
