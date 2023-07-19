
# script to run DEG plotting in parallel
library(optparse)

option_list = list(
  make_option(
    c('--invis'), type='character', default=NULL,
    help='Path to the input Seurat object containing Visium data', metavar='character'
  ),
  make_option(
    c('--insc'), type='character', default=NULL,
    help='Path to the input Seurat object containing SC/SN data', metavar='character'
  ),
  make_option(
    c('--outdir'), type='character', default='./',
    help='Directory to place output files'
  ),
  make_option(
    c('--index'), type='numeric', default=1,
    help='Index from SLURM task array'
  ),
  make_option(
    c('--vis-sample-col'), type='character', default="Sample",
    help='Name of the column in seurat_vis@meta.data containing sample info.'
  ),
  make_option(
    c('--sc-sample-col'), type='character', default="Sample",
    help='Name of the column in seurat_obj@meta.data containing sample info.'
  ),
  make_option(
    c('--cell-id-col'), type='character', default="seurat_clusters",
    help='Name of the column in seurat_obj@meta.data containing cluster/celltype/annotation info for integration.'
  ),
  make_option(
    c('--iterations'), type='numeric', default=3,
    help='Number of times to run celltrek cell charting.'
  ),
  make_option(
    c('--overwrite'), type='character', default="FALSE",
    help='Whether or not to overwrite the output files.'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)


################################################################################
# Load R packages
################################################################################

library(Seurat)
library(Matrix)
library(CellTrek)
library(tictoc)

################################################################################
# setup variables based on optparse
################################################################################

# load seurat objects
print('Loading Seurat objects...')
vis_file <- as.character(opt[['invis']])
sc_file <- as.character(opt[['insc']])
seurat_vis <- readRDS(file=vis_file)
seurat_obj <- readRDS(file=sc_file)
print('Done!')

index <- opt[['index']]
cell_id_col <- opt[['cell-id-col']]
vis_sample_col <- opt[['vis-sample-col']]
sc_sample_col <- opt[['sc-sample-col']]
n_iterations <- opt[['iterations']]
data_dir <- opt[['outdir']]
overwrite <- as.logical(opt[['overwrite']])

# create output dir:
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}

################################################################################
# setup sample info and subset seurat objects for cell mapping
################################################################################

# list of vis & sc samples
vis_samples <- unique(as.character(seurat_vis@meta.data[[vis_sample_col]]))
sc_samples <- unique(as.character(seurat_obj@meta.data[[sc_sample_col]]))

# get index of the visium & sc samples
vis_sample_ind <- ceiling(index / length(sc_samples))
sc_sample_ind <- 1 + (index %% length(sc_samples))

# get the name of the visium & sc samples
cur_vis_sample <- vis_samples[vis_sample_ind]
cur_sc_sample <- sc_samples[sc_sample_ind]

# output name is the concatenation of the sample names:
outname <- paste0('vis_', cur_vis_sample, '-sc_', cur_sc_sample)

if(file.exists(paste0(data_dir, outname, '_coembed.rds')) & !overwrite){
  stop(paste0("Output file already exists and overwrite was set to FALSE."))
}

# get the current visium sample and make sure to set up the image
cur_vis <- seurat_vis[,seurat_vis@meta.data[[vis_sample_col]] == cur_vis_sample]
cur_image <- names(cur_vis@images)[sapply(names(cur_vis@images), function(x){nrow(cur_vis@images[[x]]@coordinates) > 0})]
cur_vis@images <- list(cur_image = cur_vis@images[[cur_image]])


# need to fix the image row:
cur_row <- cur_vis@images$cur_image@coordinates$row
row_min <- min(cur_row);
cur_row <- cur_row - row_min
row_max <- max(cur_row)
cur_row <- cur_row - row_max
cur_row <- cur_row * -1
cur_row <- cur_row + row_min

cur_vis$row <- cur_row
cur_vis@images$cur_image@coordinates$row <- cur_row
cur_vis$col <- cur_vis@images$cur_image@coordinates$col
cur_vis$imagerow <- cur_vis@images$cur_image@coordinates$imagerow
cur_vis$imagecol <- cur_vis@images$cur_image@coordinates$imagecol

# get the current sc seurat object
cur_seurat <- seurat_obj[,seurat_obj@meta.data[[sc_sample_col]] == cur_sc_sample]
cur_seurat$orig.ident <- cur_seurat@meta.data[[sc_sample_col]]

# add barcode column to vis & sc:
cur_vis$bc <- colnames(cur_vis)
cur_seurat$bc <- colnames(cur_seurat)


################################################################################
# Integration of sc + vis datasets:
################################################################################
print('Integrating datasets...')
ct_train <- CellTrek::traint(
  st_data=cur_vis,
  sc_data=cur_seurat,
  sc_assay='RNA',
  cell_names=cell_id_col
)
saveRDS(ct_train, file=paste0(data_dir, outname, '_coembed.rds'))
print('Done!')

################################################################################
# celltrek mapping iteratively
################################################################################

print('Celltrek mapping...')
obj_list <- list()
mapped_bcs <- c(); keep <- colnames(ct_train)
for(i in 1:n_iterations){
  print(i)
  tic()
  ct_obj <- CellTrek::celltrek(
    st_sc_int=ct_train[,keep],
    int_assay='traint',
    sc_assay = 'RNA',
    reduction='pca',
    intp=T, intp_pnt=10000,
    intp_lin=F, nPCs=30,
    ntree=1000,
    dist_thresh=0.75,
    top_spot=1, spot_n=5,
    repel_r=20,
    repel_iter=20,
    keep_model=T,
    verbose=TRUE
  )$celltrek
  time_elapsed <- toc()

  # which cells were not mapped??
  mapped_bcs <- c(mapped_bcs, unique(ct_obj$bc))
  print(length(mapped_bcs) / ncol(cur_seurat))

  keep <- ifelse(ct_train$type == 'st', TRUE, !(as.character(ct_train$bc) %in% mapped_bcs))
  #table(keep)

  obj_list[[i]] <- ct_obj

}

ct_coords <- do.call(rbind, lapply(obj_list, function(x){
  x@reductions$celltrek@cell.embeddings
}))

ix <- match(rownames(ct_coords), colnames(ct_train))
ct_meta <- ct_train@meta.data[ix,colnames(cur_seurat@meta.data)]
ct_meta <- cbind(ct_meta, ct_coords)

# only take the barcode and the coordinates to save:
ct_meta <- dplyr::select(ct_meta, c(bc, celltrek_1, celltrek_2))


write.csv(ct_meta, quote=FALSE, file=paste0(data_dir, outname, '_mapped_coords.csv'))
