
library(optparse)

option_list = list(
  make_option(
    c('--seurat'), type='character', default=NULL,
    help='Cell Dataset .rds file', metavar='character'
  ),
  make_option(
    c('--outdir'), type='character', default='./',
    help='Directory to place output files', metavar='character'
  ),
  make_option(
    c('--type'), type='character', default='markers',
    help='Which type of test to run? Choose markers or condition. Markers does 1 vs all test, condition compares 2 conditions in a certain cluster.', metavar='character'
  ),
  make_option(
    c('--cluster'), type='character', default='seurat_clusters',
    help='Column in seurat metadata that indicates clusters'
  ),
  make_option(
    c('--name'), type='character', default='DEGs',
    help='Name to append to output files.'
  ),
  make_option(
    c('--index'), type='numeric', default=NULL,
    help='SLURM task array number goes here, selects which group to process on this task'
  ),
  make_option(
    c('--condition'), type='character', default=NULL,
    help='Seurat metadata column corresponding to the condition to test.'
  ),
  make_option(
    c('--subset-col'), type='character', default=NULL,
    help='Seurat metadata column corresponding to the group to subset.'
  ),
  make_option(
    c('--subset-groups'), type='character', default=NULL,
    help='Comma delimited list of groups to subset.'
  ),
  make_option(
    c('--group1'), type='character', default=NULL,
    help='Name of the first group for condition comparison.'
  ),
  make_option(
    c('--group2'), type='character', default=NULL,
    help='Name of the second group for condition comparison.'
  ),
  make_option(
    c('--test'), type='character', default='wilcox',
    help='Name of the test for Seurat FindMarkers, default is wilcox.'
  ),
  make_option(
    c('--pos'), type='logical', default=FALSE,
    help='Only test genes with a positive fold change?'
  ),
  make_option(
    c('--pct'), type='numeric', default=0,
    help='Value between 0 and 1, minimum % of cells expressing a gene to be included in the DGE test.'
  ),
  make_option(
    c('--logfc'), type='numeric', default=0,
    help='Threhsold for fold change. Does not test genes with fold change lower than this. Good for quickly finding markers.'
  ),
  make_option(
    c('--verbose'), type='logical', default=TRUE,
    help='Print info during DGE Test?'
  ),
  make_option(
    c('--slot'), type='character', default='data',
    help='Slot in Seurat object to pull expression data from.'
  ),
  make_option(
    c('--assay'), type='character', default='RNA',
    help='Assay in Seurat object to pull expression data from.'
  ),
  make_option(
    c('--latent'), type='character', default=NULL,
    help='Latent variables to account for in DGE models. Input multiple variables as a comma separated list like this: lvar1,lvar2,lvar3. Each variable should be a column in the Seurat metadata.'
  ),
  make_option(
    c('--cores'), type='numeric', default=8,
    help='How many cores?'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

# try to create the output directory if it doesn't exist:
if(!dir.exists(opt$outdir)){
  tt <- tryCatch(
    dir.create(opt$outdir),
    error=function(e) e,
    warning=function(w) w
  )
  if(is(tt,"warning")){
    stop(paste0("Permission denied, cannot create output directory ", opt$outdir))
  }
}


library(Seurat);
library(tidyverse);
library(MAST);
library(future);

# set up parallelization
plan('multicore', workers=opt$cores)

# load Seurat object:
seurat_obj <- readRDS(opt$seurat)

print('Data loaded successfully!')

# subset the seurat object based on the input:
subset_col <- opt[["subset-col"]]
subset_list <- opt[["subset-groups"]]
if(!is.null(subset_col)){

  # split into sublists:
  if(grepl(';', subset_list)){
    subset_list <- strsplit(subset_list, ';')[[1]]
    subset_col <- strsplit(subset_col, ';')[[1]]
    if(length(subset_col) != length(subset_list)){
      stop("Conflicting number of groups to subset, should be the same number of ';' characters.")
    }
  }

  for(i in 1:length(subset_list)){
    cur_subset <- strsplit(subset_list[i], ',')[[1]]
    print(paste0("subsetting by ", subset_col[i], ": ", paste0(cur_subset, collapse=',')))
    seurat_obj <- seurat_obj[,seurat_obj@meta.data[[subset_col[i]]] %in% cur_subset]
  }
}

# get a list of all groups:
cell_groups <- seurat_obj@meta.data[,opt$cluster] %>% unique %>% as.character
print(cell_groups)

# get current cluster based on the index:
cur_group <- cell_groups[opt$index]

# get a list of latent vars if present:
if(!is.null(opt$latent)){
  latent = str_split(opt$latent, ',')[[1]]
} else{
  latent = NULL
}

print(latent)

################################################################################
# is this a condition test?
################################################################################

if(opt$type == 'markers'){
  print('markers test')


  # reset idents based on this group:
  Idents(seurat_obj) <- ifelse(
    seurat_obj@meta.data[,opt$cluster] == cur_group,
    cur_group,
    'Rest'
  )

  print(table(Idents(seurat_obj)))

  # run DGE test!!!
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = cur_group,
    ident.2 = "Rest",
    slot = opt$slot,
    assay = opt$assay,
    test.use = opt$test,
    min.pct = opt$pct,
    logfc.threshold = opt$logfc,
    only.pos = opt$pos,
    latent.vars = latent
  )

  # add a column for the gene and for the cluster:
  markers$gene <- rownames(markers)
  markers$group <- cur_group

  # write to output file:
  write.csv(markers, file=paste0(opt$outdir, '/', opt$name, '_', cur_group, '.csv'), quote=FALSE, row.names=FALSE)

  print(head(markers))


} else if(opt$type == 'conditions'){
  print('conditions test')

  # subset seurat object to just this cluster!
  seurat_obj <- seurat_obj[,seurat_obj@meta.data[,opt$cluster] == cur_group]

  # reset idents:
  Idents(seurat_obj) <- seurat_obj@meta.data[,opt$condition]

    # run DGE test!!!
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = opt$group1,
      ident.2 = opt$group2,
      slot = opt$slot,
      assay = opt$assay,
      test.use = opt$test,
      min.pct = opt$pct,
      logfc.threshold = opt$logfc,
      only.pos = opt$pos,
      latent.vars = latent
    )

    # add a column for the gene and for the cluster:
    markers$gene <- rownames(markers)
    markers$group <- cur_group
    markers$ident1 <- opt$group1
    markers$ident2 <- opt$group2

    # write to output file:
    filename <- paste0(opt$outdir, '/', opt$name, '_', cur_group, '.csv')
    filename <- gsub(" ", "_", filename)
    write.csv(markers, file=filename, quote=FALSE, row.names=FALSE)

} else{
  stop('invalid option')
}
