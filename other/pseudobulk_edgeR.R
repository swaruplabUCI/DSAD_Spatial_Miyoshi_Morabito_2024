

pseudobulk_edgeR <- function(
  seurat_obj,
  cell_type_col = 'cell_type',
  label_col = 'Diagnosis',
  replicate_col = 'Sample',
  covariates = NULL,
  slot = 'counts',
  assay = 'RNA',
  cells_use = NULL
){

  de_family <- 'pseudobulk'
  de_method <- 'edgeR'
  de_type <- 'LRT'

  # get expression matrix from seurat object
  X <- GetAssayData(seurat_obj, slot=slot, assay=assay)
  meta <- seurat_obj@meta.data

  # are we using a subset of the cells?
  if(!is.null(cells_use)){
    X <- X[,cells_use]
    meta <- meta[cells_use,]
  }

  print(dim(X))

  # set up sample level metadata:
  sample_vars <- c(replicate_col, label_col, covariates)
  sample_meta <- meta[,sample_vars] %>% distinct()
  print(sample_vars)

  # make group sample column
  sample_meta$group_sample <- paste0(
    as.character(sample_meta[[replicate_col]]), ':',
    as.character(sample_meta[[label_col]])
  )
  rownames(sample_meta) <- 1:nrow(sample_meta)

  # use libra to make pseudobulk replicates:
  matrices <- Libra::to_pseudobulk(
    X, meta = meta,
    cell_type_col = cell_type_col,
    replicate_col = replicate_col,
    label_col = label_col
  )

  # initialize progress bar:

  # loop over all cell groups
  cell_groups <- names(matrices)
  results <- data.frame()
  pb <- utils::txtProgressBar(min = 0, max = length(cell_groups), style = 3, width = 50, char = "=")
  counter <- 1
  for(cur_group in cell_groups){

    # update progress bar:
    setTxtProgressBar(pb, counter)

    # get matrix for this celltype
    x <- matrices[[cur_group]]

    # set up sample metadata
    targets = data.frame(group_sample = colnames(x)) %>%
          mutate(group = gsub(".*\\:", "", group_sample))

    # merge with sample meta:
    gs <- targets$group_sample
    targets <- merge(sample_meta, targets, by='group_sample')
    ix <- match(gs,targets$group_sample)
    targets <- targets[ix,]

    # keep factor level from seurat obj
    if (is.factor(meta[[label_col]])) {
      targets$group %<>% factor(levels = levels(meta[[label_col]]))
    }

    # set up formula with covariates:
    if(!is.null(covariates)){
      form <- as.formula(paste(c("~", paste(c(label_col, covariates), collapse= ' + ')), collapse=' '))
    } else{
      form <- as.formula(paste(c("~", label_col), collapse=' '))
    }

    # set up design matrix:
    design = model.matrix(form, data = targets)

    # run edgeR
    y <- edgeR::DGEList(counts = x, group = targets[[label_col]]) %>%
      edgeR::calcNormFactors(method = 'TMM') %>%
      edgeR::estimateDisp(design)

    # fit glm
    fit <- edgeR::glmFit(y, design = design)
    test <- edgeR::glmLRT(fit)

    # set up results
    res = topTags(test, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column('gene') %>%
      mutate(de_family = 'pseudobulk',
             de_method = de_method,
             de_type = de_type,
             cell_type = cur_group)

     colnames(res) %<>%
       fct_recode('p_val' = 'p.value',  ## DESeq2
                  'p_val' = 'pvalue',  ## DESeq2
                  'p_val' = 'p.value',  ## t/wilcox
                  'p_val' = 'P.Value',  ## limma
                  'p_val' = 'PValue'  , ## edgeR
                  'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                  'p_val_adj' = 'adj.P.Val',      ## limma
                  'p_val_adj' = 'FDR',            ## edgeER
                  'avg_logFC' = 'log2FoldChange', ## DESEeq2
                  'avg_logFC' = 'logFC', ## limma/edgeR
                  'avg_logFC' = 'avg_log2FC' # Seurat V4
       ) %>%
       as.character()


    # remove unnecessary cols and reorder:
    res %<>% dplyr::select(c(cell_type, gene, avg_logFC, p_val, p_val_adj, de_family, de_method, de_type))

    # add to ongoing results:
    results <- rbind(results, res)

    counter <- counter + 1

  }

  # close progress bar
  close(pb)

  results

}
