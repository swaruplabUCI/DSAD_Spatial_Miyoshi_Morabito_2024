

FeatureEmbedding <- function(
  seurat_obj,
  features,
  raster = TRUE,
  facet_by = NULL,
  slot = 'data',
  assay = NULL,
  reduction = 'umap',
  dpi = 200,
  dpi_scale=0.5,
  ncol = 4,
  combine=TRUE,
  order_points=TRUE, # 'shuffle'
  point_size = 0.5,
  plot_max = 'q100',
  plot_min = 'q0',
  same_range = FALSE,
  colfunc = viridis::inferno,
  rev_colors = TRUE
){

  input_plot_max <- plot_max
  input_plot_min <- plot_min

  plot_list <- list()

  if(same_range){
    tmp <- GetAssayData(seurat_obj, slot=slot, assay=assay)[features,]
    plot_range <- range(tmp)
    print(plot_range)

    if(is.null(plot_max)){plot_max <- plot_range[2]}
    if(is.null(plot_min)){plot_min <- plot_range[1]}
    print(plot_max)
    print(plot_min)
  }


  for(feature in features){

    plot_df <- seurat_obj@meta.data
    plot_df$plot_x_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$plot_y_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,2]


    # check if the feature is in the meta-data
    if(feature %in% rownames(seurat_obj)){
      if(is.null(assay)){assay <- seurat_obj@active.assay}
      plot_df$plotfeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    } else if(feature %in% colnames(plot_df)){
      if(!is.numeric(plot_df[[feature]])){
        stop("Specified feature is not numeric. Try plotting with VisDimPlot?")
      }
      plot_df$plotfeature <- plot_df[[feature]]
    } else{
      stop("feature not found in rownames(seurat_obj) or in colnames(seurat_obj@meta.data).")
    }

    if(!same_range){
      plot_max <- input_plot_max
      plot_min <- input_plot_min

      plot_range <- range(plot_df$plotfeature)
      if(!is.null(plot_max)){
        if(is.character(plot_max)){
          quant <- as.numeric(gsub('q', '', plot_max)) / 100
          plot_max <- as.numeric(quantile(plot_df$plotfeature, quant))
        }
        plot_range[2] <- plot_max
        print(plot_max)
        plot_df$plotfeature <- ifelse(
          plot_df$plotfeature > plot_max,
          plot_max,
          plot_df$plotfeature
        )
      }

      if(is.character(plot_min)){
        quant <- as.numeric(gsub('q', '', plot_min)) / 100
        plot_min <- as.numeric(quantile(plot_df$plotfeature, quant))
      }
      plot_range[1] <- plot_min
    }
    #
    # shuffle points:
    if(order_points == TRUE){
      plot_df <- plot_df %>% dplyr::arrange(plotfeature)
    } else if(order_points == "shuffle"){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }


    # initialize gpgplot
    p <- plot_df %>% subset(plotfeature > plot_min) %>%
      ggplot(aes_(x=~plot_x_coord, y=~plot_y_coord,color=~plotfeature))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p +
        ggrastr::rasterise(
          geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size),
          dpi=dpi, scale=dpi_scale
        )  +
        ggrastr::rasterise(
          geom_point(size=point_size),
          dpi=dpi, scale=dpi_scale
        )
    } else{
      p <- p +
        geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size) +
        geom_point(size=point_size)
    }

    # add extras to plot:
    colors <- colfunc(256)
    if(rev_colors){colors <- rev(colors)}
    p <- p +
      labs(color = feature) +
      scale_color_gradientn(colors=colors, limits = plot_range) +
      coord_equal() +
      theme(
        plot.title = element_text(hjust=0.5, face='plain'),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
      )

    plot_list[[feature]] <- p

  }

  # if there's only one feature plotted
  if(length(plot_list) == 1){
    p <- plot_list[[1]]

    # do we want to facet?
    if(!is.null(facet_by)){
      p <- p + facet_wrap( ~ get(facet_by), ncol=ncol)
    }

    return(p)
  }

  if(combine){
    patch <- wrap_plots(plot_list, ncol=ncol)
    return(patch)
  } else{
    return(plot_list)
  }


}




custom_vln <- function(
  seurat_obj,
  features,
  group.by,
  split.by=NULL,
  groups = NULL,
  selected_split = NULL,
  plot.margin = margin(0,0,0,0, "cm"),
  slot = 'data',
  assay = NULL,
  raster_dpi=200,
  add_boxplot = TRUE,
  pt.size = 0,
  add.noise = TRUE,
  line.size=NULL,
  adjust=1,
  quantiles=c(0.5),
  stat_method = 'wilcox.test',
  comparisons = NULL,
  pval_y_adjust=0.7,
  ref_group = '.all.',
  group_color_df = NULL,
  split_colors = NULL,
  add_colorbar = TRUE,
  plot_ymin = 0,
  plot_ymax = NULL
){

  if(is.null(assay)){assay <- seurat_obj@active.assay}

  seurat_meta <- seurat_obj@meta.data

  if(class(seurat_meta[[group.by]]) != 'factor'){
    seurat_meta[[group.by]] <- as.factor(seurat_meta[[group.by]])
  }

  # get list of groups from seurat object
  if(is.null(groups)){
      groups <- levels(seurat_obj@meta.data[[group.by]])
  } else{
    if(!all(groups %in% seurat_obj@meta.data[[group.by]])){
      stop('Some groups not present in seurat_obj@meta.data[[group.by]]')
    }
  }

  if(!is.null(selected_split)){
    if(!all(selected_split %in% seurat_obj@meta.data[[split.by]])){
      stop('Some selected_split not present in seurat_obj@meta.data[[split.by]] ')
    }
  }

  # number of individual plots to make
  n_plots <- length(features) * length(groups)

  # colors for groups
  if(is.null(group_color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(seurat_meta[[group.by]])),
      group_name = levels(seurat_meta[[group.by]])
    )

    p <- seurat_meta %>%
      ggplot(aes_string(x=1, y=1, color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    group_color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    group_color_df$group <- factor_df$group_name

  } else if(is.character(group_color_df)){
    # check if it's a named character vec:
    if(is.null(names(group_color_df))){
      stop("color_df must be a dataframe with a group column and a colour column, or color_df can be a named character vector where the values are the colors and the names are the corresponding groups.")
    }
    group_color_df <- data.frame(
      colour = as.character(group_color_df),
      group = names(group_color_df)
    )
  }

  # only keep groups that we are using
  group_color_df <- subset(group_color_df, group %in% groups)
  group_color_df$var <- 1

    # color scheme
  group_colors <- group_color_df$colour
  names(group_colors) <- as.character(group_color_df$group)

  # makr color bar in ggplot
  if(add_colorbar){
    colorbar <- group_color_df %>%
      ggplot(aes(x=group, y=var, fill=group)) +
      geom_tile() +
      scale_fill_manual(values=group_color_df$colour) +
      NoLegend() +
      theme(
        plot.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.margin=margin(0,0,0,0)
      ) + RotatedAxis()
  }

  # initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_plots,style = 3, width = 50, char = "=")

  patch_list <- list()
  i <- 1
  for(feature in features){

    plot_df <- seurat_meta

    if(feature %in% colnames(plot_df)){
      plot_df$PlotFeature <- plot_df[[feature]]
    }
    else{
      plot_df$PlotFeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    }

    plot_range <- range(plot_df$PlotFeature)

    if(add.noise){
      noise <- rnorm(n = length(x = plot_df$PlotFeature)) / 100000
      plot_df$PlotFeature <- plot_df$PlotFeature + noise
    }

    # subset selected groups to compare only:
    if(!is.null(selected_split)){
      plot_df <- plot_df[plot_df[[split.by]] %in% selected_split,]
    }

    plot_list <- list()

    for(cur_group in groups){
      cur_df <- plot_df[as.character(plot_df[[group.by]]) == cur_group,]

      if(is.null(split.by)){
        p <- cur_df %>% ggplot(aes_string(x=group.by, y='PlotFeature', fill=group.by)) +
        scale_fill_manual(values=group_colors, na.value = 'grey90')
      } else{
        p <- cur_df %>% ggplot(aes_string(x=split.by, y='PlotFeature', fill=split.by))
        if(!is.null(split_colors)){
          p <- p + scale_fill_manual(values=split_colors)
        }
      }

      # add violin
      p <- p + geom_violin(
          trim=FALSE, adjust=adjust,
          scale='width', draw_quantiles = quantiles,
          color='black', lwd=0.5
        )

      if(add_boxplot){
        p <- p + geom_boxplot(width=0.1, fill='white', outlier.shape=NA)
      }

      if(pt.size > 0){
        p <- p +  ggrastr::rasterise(ggbeeswarm::geom_beeswarm(size=0.5), dpi=raster_dpi)
      }

      # code for stacked vln plot theme from CellChat:
      p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
        theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
      p <- p + theme(plot.title= element_blank(), # legend.position = "none",
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_text(size = rel(1), angle = 0),
                     axis.text.y = element_text(size = rel(1)),
                     plot.margin = plot.margin ) +
        theme(axis.text.y = element_text(size = 8))
      p <- p + theme(element_line(size=line.size))

      if(length(plot_list) > 0){
        p <- p + theme(
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
        )
      }


      # y limits:
      if(!is.null(plot_ymax)){
        ymax <- plot_ymax
      } else{
        ymax <- ceiling(plot_range[2]); 
      }
      ymin <- floor(plot_range[1])
      p <- p + expand_limits(y=c(ymin, ymax)) +
        scale_y_continuous(expand=c(0,0), breaks = c(ymax), limits=c(plot_ymin, NA))

      # stats:
      #if(is.na(comparisons)){
        if(!is.null(split.by)){
          if(is.null(comparisons)){
            
            p <- p + stat_compare_means(method=stat_method, aes(label = ..p.signif..), label.y=ymax*pval_y_adjust,hide.ns=TRUE)
          } else{
            print(comparisons)
            p <- p + stat_compare_means(
              method=stat_method,
              label='p.signif',
              comparisons=comparisons,
              label.y=ymax-(ymax*0.25),
              hide.ns=TRUE
             # ref.group=ref_group
            )
          }
        }
      #}


      p <- p + ylab(feature)

      # different settings for the last gene:
      if(!add_colorbar & feature == features[length(features)]){
        p <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
          RotatedAxis()
      }

      plot_list[[cur_group]] <- p

      # update progress bar
      setTxtProgressBar(pb, i)
      i <- i+1

    }
    patch_list[[feature]] <- wrap_plots(plot_list, ncol=length(groups))

  }

  # close progress bar
  close(pb)

  if(length(patch_list) == 1){
    return(patch_list[[1]])
  }

  plot_heights <- rep(50 / length(features), length(features))

  if(add_colorbar){
    plot_heights <- c(plot_heights, 1)
    patch_list <- c(patch_list, list(colorbar))
  }

  wrap_plots(patch_list, ncol=1) + plot_layout(heights=plot_heights, guides='collect')
}




PlotDEGsChromosome <- function(
  degs,
  gene_table, # matches the gene names in teh seurat object with ensembl ids
  EnsDb,
  chr_names,
  gene_biotypes = c("protein_coding"),
  plot_limit = 0.5,
  logfc_thresh = 0.25,
  raster = TRUE,
  raster_dpi = 200,
  color_df = NULL,
  return_table = FALSE,
  high_color = 'blue',
  mid_color = 'lightgrey',
  low_color = 'yellow'
){

  degs$avg_log2FC <- ifelse(
    abs(degs$avg_log2FC) > plot_limit,
    sign(degs$avg_log2FC) * plot_limit,
    degs$avg_log2FC
  )

  plot_list <- list()
  chr_n_genes <- rep(0, length(chr_names))
  names(chr_n_genes) <- chr_names

  # genes that are differentially expressed in at least one cell type:
  print(logfc_thresh)
  genes.keep <- unique(degs[abs(degs$avg_log2FC) >= logfc_thresh & degs$p_val_adj < 0.05, 'gene'])
  degs <- subset(degs, gene %in% genes.keep)

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){
    p <- degs %>%
      ggplot(aes(x=avg_log2FC, y=p_val_adj, color=group)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
  }

  # make the colorbar as its own heatmap
  color_df <- subset(color_df, group %in% unique(degs$group))
  color_df$var <- 1
  colorbar <- color_df %>%
    ggplot(aes(y=group, x=var, fill=group)) +
    geom_tile() +
    scale_fill_manual(values=color_df$colour) +
    NoLegend() +
    theme(
      plot.title=element_blank(),
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      plot.margin=margin(0,0,0,0)
    )
  
  # output dataframe:
  out_degs <- data.frame()

  for(cur_chr in chr_names){
    print(cur_chr)

    # subset the ensembl db based on
    gene.coords <- ensembldb::genes(EnsDb, filter = ~ gene_biotype == "protein_coding") %>% subset(seqnames == cur_chr)

    # subset by genes in our snRNA data
    gene.coords <- gene.coords[gene.coords$gene_id %in% gene_table$gene_id]
    ix <- match(gene.coords$gene_id, gene_table$gene_id)
    gene.coords$gene_name <- gene_table$gene_name[ix]

    # remove entries with duplicate names
    dupes <- names(which(table(gene.coords$gene_name) > 1))
    dupes <- unlist(lapply(dupes, function(x){which(gene.coords$gene_name == x)[-1]}))
    if(!is.null(dupes)){
      gene.coords <- gene.coords[-dupes]
    }

    # subset the DEGs:
    cur_degs <- degs %>% subset(gene %in% gene.coords$gene_name)
    gene.coords <- gene.coords %>% subset(gene_name %in% cur_degs$gene)

    # add chromosome to the table
    cur_degs$chr <- cur_chr

    if(nrow(cur_degs) == 0){
      warning('blank plot')
      plot_list[[paste0(cur_g, '_', cur_chr)]] <- ggplot() + theme_void()
      next
    }

    # add order based on genomic coordinates:
    gene.coords$order <- 1:length(gene.coords)
    chr_n_genes[cur_chr] <- chr_n_genes[cur_chr] + length(gene.coords)

    # add the order to the deg table:
    cur_degs$order <- match(cur_degs$gene, gene.coords$gene_name)

    cur_degs$avg_log2FC <- ifelse(
      cur_degs$p_val_adj > 0.05, NA,
      cur_degs$avg_log2FC
    )

    print(dim(cur_degs))
    print(length(unique(cur_degs$gene)))

    # plot as a heatmap:
    p <- cur_degs %>%
      ggplot(aes(x = order, y = group, fill = avg_log2FC))

    if(raster){
      p <- p +  ggrastr::rasterise(geom_tile(), dpi=raster_dpi)
    } else {
      p <- p + geom_tile()
    }

    p <- p +
      scale_fill_steps2(high=high_color, mid=mid_color, low=low_color, limits = c(-plot_limit, plot_limit), n.breaks=5, na.value='lightgrey') +
      xlab(cur_chr) +
      theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0)
      )

    if(length(plot_list) > -1){
      p <- p + theme(
        axis.text.y = element_blank(),
      )
    }

    plot_list[[paste0(cur_chr)]] <- p

    out_degs <- rbind(out_degs, cur_degs)

  }

  if(return_table){
    return(out_degs)
  }

  # chr_n_genes <- chr_n_genes / length(cluster_groups)
  plot_widths <- 50 * (chr_n_genes / sum(chr_n_genes))
  plot_widths <- c(1, plot_widths)
  plot_list <- c(list(colorbar), plot_list)
  print(length(plot_list))
  patch <- wrap_plots(plot_list, ncol=length(chr_names)+1) + plot_layout(widths = plot_widths, guides='collect')
  patch

}



PlotEmbedding <- function(
  seurat_obj,
  group.by,
  reduction = 'umap',
  split.by = NULL,
  plot_under = FALSE,
  label = TRUE,
  point_size = 1,
  legend_point_size = 5,
  text_size = 3,
  raster = TRUE,
  order_points = "shuffle",
  x = NULL,
  y = NULL,
  raster_dpi = 400,
  raster_scale = 1,
  selected = NULL,
  plot_theme = NULL,
  color_df = NULL,
  plot_ratio=TRUE
){

  plot_df <- seurat_obj@meta.data
  if(!is.null(x) & !is.null(y)){
    plot_df$x <- plot_df[[x]]
    plot_df$y <- plot_df[[y]]
  } else if(reduction %in% names(seurat_obj@reductions)){
    plot_df$x <- seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$y <- seurat_obj@reductions[[reduction]]@cell.embeddings[,2]
  }


  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  if(label){
    for(cur_cluster in unique(plot_df[[group.by]])){
      cur_meta <- plot_df[plot_df[[group.by]] == cur_cluster,]
      df <- data.frame(
        cluster = cur_cluster,
        x = mean(cur_meta$x),
        y = mean(cur_meta$y)
      )
      centroid_df <- rbind(centroid_df, df)
    }
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x='x', y='y', color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name

  }

  # only show selected groups
  if(!is.null(selected)){
    plot_df[[group.by]][!(plot_df[[group.by]] %in% selected)] <- NA

    if(is.factor(plot_df[[group.by]])){
      plot_df[[group.by]] <- droplevels(plot_df[[group.by]])
    }

    if(label){
      centroid_df <- subset(centroid_df, cluster %in% selected)
    }
  }

  print(levels(plot_df[[group.by]]))

  # shuffle points:
  if(order_points == "shuffle"){
    plot_df <- plot_df[sample(nrow(plot_df)),]
  }

  # plot a single embedding
  if(is.null(split.by)){
    p <- .PlotSingleEmbedding(plot_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, plot_theme=plot_theme, plot_ratio=plot_ratio)
  } else{

    split_groups <- unique(plot_df[[split.by]])
    plot_list <- lapply(split_groups, function(cur_split){
      print(cur_split)
      cur_df <- plot_df[plot_df[[split.by]] == cur_split,]
      split_df <- plot_df[plot_df[[split.by]] != cur_split,]
      .PlotSingleEmbedding(cur_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, split_df, cur_split, plot_theme, plot_under, plot_ratio=plot_ratio)
    })

    return(plot_list)

  }

  p

}

.PlotSingleEmbedding <- function(
  plot_df,
  group.by,
  label,
  raster,
  raster_dpi,
  raster_scale,
  point_size,
  legend_point_size,
  text_size,
  color_df,
  centroid_df,
  split_df = NULL,
  cur_split = NULL,
  plot_theme = NULL,
  plot_under=FALSE,
  plot_ratio = TRUE
){

  p <- plot_df %>%
    ggplot(aes_string(x='x', y='y', color=group.by))

  # add the under layer for the split plots
  if(!is.null(split_df) & plot_under){
    # add points
    if(raster){
      p <-  p + ggrastr::rasterise(geom_point(data = split_df, size=point_size/2, color='lightgrey'), dpi=raster_dpi, scale=raster_scale)
    } else{
      p <- p + geom_point(data = split_df, size=point_size/2, color='lightgrey')
    }
  }

  # add points
  if(raster){
    p <-  p + ggrastr::rasterise(geom_point(size=point_size), dpi=raster_dpi,  scale=raster_scale)
  } else{
    p <- p + geom_point(size=point_size)
  }

  # add labels
  if(label){
    p <- p + ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=text_size)
  }

  # color scheme
  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)
  p <- p + scale_color_manual(values=plot_colors, na.value = 'grey90')

  # add title:
  if(is.null(cur_split)){
    p <- p + ggtitle(group.by) + theme(plot.title=element_text(hjust=0.5))
  } else{
    p <- p + ggtitle(cur_split) + theme(plot.title=element_text(hjust=0.5))
  }

  # add theme:
  if(!is.null(plot_theme)){
    p <- p + plot_theme
  }

  # fixed coords
  if(plot_ratio){
    p <- p + coord_equal()
  }

  # adjust legend point size
  p <- p + guides( color = guide_legend(override.aes = list(size=legend_point_size)))

  p

}

VisDimPlot <- function(
  seurat_obj,
  group.by,
  sample_col = "Sample",
  sample_labels = NULL,
  samples_to_plot = NULL,
  raster = TRUE,
  dpi = 800,
  text_size = 8,
  ncol=4,
  color_df = NULL,
  combine=FALSE
){

  # check that row and col are in the seurat metadata:
  if(!all(c("col", "row") %in% colnames(seurat_obj@meta.data))){
    stop("Spatial coordinates must be present in seurat@meta.data, named row and col.")
  }

  plot_df <- seurat_obj@meta.data

  # check if the group.by is in the meta-data
  if(group.by %in% colnames(plot_df)){
    if(is.numeric(plot_df[[group.by]])){
      stop("Specified group.by a numeric variable. Try plotting with SampleFeaturePlot?")
    }
  } else{
    stop("group.by not found in colnames(seurat_obj@meta.data).")
  }


  # get list of samples to plot:
  if(is.null(samples_to_plot)){
    samples_to_plot <- unique(plot_df[[sample_col]])
  } else{
    if(!all(samples_to_plot %in% unique(plot_df[[sample_col]]))){
      stop("Some specified samples are not found in seurat_obj@meta.data[[sample_col]].")
    }
  }

  # check sample labels:
  if(!is.null(sample_labels)){
    if(!all(sample_labels %in% colnames(seurat_obj@meta.data))){
      stop("Some sample_labels are not found in colnames(seurat_obj@meta.data)")
    }
  }

  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]


  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]

  # group by sample_labels
  if(!is.null(sample_labels)){
    plot_df <- plot_df %>% group_by(across(all_of(sample_labels))) %>%
      as.data.frame()
  }

  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){
    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x=1, y=1, color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name
  } else if(is.character(color_df)){
    # check if it's a named character vec:
    if(is.null(names(color_df))){
      stop("color_df must be a dataframe with a group column and a colour column, or color_df can be a named character vector where the values are the colors and the names are the corresponding groups.")
    }
    color_df <- data.frame(
      colour = as.character(color_df),
      group = names(color_df)
    )
  }

  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)


  # get the ordering of samples to plot:
  order_df <- plot_df %>%
    group_by(across(all_of(sample_labels)))  %>%
    arrange(.by_group=TRUE) %>%
    dplyr::select(all_of(c(sample_labels, sample_col))) %>%
    distinct()

  samples_to_plot <- as.character(order_df[[sample_col]])


  plot_list <- list()
  for(cur_sample in samples_to_plot){
    print(cur_sample)

    # get data for this sample and get hex vertex coordinates
    cur_df <- plot_df[plot_df[[sample_col]] == cur_sample,]
    vertices <- make_hex_spots(cur_df, cur_df[[group.by]])

    # get sample labels:
    if(!is.null(sample_labels)){
      cur_labels <- unlist(lapply(sample_labels, function(x){unique(as.character(cur_df[[x]]))}))
      plot_label <- paste0(cur_labels, collapse=', ')
    } else{
      plot_label <- cur_sample
    }

    # initialize gpgplot
    p <- vertices  %>%
      ggplot(aes_(x=~x.vertex, y=~y.vertex,group=~spot, fill=~fill))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p + ggrastr::rasterise(geom_polygon(size=0), dpi=dpi)
    } else{
      p <- p + geom_polygon(size=0)
    }

    # add extras to plot:
    p <- p +
      labs(fill = group.by) +
      scale_fill_manual(values=plot_colors, na.value = 'grey90') +
      coord_equal() +
      ggtitle(plot_label) +
      theme(
        plot.title = element_text(hjust=0.5, face='plain', size=text_size),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()

      )

    plot_list[[cur_sample]] <- p

  }

  if(length(plot_list) == 1){
    return(plot_list[[1]])
  } else if(!combine){
    return(plot_list)
  } else{
    patch <- wrap_plots(plot_list, ncol=ncol, widths=1, heights=1) +
      plot_layout(guides='collect') +
      plot_annotation(title = group.by) &  theme(plot.title = element_text(hjust=0.5),)
  }
  patch

}

SampleFeaturePlot <- function(
  seurat_obj,
  feature,
  sample_col = "Sample",
  sample_labels = NULL,
  samples_to_plot = NULL,
  raster = TRUE,
  slot = 'data',
  assay = NULL,
  dpi = 800,
  ncol = 4,
  plot_max = 'q100',
  plot_min = 'q0',
  text_size = 8,
  colfunc = colorRampPalette(brewer.pal(11, 'Spectral' )),
  rev_colors = FALSE,
  combine=TRUE
){

  # check that row and col are in the seurat metadata:
  if(!all(c("col", "row") %in% colnames(seurat_obj@meta.data))){
    stop("Spatial coordinates must be present in seurat@meta.data, named row and col.")
  }


  plot_df <- seurat_obj@meta.data

  # check if the feature is in the meta-data
  if(feature %in% rownames(seurat_obj)){
    if(is.null(assay)){assay <- seurat_obj@active.assay}
    plot_df$SpatialPlotFeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
  } else if(feature %in% colnames(plot_df)){
    if(!is.numeric(plot_df[[feature]])){
      stop("Specified feature is not numeric. Try plotting with VisDimPlot?")
    }
    plot_df$SpatialPlotFeature <- plot_df[[feature]]
  } else{
    stop("feature not found in rownames(seurat_obj) or in colnames(seurat_obj@meta.data).")
  }

  # get list of samples to plot:
  if(is.null(samples_to_plot)){
    samples_to_plot <- unique(plot_df[[sample_col]])
  } else{
    if(!all(samples_to_plot %in% unique(plot_df[[sample_col]]))){
      stop("Some specified samples are not found in seurat_obj@meta.data[[sample_col]].")
    }
  }

  # check sample labels:
  if(!is.null(sample_labels)){
    if(!all(sample_labels %in% colnames(seurat_obj@meta.data))){
      stop("Some sample_labels are not found in colnames(seurat_obj@meta.data)")
    }
  }

  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]

  # group by sample_labels
  if(!is.null(sample_labels)){
    plot_df <- plot_df %>% group_by(across(all_of(sample_labels))) %>%
      as.data.frame()
  }

  plot_range <- range(plot_df$SpatialPlotFeature)
  if(!is.null(plot_max)){
    if(is.character(plot_max)){
      quant <- as.numeric(gsub('q', '', plot_max)) / 100
      plot_max <- as.numeric(quantile(plot_df$SpatialPlotFeature, quant))
    }
    plot_range[2] <- plot_max
    print(plot_max)
    plot_df$SpatialPlotFeature <- ifelse(
      plot_df$SpatialPlotFeature > plot_max,
      plot_max,
      plot_df$SpatialPlotFeature
    )
  }

  if(is.character(plot_min)){
    quant <- as.numeric(gsub('q', '', plot_min)) / 100
    plot_min <- as.numeric(quantile(plot_df$SpatialPlotFeature, quant))
  }
  plot_range[1] <- plot_min

  # get the ordering of samples to plot:
  order_df <- plot_df %>%
    group_by(across(all_of(sample_labels)))  %>%
    arrange(.by_group=TRUE) %>%
    dplyr::select(all_of(c(sample_labels, sample_col))) %>%
    distinct()

  samples_to_plot <- as.character(order_df[[sample_col]])

  plot_list <- list()
  for(cur_sample in samples_to_plot){
    print(cur_sample)

    # get data for this sample and get hex vertex coordinates
    cur_df <- plot_df[plot_df[[sample_col]] == cur_sample,]
    vertices <- make_hex_spots(cur_df, cur_df$SpatialPlotFeature)

    # get sample labels:
    if(!is.null(sample_labels)){
      cur_labels <- unlist(lapply(sample_labels, function(x){unique(as.character(cur_df[[x]]))}))
      plot_label <- paste0(cur_labels, collapse=', ')
    } else{
      plot_label <- cur_sample
    }

    # initialize gpgplot
    p <- vertices %>% subset(fill > plot_min) %>%
      ggplot(aes_(x=~x.vertex, y=~y.vertex,group=~spot, fill=~fill))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p +
        ggrastr::rasterise(
          geom_polygon(data = subset(vertices, fill <= plot_min), fill='lightgrey', size=0),
        dpi=dpi/2
      ) +
      ggrastr::rasterise(geom_polygon(size=0), dpi=dpi)
    } else{
      p <- p +
        geom_polygon(data = subset(vertices, fill <= plot_min), fill='lightgrey', size=0) +
        geom_polygon(size=0)
    }

    # add extras to plot:
    #colfunc <- colorRampPalette(brewer.pal(11, 'Spectral' ))
    colors <- colfunc(256)
    if(rev_colors){colors <- rev(colors)}
    p <- p +
      labs(fill = feature) +
      scale_fill_gradientn(colors=colors, limits = plot_range) +
      coord_equal() +
      ggtitle(plot_label) +
      theme(
        plot.title = element_text(hjust=0.5, face='plain', size=text_size),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()

      )

    plot_list[[cur_sample]] <- p

  }


  if(length(plot_list) == 1){
    return(plot_list[[1]])
  } else if(!combine){
    return(plot_list)
  }
    else {
    patch <- wrap_plots(plot_list, ncol=ncol, widths=1, heights=1) +
      plot_layout(guides='collect') +
      plot_annotation(title = feature) &  theme(plot.title = element_text(hjust=0.5),)
  }
  patch
}











################################################################################
# Helper functions
################################################################################



make_hex_spots <- function(cdata, fill){
  r <- 1/2
   R <- (2/sqrt(3)) * r

  spot_positions <- select_spot_positions(cdata, fill = fill)
  spot_positions <- adjust_hex_centers(spot_positions)
  vertex_offsets <- data.frame(x.offset = c(0, r, r, 0, -r,
      -r), y.offset = c(-R, -R/2, R/2, R, R/2, -R/2))
  spot_vertices <- make_spot_vertices(spot_positions, vertex_offsets)
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  spot_vertices
}

select_spot_positions <- function (cdata, x = "col", y = "row", fill = "spatial.cluster")
{
    assertthat::assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
    }
    else if (is.vector(fill) || is.factor(fill)) {

        assertthat::assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y)]
        colnames(spot_positions) <- c("x.pos", "y.pos")
        spot_positions$fill <- fill
    }
    spot_positions$spot <- rownames(spot_positions)
    spot_positions
}

adjust_hex_centers <- function (spot_positions)
{
    r <- 1/2
    R <- (2/sqrt(3)) * r

    spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) +
        1
    spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) +
        1
    spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
    spot_positions$x.pos <- (spot_positions$x.pos + 1)/2

    spot_positions
}

make_spot_vertices <- function (spot_positions, vertex_offsets)
{
    spot_vertices <- merge(spot_positions, vertex_offsets)
    spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
    spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
    as.data.frame(spot_vertices)
}




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



netVisual_embeddingPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, color.use = NULL, point.shape = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                        xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    # pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
    #  pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  # color dots (light inside color and dark border) based on clustering and no labels
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = clusters, colour = clusters, shape = group)) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) #+ scale_alpha(group, range = c(0.1, 1))
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5, max.overlaps=Inf) + scale_alpha_discrete(range = c(1, 0.6))
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}
