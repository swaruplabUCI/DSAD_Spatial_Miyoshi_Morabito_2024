
# script to run DEG plotting in parallel
library(optparse)

option_list = list(
  make_option(
    c('--infile'), type='character', default=NULL,
    help='DEG results as a delimited file (.csv by default)', metavar='character'
  ),
  make_option(
    c('--outdir'), type='character', default='./',
    help='Directory to place output files'
  ),
  make_option(
    c('--sep'), type='character', default=',',
    help='Delimiter for DEG results file, comma delimited by default'
  ),
  make_option(
    c('--name'), type='character', default=NULL,
    help='name for output files'
  ),
  make_option(
    c('--seurat'), type='character', default=NULL,
    help='Cell Dataset .rds file'
  ),
  make_option(
    c('--cluster'), type='character', default='seurat_clusters',
    help='Column in seurat metadata that indicates clusters'
  ),
  make_option(
    c('--volcano'), type='logical', default=TRUE,
    help='logical (TRUE or FALSE) indicating whether or not to produce volcano plots'
  ),
  make_option(
    c('--nlabel-volcano'), type='numeric', default=5,
    help='number of genes (up & down-reg) to label on the volcano plot'
  ),
  make_option(
    c('--volcano-color1'), type='character', default="darkgoldenrod3",
    help='Color of the first group on the volcano plot (right side of plot)'
  ),
  make_option(
    c('--volcano-color2'), type='character', default="hotpink3",
    help='Color of the second group on the volcano plot (left side of plot)'
  ),
  make_option(
    c('--volcano-w'), type='character', default=5,
    help='Width of volcano plot'
  ),
  make_option(
    c('--volcano-h'), type='numeric', default=5,
    help='Height of volcano plot'
  ),
  make_option(
    c('--enrichr'), type='logical', default=TRUE,
    help='logical (TRUE or FALSE) indicating whether or not to run enrichR'
  ),
  make_option(
    c('--raster'), type='logical', default=TRUE,
    help='logical (TRUE or FALSE) indicating whether or not to rasterise data points in the volcano plots.'
  ),
  make_option(
    c('--dpi'), type='numeric', default=800,
    help='pixels per inch for rasterised plots.'
  ),
  make_option(
    c('--ngenes-enrichr'), type='numeric', default=200,
    help='How many DEGs by fold-change to include for enrichment testing?'
  ),
  make_option(
    c('--nterms-enrichr'), type='numeric', default=15,
    help='How many enriched terms to plot?'
  ),
  make_option(
    c('--dbs'), type='character', default="GO_Biological_Process_2021,GO_Cellular_Component_2021,GO_Molecular_Function_2021,WikiPathways_2019_Mouse,KEGG_2019_Mouse",
    help='Comma delimited list of enrichR databases.'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

################################################################################
# setup variables based on optparse
################################################################################

# volcano plot colors
color1 = opt["volcano-color1"]
color2 = opt["volcano-color2"]
nlabel = as.numeric(opt[["nlabel-volcano"]])
volcano_width = as.numeric(opt[["volcano-w"]])
volcano_height = as.numeric(opt[["volcano-h"]])
make_volcano = as.logical(opt[["volcano"]])
name = opt$name
clusters = opt$cluster

# enrichR dbs
dbs = unlist(strsplit(opt$dbs, ','))
ngenes_enrichr = as.numeric(opt[["ngenes-enrichr"]])
nterms = as.numeric(opt[["nterms-enrichr"]])

# create output dir:
dir.create(opt$outdir)

################################################################################
# load packages
################################################################################

library(Seurat);
library(tidyverse);
library(cowplot);
library(viridis);
library(ggpubr);
library(ggrepel);
library(patchwork);
library(RColorBrewer);
library(enrichR);
theme_set(theme_cowplot())

################################################################################
# load DEG file & Seurat object
################################################################################

# load DEGs
cur_degs <- read.table(
  opt$infile, sep = opt$sep, header=TRUE
)
print('DEGs loaded successfully!')

# load Seurat object:
seurat_obj <- readRDS(opt$seurat)
print('Seurat data loaded successfully!')


################################################################################
# get cluster colors from seurat obj
################################################################################

cur_degs$gene <- as.character(cur_degs$gene)

if('ident1' %in% colnames(cur_degs)){
  group1 <- as.character(unique(cur_degs$ident1)); group2 <- as.character(unique(cur_degs$ident2))
} else{
  group1 <- 'Upregulated'; group2 <- 'Downregulated'
}

# get color scheme:
p <- DimPlot(seurat_obj, group.by=clusters, reduction='umap', label=T)
g <- ggplot_build(p)
colors <- g$data[[1]]["colour"]
groups <- g$data[[1]]['group']
color_df <- unique(data.frame(colors, groups)) %>% arrange(group)
color_df$group <- unique(seurat_obj@meta.data[,clusters])
color_scheme <- color_df$colour
names(color_scheme) <- color_df$group

################################################################################
# setup DEG table for volcano plots
################################################################################

# set groups asa a factor
cur_degs$group <- factor(as.character(cur_degs$group), levels=unique(cur_degs$group)[order(as.character(unique(cur_degs$group)))])

# label the top and bottom significant genes by log fold change
cur_degs <- Reduce(rbind, lapply(unique(cur_degs$group), function(x){
  cur <- subset(cur_degs, group == x)

  top_thresh <- cur %>% subset(p_val_adj <= 0.05) %>% top_n(nlabel, wt=avg_log2FC) %>% .$avg_log2FC %>% min
  bottom_thresh <- cur %>% subset(p_val_adj <= 0.05) %>% top_n(-1*nlabel, wt=avg_log2FC) %>% .$avg_log2FC %>% max

  cur$anno <- ifelse(cur$p_val_adj <= 0.05 & cur$avg_log2FC >= top_thresh, cur$gene, NA)
  cur$anno <- ifelse(cur$p_val_adj <= 0.05 & cur$avg_log2FC <= bottom_thresh, cur$gene, cur$anno)
  cur$color <- ifelse(cur$p_val_adj > 0.05, 'gray', ifelse(cur$avg_log2FC > 0, color1, color2))
  cur
}))

################################################################################
# build volcano plot
# 1: plot significance line
# 2: plot all points
# 3: plot points with labels (black outline)
# 4: plot the text labels
# 5+: formatting
################################################################################

# are we making a volcano?
if(make_volcano){

  p_val_min <- sort(unique(cur_degs$p_val_adj))[2]
  cur_degs$p_val_adj <- ifelse(cur_degs$p_val_adj == 0, p_val_min, cur_degs$p_val_adj)

  # individual volcano plots:
  plot_list <- list()
  for(cluster in as.character(unique(cur_degs$group))){
    print(cluster)
    plot_degs <-  cur_degs %>% subset(group == cluster)
    p <- plot_degs %>%
       ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj))) +
       geom_hline(yintercept=-log10(0.05), linetype='dashed')

    if(opt$raster){
      p <- p + ggrastr::rasterise(geom_point(alpha=0.5, color=plot_degs$color), dpi=opt$dpi)
    } else{
      p <- p + geom_point(alpha=0.5, color=plot_degs$color)
    }

    p <- p +
       geom_point(inherit.aes=FALSE, data=subset(plot_degs, !is.na(anno)), aes(avg_log2FC, -log10(p_val_adj)),fill=subset(plot_degs, !is.na(anno)) %>% .$color, shape=21, size=2, color='black') +
       geom_text_repel(aes(label=anno), color='black', fontface='italic',  min.segment.length=0) +
       xlim(-1*max(abs(plot_degs$avg_log2FC))-0.1, max(abs(plot_degs$avg_log2FC))+0.1) +
       ggtitle(paste0(cluster)) +
       theme(
         panel.grid.major = element_blank(),
         plot.title = element_text(hjust = 0.5),
         legend.position='bottom'
       )

      plot_list[[cluster]] <- p
  }

  pdf(paste0(opt$outdir, '/volcano_', opt$name, '.pdf'), width=volcano_width, height=volcano_height, useDingbats=FALSE)
  for(p in plot_list){print(p)}
  dev.off()

} else{
  print('Skipping Volcano')
}


################################################################################
# run enrichR
################################################################################

if(opt$enrichr){

  dir.create(paste0(opt$outdir, '/enrichr/'))

  combined_output <- data.frame()
  for(cur_group in unique(cur_degs$group)){

    print(cur_group)
    genes.up <- cur_degs %>% subset(group == cur_group & p_val_adj <= 0.05 & avg_log2FC >= 0) %>% top_n(ngenes_enrichr, wt=avg_log2FC) %>% .$gene %>% as.character
    genes.down <- cur_degs %>% subset(group == cur_group & p_val_adj <= 0.05 & avg_log2FC < 0) %>% top_n(-ngenes_enrichr, wt=avg_log2FC) %>% .$gene %>% as.character

    if(length(genes.up) < 10 | length(genes.down) < 5){next}

    enriched_up <- enrichr(genes.up,dbs)
    enriched_down <- enrichr(genes.down,dbs)

    for(db in dbs){
      cur_df_up <- enriched_up[[db]]
      cur_df_down <- enriched_down[[db]]

      if (nrow(cur_df_up) > 1){
        cur_df_up$db <- db
        cur_df_up$group <- cur_group
        cur_df_up$upregulated <- TRUE
        combined_output <- rbind(combined_output, cur_df_up)
      }

      if (nrow(cur_df_down) > 1){
        cur_df_down$db <- db
        cur_df_down$group <- cur_group
        cur_df_down$upregulated <- FALSE
        combined_output <- rbind(combined_output, cur_df_down)
      }

    }
  }

  write.table(combined_output, file=paste0(opt$outdir, '/', name, '_GO_terms.tsv'), quote=FALSE, row.names=FALSE, sep='\t')


  # plot results
  wrapText <- function(x, len) {
      sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
  }

  for(cur_group in unique(cur_degs$group)){

    print(cur_group)

    cur_terms <- subset(combined_output, group == cur_group)
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)

    cur_color <- color_df %>% subset(group == cur_group) %>% .$colour

    # plot top terms as barplot
    plot_list <- list()
    for(cur_db in dbs){

      print(cur_db)
      plot_df <- cur_terms %>% subset(db == cur_db & upregulated == TRUE) %>% top_n(nterms, wt=Combined.Score)
      plot_df$Term <- do.call(rbind, strsplit(plot_df$Term, '[(]GO:'))[,1]
      plot_df$wrap <- wrapText(plot_df$Term, 50)

      p1 <- ggplot(plot_df, aes(x=Combined.Score, y=reorder(wrap, Combined.Score)))+
      geom_bar(stat='identity', position='identity', color='white', fill=cur_color, alpha=0.75) +
      geom_text(aes(label=wrap), x=1, color="black", size=3.5, hjust='left') +
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) + ylab('Term') + xlab('Combined.Score') + ggtitle(group1)

      if(nrow(plot_df) == 0){
        p1 <- ggplot() + theme_void()
      }

      print(cur_db)
      plot_df <- cur_terms %>% subset(db == cur_db & upregulated == FALSE) %>% top_n(nterms, wt=Combined.Score)
      plot_df$Term <- do.call(rbind, strsplit(plot_df$Term, '[(]GO:'))[,1]
      plot_df$wrap <- wrapText(plot_df$Term, 50)


      p2 <- ggplot(plot_df, aes(x=Combined.Score, y=reorder(wrap, Combined.Score)))+
      geom_bar(stat='identity', position='identity', color='white', fill=cur_color, alpha=0.75) +
      geom_text(aes(label=wrap), x=1, color="black", size=3.5, hjust='left') +
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) + ylab('') + xlab('Combined.Score') + ggtitle(group2)


      if(nrow(plot_df) == 0){
        p2 <- ggplot() + theme_void()
      }

      patch <- p1 + p2 + plot_annotation(title=paste(cur_group, cur_db ), theme = theme(plot.title = element_text(hjust = 0.5)))

      plot_list[[cur_db]] <- patch

    }

    if(length(plot_list) == 0){next}

    # plot height should be the number of terms / 2, with a min of 4
    if(nterms > 8){
      plot_height = nterms/2
    } else{
      plot_height = 4
    }

    pdf(paste0(opt$outdir, '/enrichr/', name, '_', cur_group, '.pdf'), width=8, height=plot_height)
    for(plot in plot_list){
      print(plot)
    }
    dev.off()

  }

} else{ print('Skipping EnrichR') }
