
# source activate scrublet
import scanpy as sc
import scrublet as scr
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

os.chdir('/dfs3b/swaruplab/smorabit/analysis/ADDS_2021/doublet_detection')

# path to split-seq directory:
splitseq_dir = "/dfs3b/swaruplab/smorabit/data/ADDS_2021/analysis/"

# batch directories:
batch_dirs = ["{}{}/".format(splitseq_dir, d) for d in os.listdir(splitseq_dir)]

# output directories
fig_dir = 'figures/'
data_dir = 'data/'

#
big_adata_list = []
for i, batch in enumerate(os.listdir(splitseq_dir)):
  print(batch)
  adata_list = []
  for sublib in os.listdir(batch_dirs[i]):
    infile = '{}{}/all-well/DGE_unfiltered/cellbender_filtered.h5'.format(batch_dirs[i], sublib)
    cur_genes = pd.read_table('{}{}/all-well/DGE_unfiltered/genes.tsv'.format(batch_dirs[i], sublib), sep='\t', header=None)
    print(infile)
    cur_adata = sc.read_10x_h5(infile, genome='background_removed')
    cur_adata.var_names_make_unique()
    # add metadata
    meta_infile = '{}{}/all-well/DGE_unfiltered/cell_metadata.csv'.format(batch_dirs[i], sublib)
    cur_meta = pd.read_table(meta_infile, sep=',')
    cur_meta.index = cur_meta['cell_barcode']
    cur_meta.index.name = None
    cur_meta = cur_adata.obs.merge(cur_meta, how='left', left_index=True, right_on='cell_barcode')
    cur_adata.obs = cur_meta
    cur_adata.obs['Sublibrary'] = sublib
    cur_adata.obs['Batch'] = batch
    # add gene data
    cur_genes = cur_adata.var.merge(cur_genes, how='left', left_index=True, right_on=1)
    cur_adata.var.index = cur_genes[2]
    cur_adata.var['ensembl_id'] = list(cur_genes[1])
    cur_adata.var.index.name = None
    cur_adata.var = cur_adata.var.drop('gene_ids', axis=1)
    cur_adata.var_names_make_unique()
    adata_list.append(cur_adata)
  # merge list of anndata objects:
  adata = adata_list[0].concatenate(adata_list[1:])
  sc.pp.filter_cells(adata, min_counts = 100)
  scrub = scr.Scrublet(adata.X)
  adata.obs['doublet_scores'], predicted_doublets = scrub.scrub_doublets()
  adata.obs['doublets'] = np.where(predicted_doublets== False, 'Singlet', 'Doublet')
  ax = scrub.plot_histogram()
  plt.savefig('{}{}_scrublet_histogram.pdf'.format(fig_dir, batch))
  adata.obs.to_csv('{}{}_scrublet.csv'.format(data_dir, batch))
  big_adata_list.append(adata)

# concatenate adatas:
adata = big_adata_list[0].concatenate(big_adata_list[1:])

# add the sample metadata:
adata.obs['SampleID'] = list(adata.obs['sample'].astype(str))
sample_meta = pd.read_table("~/swaruplab/smorabit/analysis/ADDS_2021/data/AD-DS_Cases.csv", sep=',')
sample_meta['SampleID'] = list(sample_meta['SampleID'].astype(str))
sample_meta_cells = adata.obs.merge(sample_meta, how='left', left_on = 'SampleID', right_on = 'SampleID')
sample_meta_cells.index = adata.obs.index
sample_meta_cells = sample_meta_cells.drop(['Unnamed: 0', 'sample', 'Unnamed: 20'], axis=1)
adata.obs = sample_meta_cells

# save the output:
adata.shape
adata.var = adata.var[['ensembl_id-0-0']]
adata.var['gene_name'] = adata.var.index
adata.write_h5ad('{}/ADDS_cellbender_scrublet_unfiltered.h5ad'.format(data_dir))
