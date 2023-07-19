#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/smorabit/analysis/ADDS_2021/

# get a list of all samples:
#splitseq_files="sublibrary_dirs.csv"
splitseq_files="Batch5_sublibrary_dirs.csv"

#for i in {1..32}; do
for i in {1..8}; do

    # get sublibrary directory
    sublib_dir=$(head -$i $splitseq_files | tail -1)
    echo $sublib_dir
    
    # set input and output file
    infile=$sublib_dir
    outfile=$sublib_dir/cellbender.h5

    # run cellbender
    cellbender remove-background \
       --input $infile \
       --output $outfile \
       --expected-cells 10000 \
       --total-droplets-included 25000 \
       --epochs 150 \
       --cuda

done;



# test one:

# test run cellbender on split-seq data:
# potential problem: is the matrix the wrong orientation???
# from the cellbender log:
# Including 252018 genes that have nonzero counts.
#
# I made a new .mtx that is flipped
# cellbender remove-background \
#    --input test_cellbender/Sublibrary1_S1/all-well/DGE_unfiltered/ \
#    --output test_cellbender/Sublibrary1_S1/all-well/DGE_unfiltered/cellbender.h5 \
#    --expected-cells 10000 \
#    --total-droplets-included 25000 \
#    --epochs 150 \
#    --cuda
