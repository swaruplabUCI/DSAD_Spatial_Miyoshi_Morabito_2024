#!bin/bash

# directory containing cellranger outputs
splitseq_dir=$1

# output
outdir=$2

mkdir $outdir

# copy web reports:
for dir in $splitseq_dir/*
do
  name=$(basename $dir)
  cp $dir"/all-well_analysis_summary.html" $outdir/$name"_all-well_analysis_summary.html"
done
