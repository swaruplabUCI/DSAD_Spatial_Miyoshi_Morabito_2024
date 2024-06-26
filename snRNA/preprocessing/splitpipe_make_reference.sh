#!/bin/bash

# download reference:
wget ftp://​ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://​ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz

# load required software
conda activate splitseq
module load samtools
module load star

# make reference
splitpipe/split-pipe \
--mode mkref \
--genome hg38 \
--fasta /dfs3b/swaruplab/smorabit/data/ADDS_2021/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes /dfs3b/swaruplab/smorabit/data/ADDS_2021/genomes/Homo_sapiens.GRCh38.93.gtf.gz \
--output_dir ./genomes/hg38
