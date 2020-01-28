#!/bin/bash

# Download and unzip fasta and gtf
wget ftp://ftp.ensembl.org/pub/release-98/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/gtf/gallus_gallus/Gallus_gallus.GRCg6a.98.gtf.gz
gunzip Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
gunzip Gallus_gallus.GRCg6a.98.gtf.gz

# make a filtered version of the gtf without any pseudogenes etc.
cellranger mkgtf Gallus_gallus.GRCg6a.98.gtf Gallus_gallus.cellranger.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:IG_V_gene

# make a reference for cellranger with STAR indices etc.
cellranger mkref --nthreads 16 --genome GalGal5_cellranger \
    --fasta Gallus_gallus.GRCg6a.dna.toplevel.fa.gz \
    --genes Gallus_gallus.cellranger.gtf

# make a reference for single-nucleus that includes pre-mRNA
awk 'BEGIN {FS="\t"; OFS="\t"} $3 == "transcript" { $3 = "exon"; print }' \
    GalGal5_cellranger/genes/genes.gtf > GalGal5.premRNA.gtf
cellranger mkref --nthreads 16 --genome GalGal5_premRNA \
    --fasta GalGal5_cellranger/fasta/genome.fa \
    --genes GalGal5.premRNA.gtf
