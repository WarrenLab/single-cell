#!/bin/bash

# this script is made to be invoked by sbatch, with the following
# variables exported:
# - fasta : path to fasta file of reference genome
# - gff | gtf : annotations of reference genome
# - genome : genome name
# - pre_mRNA : set to any non-empty string for a pre-mRNA ref,
#              otherwise don't set

# load modules
module load biocompute/biocompute-modules
module load gffread/gffread-0.11.3
source /cluster/biocompute/software/cellranger/cellranger-3.1.0/sourceme.bash

# convert gff to gtf, if necessary
if [ -n "$gtf" ]; then
    gtf="${gff%.gff}.gtf"
    gffread ${gff} -g ${fasta} -T -o $gtf
fi

# convert gtf to pre-mRNA format, if applicable
if [ -z "$pre_mRNA" ]; then
    new_gtf="${gtf%.gtf}.premRNA.gtf"
    awk 'BEGIN{FS="\\t"; OFS="\\t"} $3 == "transcript"{ $3="exon"; print}' \
        < ${gtf} > ${new_gtf}
    gtf=${new_gtf}
fi

# prepend "MT-" to names of genes on chrMT, if necessary
new_gtf="${gtf%.gtf}.mt_prepended.gtf"
mt_prepend_gtf.py ${gtf} > ${new_gtf}
gtf=${new_gtf}

# Convert to style preferred by cellranger
#  (in some cases, it only changes LF to CRLF
new_gtf="${gtf%.gtf}.final.gtf"
cellranger mkgtf ${gtf} ${new_gtf}
    --attribute=gene_biotype:protein_coding  \
    --attribute=gene_biotype:lincRNA         \
    --attribute=gene_biotype:antisense       \
    --attribute=gene_biotype:IG_LV_gene      \
    --attribute=gene_biotype:IG_V_gene       \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene       \
    --attribute=gene_biotype:IG_J_gene       \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene       \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene       \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene       \
    --attribute=gene_biotype:TR_J_gene       \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene
gtf=${new_gtf}

# Make reference that cellranger understands
cellranger mkref \
    --nthreads=$SLURM_CPUS_PER_TASK \
    --genome=${genome} \
    --fasta=${fasta} \
    --genes=${gtf}

# Rename STAR log file
mv Log.out ${gtf%.gtf}.STAR.log


