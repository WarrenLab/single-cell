#!/usr/bin/env python3
"""Create a single-cell genome reference using cellranger.

This follows the recommendations of 10X (see
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references),
including the optional creation of a "pre-mRNA" version of a reference
for nuclei data.

It also uses "mt_prepend_gtf.py" to prepend "MT-" to each gene name
from the "MT" chromosome, if that gene name does not already start with
those three characters.
"""

import argparse
import os
import subprocess


def local_replace_ext(file_name, replacement):
    """
    Replace a file name's extension.

    Args:
        file_name: a path-like object containing the name of the file
            whose extension will be replaced
        replacement: string containing replacement extension

    Returns:
        string containing modified filename
    """
    local_file_name = os.path.basename(file_name)
    return os.path.splitext(local_file_name)[0] + replacement


# module loading (system specific. This works on MU's Lewis cluster)
MODULE_LOAD_TEXT = """module load biocompute/biocompute-modules
module load gffread/gffread-0.11.3
source /cluster/biocompute/software/cellranger/cellranger-3.1.0/sourceme.bash"""


def parse_args():
    """parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--gtf',
        type=str,
        help='name of GTF file from ensembl',
    )
    group.add_argument(
        '--gff',
        type=str,
        help='name of GFF3 file to be converted',
    )

    parser.add_argument(
        '--fasta',
        type=str,
        help='name of reference FASTA file to be converted',
        required=True,
    )
    parser.add_argument(
        '--cpus',
        type=int,
        help='cpu cores to request',
        default=24,
    )
    parser.add_argument(
        '--mem',
        type=str,
        help='Total RAM to allocate in GB',
        default='300',
    )
    parser.add_argument(
        '--partition',
        type=str,
        help='partition to use',
        default='Lewis,hpc5,hpc6,BioCompute',
    )
    parser.add_argument(
        '--genome',
        type=str,
        help='name of genome',
    )
    parser.add_argument(
        '--no-run',
        help="Do not run script immediately after creating it.",
        action='store_false',
        dest='run',
    )
    parser.add_argument(
        '--run',
        help="Run script immediately after creating it.",
        action='store_true',
    )
    parser.add_argument(
        '--pre-mRNA',
        help='Convert GTF transcripts to behave as one large exon '
             '(thus pre-mRNA reads map to the transcriptome, may be '
             'useful for nuclei data)',
        action='store_true',
        default=False,
        dest='pre_mrna',
    )
    parser.add_argument(
        '--account',
        help="SLURM account to use",
        default='warrenlab',
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='run in debug mode',
    )

    return parser.parse_args()


def main():
    """main method of script"""
    args = parse_args()

    # Specify that the MT genes may be renamed
    genome = args.genome + ".MT-renamed"

    job_name = f"create_gtf__{args.genome}"
    if args.pre_mrna:
        job_name += ".pre-mRNA"
        genome += ".pre-mRNA"
    gtf_first = ''
    filter_command = ''
    filter_comment = ''
    gxf = ''
    if args.gtf:
        gtf_first = args.gtf
        gxf = args.gtf
    elif args.gff:
        gtf_first = local_replace_ext(args.gff, ".converted_from_gff.gtf")
        filter_comment += '# convert GFF to a GTF file and'
        filter_command += f"gffread {args.gff} -g {args.fasta} -T -o {gtf_first} "
        gxf = args.gff

    base_gxf = ""

    if args.pre_mrna:
        base_gxf = local_replace_ext(gxf, ".pre-mRNA.x")
    else:
        base_gxf = local_replace_ext(gxf, ".x")

    gtf_temp = local_replace_ext(base_gxf, ".only_lines_with_gene_id.gtf")
    gtf_mt_prepended = local_replace_ext(gtf_temp, ".MT-renamed.gtf")
    gtf_final = local_replace_ext(gtf_mt_prepended, ".filtered.gtf")

    log_final = local_replace_ext(base_gxf, ".STAR.log")

    if args.pre_mrna:
        filter_comment += (
            '# convert transcripts into exons, keeping only '
            + 'those containing "gene_id"')
        filter_command +=  \
        """awk 'BEGIN{FS="\\t"; OFS="\\t"} $3 == "transcript"{ $3="exon"; print}' """ + \
        f"{gtf_first} | grep -w gene_id > {gtf_temp} "
    else:
        filter_comment += '# keep only those GTF records containing "gene_id" '
        filter_command += f"grep -w gene_id {gtf_first} > {gtf_temp} "

    # For every gene name in the "MT" chromosome, prepend "MT-" if it
    # isn't already present
    prepend_mt_command = f"mt_prepend_gtf.py {gtf_temp} > {gtf_mt_prepended}"

    # partition here are system specific
    filled_template = """#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}G
#SBATCH --partition={partition}
#SBATCH --account={account}
#SBATCH --job-name={job_name}
#SBATCH -o {job_name}_%j_o.out

#load modules
{MODULE_LOAD_TEXT}

{filter_comment}
{filter_command}
{prepend_MT_command}

# Convert to style preferred by cellranger
#  (in some cases, it only changes LF to CRLF
cellranger mkgtf {gtf_MT_prepended} {gtf_final}           \
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

# Make reference that cellranger understands
cellranger mkref --nthreads=$SLURM_CPUS_PER_TASK   \\
                 --genome={genome}                 \\
                 --fasta={fasta}                   \\
                 --genes={gtf_final}

# Rename STAR log file
mv Log.out {log_final}
""".format(
        #module info
        MODULE_LOAD_TEXT=MODULE_LOAD_TEXT,

        #SLURM
        cpus=args.cpus,
        account=args.account,
        partition=args.partition,
        mem=args.mem,
        job_name=job_name,

        #Genome info, files
        genome=genome,
        fasta=args.fasta,
        gtf_MT_prepended=gtf_mt_prepended,
        gtf_final=gtf_final,

        #Previously build bash commands
        filter_command=filter_command,
        filter_comment=filter_comment,
        prepend_MT_command=prepend_mt_command,

        log_final=log_final
        )

    if not args.debug:
        filled_template = filled_template + \
            "# clean up"      + "\n" +      \
            "rm " + gtf_temp  + "\n" +      \
            "rm " + gtf_mt_prepended  + "\n"

    job_script_name = job_name + '.sbatch'
    with open(job_script_name, "w") as file_handle:
        file_handle.write(filled_template)

    # Run batch file (if requested)
    if args.run:
        output = subprocess.run(
            ['sbatch', job_script_name],
            stdout=subprocess.PIPE
        )
        print(output.stdout.decode('utf-8'))

# The only "filtering" really done above is keeping only lines
# containing "gene_id". You can do more specific filters. See
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

if __name__ == '__main__':
    main()
