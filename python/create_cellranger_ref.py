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
import subprocess
import sys

import mkref


def parse_args():
    """parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gtf', type=str, help='name of GTF file from ensembl')
    group.add_argument(
        '--gff', type=str, help='name of GFF3 file to be converted'
    )

    parser.add_argument(
        '--fasta',
        type=str,
        help='name of reference FASTA file to be converted',
        required=True,
    )
    parser.add_argument(
        '--cpus', type=int, help='cpu cores to request', default=24
    )
    parser.add_argument(
        '--mem', type=str, help='Total RAM to allocate in GB', default='300'
    )
    parser.add_argument(
        '--partition',
        type=str,
        help='partition to use',
        default='Lewis,hpc5,hpc6,BioCompute',
    )
    parser.add_argument('--genome', type=str, help='name of genome')
    parser.add_argument(
        '--no-run',
        help="Do not run script immediately after creating it.",
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
        '--account', help="SLURM account to use", default='warrenlab'
    )
    parser.add_argument(
        '--time', help='SLURM time limit', default='1-00:00:00'
    )

    return parser.parse_args()


def main():
    """main method of script"""
    args = parse_args()

    # these are the parameters that are passed directly to the sbatch
    # executable to tell it, e.g., what resources to reserve
    sbatch_params = {
        'time': args.time,
        'cpus-per-task': args.cpus,
        'mem': f'{args.mem}G',
        'partition': args.partition,
        'account': args.account,
        'job-name': f'mkref_{args.genome}',
        'out': '%x_%j.out',
        'err': '%x_%j.err',
    }
    # these are the parameters that are passed to the sbatch script
    # through the sbatch executable; they end up as bash variables of
    # the same name that can be accessed from within the bash script
    sbatch_exports = {
        'genome': args.genome + '.MT-renamed',
        'fasta': args.fasta,
    }

    if args.gtf is not None:
        sbatch_exports['gtf'] = args.gtf
    else:
        sbatch_exports['gff'] = args.gff

    if args.pre_mrna:
        sbatch_params['job_name'] += ".pre-mRNA"
        sbatch_exports['genome'] += ".pre-mRNA"
        sbatch_exports['pre_mRNA'] = 'yes'

    command = mkref.make_sbatch_command(sbatch_params, sbatch_exports)
    print(command)

    if not args.no_run:
        subprocess.run(
            command, check=True, stdout=sys.stdout, stderr=sys.stderr
        )


if __name__ == '__main__':
    main()
