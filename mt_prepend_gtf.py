#!/bin/env python3

import argparse

def main(args):
    with open(args.gtf) as gtf_fh:
        for line in gtf_fh:
            line = line.rstrip()
            if line.startswith('MT'):
                if 'gene_name "' in line:
                    if 'gene_name "MT-' not in line:
                        line = line.replace('gene_name "','gene_name "MT-')
            print(line)

# command line interface (making this a modulino)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                description='Prepend "MT-" to mitochondrial gene names (if not already present) for cases when the mitochondrial chromosome is named "MT"'
             )
    parser.add_argument(
        'gtf',
        type=str,
        help='Name of GTF file',
    )

    args = parser.parse_args()

    main(args)
