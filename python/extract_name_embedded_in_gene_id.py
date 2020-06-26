#!/bin/env python3
"""
Extract gene name from gene ID if gene name not defined and if gene
ID starts with "gene-"
"""
import argparse

import gtfez


def parse_args():
    """argument parser for main method"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'gtf',
        type=argparse.FileType('r'),
        help='location of GTF file',
    )

    args = parser.parse_args()
    return args


def main(args):
    """main method of script"""
    for line in args.gtf:
        if line.startswith('#'):
            print(line.rstrip())
            continue

        record = gtfez.GTFRecord(line)
        if ('gene_name' not in record.attributes
                and 'gene_id' in record.attributes
                and record.attributes['gene_id'].lower().startswith('gene-')):
            record.attributes['gene_name'] = record.attributes['gene_id'][5:]
        print(record)


# command line interface (making this a modulino)
if __name__ == '__main__':
    main(parse_args())

