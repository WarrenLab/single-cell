#!/bin/env python3
"""
If gene name is simply totally missing, use the gene ID for the gene name.
(Don't confuse this with the case in which there is actualy a gene name
embeded in the gene id, for which you would want to use the utility
'extract_name_embedded_in_gene_id.py'.)
"""

import argparse
import gtfez


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'gtf',
        type=argparse.FileType('r'),
        help='location of GTF file',
    )

    args = parser.parse_args()
    return args


def main(args):
    for line in args.gtf:
        if line.startswith('#'):
            print(line.rstrip())
            continue

        record = gtfez.GTFRecord(line)
        if ('gene_name' not in record.attributes) and 'gene_id' in record.attributes:
            record.attributes['gene_name'] = record.attributes['gene_id']
        print(record)


# command line interface (making this a modulino)
if __name__ == '__main__':
    main(parse_args())
