#!/bin/env python3
"""
Prepend "MT-" to mitochondrial gene names (if not already present) in
a gtf for cases when the mitochondrial chromosome is named "MT"
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
        if record.seqname.upper() == 'MT' and record.attributes is not None:
            if ('gene_name' in record.attributes
                    and not record.attributes['gene_name'].startswith('MT-')):
                record.attributes['gene_name'] = (
                    'MT-' + record.attributes['gene_name']
                )
        elif record.seqname.upper() == 'PT' and record.attributes is not None:
            if ('gene_name' in record.attributes
                    and not record.attributes['gene_name'].startswith('PT-')):
                record.attributes['gene_name'] = (
                    'PT-' + record.attributes['gene_name']
                )
        print(record)


# command line interface (making this a modulino)
if __name__ == '__main__':
    main(parse_args())

