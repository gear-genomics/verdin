#! /usr/bin/env python

from __future__ import print_function
import os
import uuid
import datetime
import csv
import collections
import subprocess 
import argparse


def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()

def primerDesign(filename, genome):
    # Parameters
    sangerLen = 500
    pcrLen = 5000
    primerToVariantSpacer = 50

    # Load variants
    variants = collections.defaultdict(list)
    with open(filename, 'rb') as csvfile:
        prireader = csv.DictReader(csvfile, delimiter='\t')
        for row in prireader:
            variants[row['type']].append((row['chr1'], int(row['pos1']), row['chr2'], int(row['pos2'])))

    # Run primer design for each variant type
    for t in variants.keys():
        if ((t.startswith("DEL")) or (t.startswith("SNV")) or (t.startswith("INS"))):
            for v in variants[t]:
                ls = max(0, v[1] - primerToVariantSpacer - pcrLen)
                le = max(ls + 1, v[1] - primerToVariantSpacer)
                seq1 = localRef(genome, v[0], ls, le)
                ls = v[3] + primerToVariantSpacer
                le = v[3] + primerToVariantSpacer + pcrLen
                seq2 = localRef(genome, v[2], ls, le)
                print(seq1, seq2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verdin')
    parser.add_argument('-v', '--variants', required=True, metavar="variants.tsv", dest='variants', help='input variants')
    parser.add_argument('-g', '--genome', required=True, metavar="genome.fa.gz", dest='genome', help='input genome')
    args = parser.parse_args()
    primerDesign(args.variants, args.genome)
