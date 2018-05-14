#! /usr/bin/env python

from __future__ import print_function
import os
import uuid
import datetime
import csv
import collections
import subprocess 
import argparse


def primer3Input(params, seq1, seq2, idcounter, f):
    nrun = 'N'*(2 * params['spacer'])
    print("SEQUENCE_ID=", idcounter, sep='', file=f)
    print("SEQUENCE_TEMPLATE=", seq1, nrun, seq2, sep='', file=f)
    print("SEQUENCE_TARGET=", len(seq1), ",", len(nrun), sep='', file=f)
    print("PRIMER_PRODUCT_SIZE_RANGE=", len(nrun), "-", len(seq1+nrun+seq2), sep='', file=f)
    print("PRIMER_TASK=generic", file=f)
    print("PRIMER_PICK_LEFT_PRIMER=1", file=f)
    print("PRIMER_PICK_INTERNAL_OLIGO=0", file=f)
    print("PRIMER_PICK_RIGHT_PRIMER=1", file=f)
    print("PRIMER_MAX_NS_ACCEPTED=1", file=f)
    print("PRIMER_NUM_RETURN=", params['nprimer'], sep='', file=f)
    print("PRIMER_EXPLAIN_FLAG=1", file=f)
    print("PRIMER_MAX_TM=", params['PRIMER_MAX_TM'], sep='', file=f)
    print("PRIMER_MIN_TM=", params['PRIMER_MIN_TM'], sep='', file=f)
    print("PRIMER_OPT_TM=", params['PRIMER_OPT_TM'], sep='', file=f)
    print("PRIMER_OPT_SIZE=", params['PRIMER_OPT_SIZE'], sep='', file=f)
    print("PRIMER_MIN_SIZE=", params['PRIMER_MIN_SIZE'], sep='', file=f)
    print("PRIMER_MAX_SIZE=", params['PRIMER_MAX_SIZE'], sep='', file=f)
    print("=", file=f)

def primer3(primer3In, primer3Out):
    subprocess.call(['primer3_core', '--output=' + primer3Out, primer3In])

def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()

def primerDesign(filename, genome, prefix):
    # Parameters
    #params = { 'sangerLen': 500, 'pcrLen': 5000, 'spacer': 50, 'nprimer': 1000, 'PRIMER_MAX_TM': 65}
    params = { 'sangerLen': 500, 'pcrLen': 200, 'spacer': 50, 'nprimer': 5, 'PRIMER_MAX_TM': 65, 'PRIMER_MIN_TM': 56, 'PRIMER_OPT_TM': 62, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 20, 'PRIMER_MAX_SIZE': 24 }

    # Load variants
    variants = collections.defaultdict(list)
    idcounter = 1
    with open(filename, 'rb') as csvfile:
        prireader = csv.DictReader(csvfile, delimiter='\t')
        for row in prireader:
            variants[row['type']].append((row['chr1'], int(row['pos1']), row['chr2'], int(row['pos2']), idcounter))
            idcounter += 1

    # Generate Primer3 input
    primer3In = prefix + ".primer3.input"
    with open(primer3In, 'w') as f:
        for t in variants.keys():
            if ((t.startswith("DEL")) or (t.startswith("SNV")) or (t.startswith("INS"))):
                for v in variants[t]:
                    ls = max(0, v[1] - params['spacer'] - params['pcrLen'])
                    le = max(ls + 1, v[1] - params['spacer'])
                    seq1 = localRef(genome, v[0], ls, le)
                    ls = v[3] + params['spacer']
                    le = v[3] + params['spacer'] + params['pcrLen']
                    seq2 = localRef(genome, v[2], ls, le)
                    primer3Input(params, seq1, seq2, v[4], f)
        f.close()

    # Generate primer candidates
    primer3Out = prefix + ".primer3.output"
    primer3(primer3In, primer3Out)

    # Silica input
    silicaIn = prefix + ".silica.input"
    with open(silicaIn, 'w') as f:
        with open(primer3Out, 'rb') as keyvalf:
            seqreader = csv.reader(keyvalf, delimiter='=')
            seqid = "NA"
            for row in seqreader:
                if row[0] == 'SEQUENCE_ID':
                    seqid = row[1]
                if row[0].endswith("_SEQUENCE"):
                    print(">" + seqid + "_" + row[0], file=f)
                    print(row[1], file=f)
        
    
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verdin')
    parser.add_argument('-v', '--variants', required=True, metavar="variants.tsv", dest='variants', help='input variants')
    parser.add_argument('-g', '--genome', required=True, metavar="genome.fa.gz", dest='genome', help='input genome')
    parser.add_argument('-p', '--prefix', required=True, metavar="outprefix", dest='prefix', help='output prefix')
    args = parser.parse_args()
    primerDesign(args.variants, args.genome, args.prefix)
