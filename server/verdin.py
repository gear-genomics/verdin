#! /usr/bin/env python

from __future__ import print_function
import os
import uuid
import datetime
import csv
import collections
import subprocess 
import argparse
import json
import numpy

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

def silicaStrict(genome, silicaIn, silicaOut1, silicaOut2):
    subprocess.call(['silica', '-q', '-m', '5', '-d', '0', '-f', 'json', '-c', '55', '-p', silicaOut1, '-o', silicaOut2, '-g', genome, silicaIn])

def silica(genome, silicaIn, silicaOut1, silicaOut2):
    subprocess.call(['silica', '-f', 'json', '-c', '55', '-p', silicaOut1, '-o', silicaOut2, '-g', genome, silicaIn])
    
def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()

def primerDesign(filename, genome, prefix):
    # Parameters
    #params = { 'sangerLen': 500, 'pcrLen': 5000, 'spacer': 50, 'nprimer': 1000, 'PRIMER_MAX_TM': 65}
    params = { 'sangerLen': 500, 'pcrLen': 500, 'spacer': 50, 'nprimer': 100, 'PRIMER_MAX_TM': 65, 'PRIMER_MIN_TM': 56, 'PRIMER_OPT_TM': 62, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 20, 'PRIMER_MAX_SIZE': 24 }

    # Load variants
    variants = collections.defaultdict(dict)
    idcounter = 0
    with open(filename, 'rb') as csvfile:
        prireader = csv.DictReader(csvfile, delimiter='\t')
        for row in prireader:
            variants[idcounter] = {'chr1': row['chr1'], 'pos1': int(row['pos1']), 'chr2': row['chr2'], 'pos2': int(row['pos2']), 'type': row['type']}
            idcounter += 1

    # Generate Primer3 input
    primer3In = prefix + ".primer3.input"
    with open(primer3In, 'w') as f:
        for idname in variants.keys():
            t = variants[idname]
            if ((t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("INS"))):
                ls = max(0, t['pos1'] - params['spacer'] - params['pcrLen'])
                le = max(ls + 1, t['pos1'] - params['spacer'])
                seq1 = localRef(genome, t['chr1'], ls, le)
                ls = t['pos2'] + params['spacer']
                le = t['pos2'] + params['spacer'] + params['pcrLen']
                seq2 = localRef(genome, t['chr2'], ls, le)
                primer3Input(params, seq1, seq2, idname, f)
        f.close()

    # Generate primer candidates
    primer3Out = prefix + ".primer3.output"
    primer3(primer3In, primer3Out)

    # Silica strict primer pruning
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
    silicaOut1 = prefix + ".silica.primer.output"
    silicaOut2 = prefix + ".silica.amplicon.output"
    silicaStrict(genome, silicaIn, silicaOut1, silicaOut2)

    # Primer count table
    pritable = numpy.full((idcounter * params['nprimer'], 2), 0)
    priseq = dict()
    with open(silicaOut1) as prijson:
        for line in prijson:
            if line.startswith('{'):
                line = line.strip().rstrip(',')
                d = json.loads(line)
                fields = d["Name"].split('_')
                idname = int(fields[0])
                lr = fields[2]
                candidate = int(fields[3])
                if lr == "LEFT":
                    pritable[idname * params['nprimer'] + candidate, 0] += 1
                    priseq[(idname * params['nprimer'] + candidate, 0)] = d["Seq"]
                elif lr == "RIGHT":
                    pritable[idname * params['nprimer'] + candidate, 1] += 1
                    priseq[(idname * params['nprimer'] + candidate, 1)] = d["Seq"]

    # Iterate all primers
    scorelst = []
    for i in range(idcounter):
        for j in range(params['nprimer']):
            if (pritable[i * params['nprimer'] + j, 0] > 0) and (pritable[i * params['nprimer'] + j, 1] > 0):
                sumhits = pritable[i * params['nprimer'] + j, 0] + pritable[i * params['nprimer'] + j, 1]
                scorelst.append((sumhits, i, j))

    # Process primers in-order
    vartodo = numpy.full(idcounter, True)
    scorelst = sorted(scorelst)
    for (score, idname, candidate) in scorelst:
        if vartodo[idname]:
            with open(silicaIn, 'w') as f:
                print(">" + str(idname) + "_PRIMER_LEFT_" + str(candidate) + "_SEQUENCE", file=f)
                print(priseq[(idname * params['nprimer'] + candidate, 0)], file=f)
                print(">" + str(idname) + "_PRIMER_RIGHT_" + str(candidate) + "_SEQUENCE", file=f)
                print(priseq[(idname * params['nprimer'] + candidate, 1)], file=f)

            # Run Silica
            silica(genome, silicaIn, silicaOut1, silicaOut2)
            pricount = collections.Counter()
            with open(silicaOut1) as prijson:
                for line in prijson:
                    if line.startswith('{'):
                        line = line.strip().rstrip(',')
                        d = json.loads(line)
                        if (int(d['Tm']) >= params['PRIMER_MIN_TM']) and (int(d['Tm']) <= params['PRIMER_MAX_TM']):
                            pricount[d['Name']] += 1
            ampcount = collections.Counter()
            ampjs = dict()
            with open(silicaOut2) as ampjson:
                for line in ampjson:
                    if line.startswith('{'):
                        line = line.strip().rstrip(',')
                        d = json.loads(line)
                        forname = d['ForName'].replace("_LEFT_","_")
                        revname = d['RevName'].replace("_RIGHT_","_")
                        if forname == revname:
                            ampcount[forname] += 1
                            idname = int(forname.split('_')[0])
                            ampjs[idname] = d
                            
            for pleft in pricount.keys():
                if "_LEFT_" in pleft:
                    pright = pleft.replace("_LEFT_", "_RIGHT_")
                    amp = pleft.replace("_LEFT_", "_")
                    t = variants[idname]
                    if ((t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("INS"))):
                        if ampcount[amp] == 1:
                            if (ampjs[idname]['Chrom'] == t['chr1']) and (ampjs[idname]['ForPos'] < t['pos1']) and (ampjs[idname]['RevPos'] > t['pos2']):
                                #print(t)
                                #print(ampjs[idname])
                                #print(pricount[pleft])
                                #print(pricount[pright])
                                #print(ampcount[amp])
                                print(idname)
                                vartodo[idname] = False
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verdin')
    parser.add_argument('-v', '--variants', required=True, metavar="variants.tsv", dest='variants', help='input variants')
    parser.add_argument('-g', '--genome', required=True, metavar="genome.fa.gz", dest='genome', help='input genome')
    parser.add_argument('-p', '--prefix', required=True, metavar="outprefix", dest='prefix', help='output prefix')
    args = parser.parse_args()
    primerDesign(args.variants, args.genome, args.prefix)
