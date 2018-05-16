#! /usr/bin/env python

from __future__ import print_function
from flask import jsonify
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
    subprocess.call(['silica', '-q', '-m', '5', '-d', '0', '-f', 'json', '-p', silicaOut1, '-o', silicaOut2, '-g', genome, silicaIn])

def silica(genome, silicaIn, silicaOut1, silicaOut2):
    subprocess.call(['silica', '-f', 'json', '-c', '50', '-p', silicaOut1, '-o', silicaOut2, '-g', genome, silicaIn])
    
def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()

def primerDesign(filename, genome, prefix):
    # Parameters
    params = { 'pcrLen': 500, 'spacer': 50, 'nprimer': 100, 'PRIMER_MAX_TM': 65, 'PRIMER_MIN_TM': 56, 'PRIMER_OPT_TM': 62, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 20, 'PRIMER_MAX_SIZE': 24 }

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
    primerlst = list()
    scorelst = sorted(scorelst)
    vartodo = numpy.full(idcounter, True)
    used = numpy.full(idcounter * params['nprimer'], False)
    while sum(vartodo):
        print("Primer selection for", sum(vartodo), "variants,", len(scorelst) - sum(used), " primer pairs left for checking.")
        varincl = numpy.full(idcounter, False)
        with open(silicaIn, 'w') as f:
            for (score, idname, candidate) in scorelst:
                if sum(varincl) == sum(vartodo):
                    break
                if used[idname * params['nprimer'] + candidate]:
                    continue
                if (vartodo[idname]) and (not varincl[idname]):
                    varincl[idname] = True
                    print(">" + str(idname) + "_PRIMER_LEFT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 0)], file=f)
                    print(">" + str(idname) + "_PRIMER_RIGHT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 1)], file=f)
                    used[idname * params['nprimer'] + candidate] = True
        if not sum(varincl):
            print("Leftover variants no primers were found", sum(vartodo))
            return primerlst

        # Run Silica
        silica(genome, silicaIn, silicaOut1, silicaOut2)

        # Get primer and amplicon counts
        pricount = collections.Counter()
        prleftjs = dict()
        prrightjs = dict()
        with open(silicaOut1) as prijson:
            for line in prijson:
                if line.startswith('{'):
                    line = line.strip().rstrip(',')
                    d = json.loads(line)
                    if (int(d['Tm']) >= params['PRIMER_MIN_TM']) and (int(d['Tm']) <= params['PRIMER_MAX_TM']):
                        pricount[d['Name']] += 1
                        fields = d["Name"].split('_')
                        idname = int(fields[0])
                        lr = fields[2]
                        candidate = int(fields[3])
                        t = variants[idname]
                        if lr == "LEFT":
                            if ((t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("INS"))):
                                if (d['Chrom'] == t['chr1']) and (d['Pos'] < t['pos1']) and (d['Pos'] + params['pcrLen'] > t['pos1']):
                                    prleftjs[idname] = d
                        if lr == "RIGHT":
                            if ((t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("INS"))):
                                if (d['Chrom'] == t['chr2']) and (d['Pos'] > t['pos2']) and (d['Pos'] < t['pos2'] + params['pcrLen']):
                                    prrightjs[idname] = d

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
                        
        # Check primers
        for pleft in pricount.keys():
            if "_LEFT_" in pleft:
                pright = pleft.replace("_LEFT_", "_RIGHT_")
                amp = pleft.replace("_LEFT_", "_")
                idname = int(pleft.split('_')[0])
                if idname in prleftjs:
                    if idname in prrightjs:
                        t = variants[idname]
                        if ((t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("INS"))):
                            if ampcount[amp] == 1:
                                if (ampjs[idname]['Chrom'] == t['chr1']) and (ampjs[idname]['ForPos'] < t['pos1']) and (ampjs[idname]['RevPos'] > t['pos2']):
                                    vartodo[idname] = False
                                elif ampcount[amp] == 0:
                                    vartodo[idname] = False
                if vartodo[idname]:
                    if 0:
                        print(t)
                        print(pricount[pleft], pricount[pright], ampcount[amp])
                        print(prleftjs[idname])
                        print(prrightjs[idname])
                        if ampcount[amp]:
                            print(ampjs[idname])
                else:
                    out = t
                    out['Primer1Chrom'] = prleftjs[idname]['Chrom']
                    out['Primer1Pos'] = prleftjs[idname]['Pos']
                    out['Primer1Seq'] = prleftjs[idname]['Seq']
                    out['Primer1Tm'] = prleftjs[idname]['Tm']
                    out['Primer1Ori'] = prleftjs[idname]['Ori']
                    out['Primer1Hits'] = pricount[pleft]
                    out['Primer2Chrom'] = prrightjs[idname]['Chrom']
                    out['Primer2Pos'] = prrightjs[idname]['Pos']
                    out['Primer2Seq'] = prrightjs[idname]['Seq']
                    out['Primer2Tm'] = prrightjs[idname]['Tm']
                    out['Primer2Ori'] = prrightjs[idname]['Ori']
                    out['Primer2Hits'] = pricount[pright]
                    primerlst.append(out)
    return primerlst
                        
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verdin')
    parser.add_argument('-v', '--variants', required=True, metavar="variants.tsv", dest='variants', help='input variants')
    parser.add_argument('-g', '--genome', required=True, metavar="genome.fa.gz", dest='genome', help='input genome')
    parser.add_argument('-p', '--prefix', required=True, metavar="outprefix", dest='prefix', help='output prefix')
    args = parser.parse_args()
    primerlst = primerDesign(args.variants, args.genome, args.prefix)
    verdinOut = args.prefix + ".verdin"
    with open(verdinOut, 'w') as vout:
        json.dump(primerlst, vout)

