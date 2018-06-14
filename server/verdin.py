#! /usr/bin/env python

from __future__ import print_function
import sys
import csv
import collections
import subprocess
import argparse
import json
import numpy
import os


def revcpl(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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

def variantsToPrimer3(params, vrs):
    with open(params['fnPrimer3In'], 'w') as f:
        for idname in vrs.keys():
            t = vrs[idname]
            if (t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("SNP")) or (t['type'].startswith("INS")) or (t['type'] == "BND_3to5"):
                ls = max(0, t['pos1'] - params['spacer'] - params['pcrLen'])
                le = max(ls + 1, t['pos1'] - params['spacer'])
                seq1 = localRef(params['genome'], t['chr1'], ls, le)
                ls = t['pos2'] + params['spacer']
                le = t['pos2'] + params['spacer'] + params['pcrLen']
                seq2 = localRef(params['genome'], t['chr2'], ls, le)
            elif (t['type'] == "INV_3to3") or (t['type'] == "BND_3to3"):
                ls = max(0, t['pos1'] - params['spacer'] - params['pcrLen'])
                le = max(ls + 1, t['pos1'] - params['spacer'])
                seq1 = localRef(params['genome'], t['chr1'], ls, le)
                ls = t['pos2'] - params['spacer'] - params['pcrLen']
                if t['type'] == "INV_3to3":
                    ls = max(t['pos1'] + params['spacer'], ls)
                le = max(ls + 1, t['pos2'] - params['spacer'])
                seq2 = revcpl(localRef(params['genome'], t['chr2'], ls, le))
            elif (t['type'] == "INV_5to5") or (t['type'] == "BND_5to5"):
                ls = t['pos1'] + params['spacer']
                le = t['pos1'] + params['spacer'] + params['pcrLen']
                if t['type'] == "INV_5to5":
                    le = max(ls + 1, min(t['pos2'] - params['spacer'], le))
                seq1 = revcpl(localRef(params['genome'], t['chr1'], ls, le))
                ls = t['pos2'] + params['spacer']
                le = t['pos2'] + params['spacer'] + params['pcrLen']
                seq2 = localRef(params['genome'], t['chr2'], ls, le)
            elif (t['type'].startswith("DUP")) or (t['type'] == "BND_5to3"):
                ls = t['pos1'] + params['spacer']
                le = t['pos1'] + params['spacer'] + params['pcrLen']
                if t['type'].startswith("DUP"):
                    le = max(ls + 1, min(t['pos2'] - params['spacer'], le))
                seq1 = revcpl(localRef(params['genome'], t['chr1'], ls, le))
                ls = t['pos2'] - params['spacer'] - params['pcrLen']
                if t['type'].startswith("DUP"):
                    ls = max(t['pos1'] + params['spacer'], ls)
                le = max(ls + 1, t['pos2'] - params['spacer'])
                seq2 = revcpl(localRef(params['genome'], t['chr2'], ls, le))
            else:
                print("Unknown variant type:", t['type'], file=sys.stderr)
                quit()
            primer3Input(params, seq1, seq2, idname, f)
        f.close()

    
def primer3(prs, vrs):
    variantsToPrimer3(prs, vrs)
    subprocess.call(['primer3_core', '--output=' + prs['fnPrimer3Out'], prs['fnPrimer3In']])
    os.remove(prs['fnPrimer3In'])

def primer3ToSilica(params):
    with open(params['fnSilicaIn'], 'w') as f:
        with open(params['fnPrimer3Out'], 'rb') as keyvalf:
            seqreader = csv.reader(keyvalf, delimiter='=')
            seqid = "NA"
            for row in seqreader:
                if row[0] == 'SEQUENCE_ID':
                    seqid = row[1]
                if row[0].endswith("_SEQUENCE"):
                    print(">" + seqid + "_" + row[0], file=f)
                    print(row[1], file=f)
    os.remove(params['fnPrimer3Out'])
    
def silicaStrict(params):
    primer3ToSilica(params)
    subprocess.call(['silica', '-q', '-m', '5', '-d', '0', '-f', 'json', '-p', params['fnSilicaOut1'], '-o', params['fnSilicaOut2'], '-g', params['genome'], params['fnSilicaIn']])


def silica(params):
    subprocess.call(['silica', '-f', 'json', '-c', '50', '-p', params['fnSilicaOut1'], '-o', params['fnSilicaOut2'], '-g', params['genome'], params['fnSilicaIn']])

def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()
    

def primerDesign(filename, genome, prefix):
    # Parameters
    params = {'fnPrimer3In': prefix + ".primer3.input", 'fnPrimer3Out': prefix + ".primer3.output", 'fnSilicaIn': prefix + ".silica.input", 'fnSilicaOut1': prefix + ".silica.primer.output", 'fnSilicaOut2': prefix + ".silica.amplicon.output", 'genome': genome, 'pcrLen': 500, 'spacer': 50, 'nprimer': 100, 'PRIMER_MAX_TM': 65, 'PRIMER_MIN_TM': 56, 'PRIMER_OPT_TM': 62, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 20, 'PRIMER_MAX_SIZE': 24}

    # Load variants
    variants = collections.defaultdict(dict)
    idcounter = 0
    with open(filename, 'rb') as csvfile:
        prireader = csv.DictReader(csvfile, delimiter='\t')
        for row in prireader:
            if not row['chr1'].startswith("#"):
                variants[idcounter] = {'chr1': row['chr1'], 'pos1': int(row['pos1']), 'chr2': row['chr2'], 'pos2': int(row['pos2']), 'type': row['type']}
                idcounter += 1
                
    # Generate primer candidates
    primer3(params, variants)

    # Silica strict primer pruning
    silicaStrict(params)

    # Primer count table
    pritable = numpy.full((idcounter * params['nprimer'], 2), 0)
    priseq = dict()
    with open(params['fnSilicaOut1']) as prijson:
        for line in prijson:
            if line.startswith('{'):
                line = line.strip().rstrip(',')
                d = json.loads(line)
                fields = d["Name"].split('_')
                idname = int(fields[0])
                lr = fields[2]
                candidate = int(fields[3])
                if lr == "LEFT":
                    if pritable[idname * params['nprimer'] + candidate, 0] == 0:
                        priseq[(idname * params['nprimer'] + candidate, 0)] = d["Seq"]
                    pritable[idname * params['nprimer'] + candidate, 0] += 1
                elif lr == "RIGHT":
                    if pritable[idname * params['nprimer'] + candidate, 1] == 0:
                        priseq[(idname * params['nprimer'] + candidate, 1)] = d["Seq"]
                    pritable[idname * params['nprimer'] + candidate, 1] += 1
                    
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
        print("Primer selection for", sum(vartodo), "variants")
        varincl = numpy.full(idcounter, False)
        with open(params['fnSilicaIn'], 'w') as f:
            for (score, idname, candidate) in scorelst:
                if used[idname * params['nprimer'] + candidate]:
                    continue
                if (vartodo[idname]) and (not varincl[idname]):
                    varincl[idname] = True
                    print(">" + str(idname) + "_PRIMER_LEFT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 0)], file=f)
                    print(">" + str(idname) + "_PRIMER_RIGHT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 1)], file=f)
                    used[idname * params['nprimer'] + candidate] = True
                    if sum(varincl) == sum(vartodo):
                        break

        if not sum(varincl):
            print("Leftover variants no primers were found", sum(vartodo))
            for idname in range(idcounter):
                if vartodo[idname]:
                    print(variants[idname])
            return primerlst

        # Run Silica
        silica(params)

        # Safely discard all primer pairs with >1 amplicon
        ampct = numpy.full(idcounter * params['nprimer'], 0)
        ampsmall = numpy.full(idcounter * params['nprimer'], 0)
        with open(params['fnSilicaOut2']) as ampjson:
            for line in ampjson:
                if line.startswith('{'):
                    line = line.strip().rstrip(',')
                    d = json.loads(line)
                    forname = d['ForName'].replace("_LEFT_", "_")
                    revname = d['RevName'].replace("_RIGHT_", "_")
                    if forname == revname:
                        fields = forname.split('_')
                        idname = int(fields[0])
                        candidate = int(fields[2])
                        # For small variants we may have one amplicon overlapping the variant
                        t = variants[idname]
                        if (t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("SNP")) or (t['type'].startswith("INS")):
                            if (d['Chrom'] == t['chr1']) and (d['ForPos'] < t['pos1']) and (d['RevPos'] > t['pos2']):
                                ampsmall[idname * params['nprimer'] + candidate] += 1
                                # Multiple amplicons over a small variant or only one (accepted)
                                if ampsmall[idname * params['nprimer'] + candidate] == 1:
                                    continue
                        ampct[idname * params['nprimer'] + candidate] += 1
        
        # Get primer and amplicon counts
        pricount = collections.Counter()
        prleftjs = dict()
        prrightjs = dict()
        with open(params['fnSilicaOut1']) as prijson:
            for line in prijson:
                if line.startswith('{'):
                    line = line.strip().rstrip(',')
                    d = json.loads(line)
                    if int(d['Tm']) < params['PRIMER_MIN_TM']:
                        continue
                    if int(d['Tm']) > params['PRIMER_MAX_TM']:
                        continue
                    pricount[d['Name']] += 1
                    fields = d["Name"].split('_')
                    idname = int(fields[0])
                    lr = fields[2]
                    candidate = int(fields[3])
                    if ampct[idname * params['nprimer'] + candidate] == 0:
                        t = variants[idname]
                        if lr == "LEFT":
                            if (t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("SNP")) or (t['type'].startswith("INS")) or (t['type'] == "BND_3to5"):
                                if (d['Chrom'] == t['chr1']) and (d['Pos'] < t['pos1']) and (d['Pos'] + params['pcrLen'] > t['pos1']) and (d['Ori'] == "forward"):
                                    prleftjs[idname] = d
                            elif (t['type'] == "INV_3to3") or (t['type'] == "BND_3to3"):
                                if (d['Chrom'] == t['chr1']) and (d['Pos'] < t['pos1']) and (d['Pos'] + params['pcrLen'] > t['pos1']) and (d['Ori'] == "forward"):
                                    prleftjs[idname] = d
                            elif (t['type'] == "INV_5to5") or (t['type'] == "BND_5to5"):
                                if (d['Chrom'] == t['chr1']) and (d['Pos'] > t['pos1']) and (d['Pos'] < t['pos1'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prleftjs[idname] = d
                            elif (t['type'].startswith("DUP")) or (t['type'] == "BND_5to3"):
                                if (d['Chrom'] == t['chr1']) and (d['Pos'] > t['pos1']) and (d['Pos'] < t['pos1'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prleftjs[idname] = d
                        if lr == "RIGHT":
                            if (t['type'].startswith("DEL")) or (t['type'].startswith("SNV")) or (t['type'].startswith("SNP")) or (t['type'].startswith("INS")) or (t['type'] == "BND_3to5"):
                                if (d['Chrom'] == t['chr2']) and (d['Pos'] > t['pos2']) and (d['Pos'] < t['pos2'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prrightjs[idname] = d
                            elif (t['type'] == "INV_3to3") or (t['type'] == "BND_3to3"):
                                if (d['Chrom'] == t['chr2']) and (d['Pos'] < t['pos2']) and (d['Pos'] + params['pcrLen'] > t['pos2']) and (d['Ori'] == "forward"):
                                    prrightjs[idname] = d
                            elif (t['type'] == "INV_5to5") or (t['type'] == "BND_5to5"):
                                if (d['Chrom'] == t['chr2']) and (d['Pos'] > t['pos2']) and (d['Pos'] < t['pos2'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prrightjs[idname] = d
                            elif (t['type'].startswith("DUP")) or (t['type'] == "BND_5to3"):
                                if (d['Chrom'] == t['chr2']) and (d['Pos'] < t['pos2']) and (d['Pos'] + params['pcrLen'] > t['pos2']) and (d['Ori'] == "forward"):
                                    prrightjs[idname] = d

        # Check primers
        for pleft in pricount.keys():
            if "_LEFT_" in pleft:
                pright = pleft.replace("_LEFT_", "_RIGHT_")
                amp = pleft.replace("_LEFT_", "_")
                fields = amp.split('_')
                idname = int(fields[0])
                candidate = int(fields[2])
                if idname in prleftjs:
                    if idname in prrightjs:
                        if ampct[idname * params['nprimer'] + candidate] == 0:
                            vartodo[idname] = False
                if not vartodo[idname]:
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
    plst = primerDesign(args.variants, args.genome, args.prefix)
    verdinOut = args.prefix + ".verdin"
    with open(verdinOut, 'w') as vout:
        json.dump(plst, vout)
