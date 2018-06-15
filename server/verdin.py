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
            if vrs[idname]['done']:
                continue
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
    subprocess.call(['silica', '-q', '-m', '5', '-d', '0', '-c', '50', '-f', 'json', '-p', params['fnSilicaOut1'], '-o', params['fnSilicaOut2'], '-g', params['genome'], params['fnSilicaIn']])

def silica(params):
    subprocess.call(['silica', '-f', 'json', '-c', '50', '-p', params['fnSilicaOut1'], '-o', params['fnSilicaOut2'], '-g', params['genome'], params['fnSilicaIn']])

def localRef(genome, chrom, start, end):
    output = subprocess.check_output(['samtools', 'faidx', genome, '{}:{}-{}'.format(chrom, start, end)])
    return ''.join(output.split('\n')[1:]).upper()

def primerDesignRun(params, vrs):
    # Number of variants
    idcounter = len(vrs.keys())
    
    # Generate primer candidates
    primer3(params, vrs)

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

    # Clean-up
    os.remove(params['fnSilicaIn'])
    os.remove(params['fnSilicaOut1'])
                    
    # Iterate all primers
    scorelst = []
    for i in range(idcounter):
        for j in range(params['nprimer']):
            if (pritable[i * params['nprimer'] + j, 0] > 0) and (pritable[i * params['nprimer'] + j, 1] > 0):
                sumhits = pritable[i * params['nprimer'] + j, 0] + pritable[i * params['nprimer'] + j, 1]
                scorelst.append((sumhits, i, j))

    # Process primers in-order
    scorelst = sorted(scorelst)
    used = numpy.full(idcounter * params['nprimer'], False)
    vartodo = 0
    for idname in range(idcounter):
        if not vrs[idname]['done']:
            vartodo += 1
    while vartodo:
        print("Primer selection for", vartodo, "variants")
        varincl = numpy.full(idcounter, False)
        with open(params['fnSilicaIn'], 'w') as f:
            for (score, idname, candidate) in scorelst:
                if used[idname * params['nprimer'] + candidate]:
                    continue
                if (not vrs[idname]['done']) and (not varincl[idname]):
                    varincl[idname] = True
                    print(">" + str(idname) + "_PRIMER_LEFT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 0)], file=f)
                    print(">" + str(idname) + "_PRIMER_RIGHT_" + str(candidate) + "_SEQUENCE", file=f)
                    print(priseq[(idname * params['nprimer'] + candidate, 1)], file=f)
                    used[idname * params['nprimer'] + candidate] = True
                    if sum(varincl) == vartodo:
                        break

        # No more candidate primers
        if not sum(varincl):
            os.remove(params['fnSilicaIn'])
            return vartodo

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
                        if (vrs[idname]['type'].startswith("DEL")) or (vrs[idname]['type'].startswith("SNV")) or (vrs[idname]['type'].startswith("SNP")) or (vrs[idname]['type'].startswith("INS")):
                            if (d['Chrom'] == vrs[idname]['chr1']) and (d['ForPos'] < vrs[idname]['pos1']) and (d['RevPos'] > vrs[idname]['pos2']):
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
                        if lr == "LEFT":
                            if (vrs[idname]['type'].startswith("DEL")) or (vrs[idname]['type'].startswith("SNV")) or (vrs[idname]['type'].startswith("SNP")) or (vrs[idname]['type'].startswith("INS")) or (vrs[idname]['type'] == "BND_3to5"):
                                if (d['Chrom'] == vrs[idname]['chr1']) and (d['Pos'] < vrs[idname]['pos1']) and (d['Pos'] + params['pcrLen'] > vrs[idname]['pos1']) and (d['Ori'] == "forward"):
                                    prleftjs[idname] = d
                            elif (vrs[idname]['type'] == "INV_3to3") or (vrs[idname]['type'] == "BND_3to3"):
                                if (d['Chrom'] == vrs[idname]['chr1']) and (d['Pos'] < vrs[idname]['pos1']) and (d['Pos'] + params['pcrLen'] > vrs[idname]['pos1']) and (d['Ori'] == "forward"):
                                    prleftjs[idname] = d
                            elif (vrs[idname]['type'] == "INV_5to5") or (vrs[idname]['type'] == "BND_5to5"):
                                if (d['Chrom'] == vrs[idname]['chr1']) and (d['Pos'] > vrs[idname]['pos1']) and (d['Pos'] < vrs[idname]['pos1'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prleftjs[idname] = d
                            elif (vrs[idname]['type'].startswith("DUP")) or (vrs[idname]['type'] == "BND_5to3"):
                                if (d['Chrom'] == vrs[idname]['chr1']) and (d['Pos'] > vrs[idname]['pos1']) and (d['Pos'] < vrs[idname]['pos1'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prleftjs[idname] = d
                        if lr == "RIGHT":
                            if (vrs[idname]['type'].startswith("DEL")) or (vrs[idname]['type'].startswith("SNV")) or (vrs[idname]['type'].startswith("SNP")) or (vrs[idname]['type'].startswith("INS")) or (vrs[idname]['type'] == "BND_3to5"):
                                if (d['Chrom'] == vrs[idname]['chr2']) and (d['Pos'] > vrs[idname]['pos2']) and (d['Pos'] < vrs[idname]['pos2'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prrightjs[idname] = d
                            elif (vrs[idname]['type'] == "INV_3to3") or (vrs[idname]['type'] == "BND_3to3"):
                                if (d['Chrom'] == vrs[idname]['chr2']) and (d['Pos'] < vrs[idname]['pos2']) and (d['Pos'] + params['pcrLen'] > vrs[idname]['pos2']) and (d['Ori'] == "forward"):
                                    prrightjs[idname] = d
                            elif (vrs[idname]['type'] == "INV_5to5") or (vrs[idname]['type'] == "BND_5to5"):
                                if (d['Chrom'] == vrs[idname]['chr2']) and (d['Pos'] > vrs[idname]['pos2']) and (d['Pos'] < vrs[idname]['pos2'] + params['pcrLen']) and (d['Ori'] == "reverse"):
                                    prrightjs[idname] = d
                            elif (vrs[idname]['type'].startswith("DUP")) or (vrs[idname]['type'] == "BND_5to3"):
                                if (d['Chrom'] == vrs[idname]['chr2']) and (d['Pos'] < vrs[idname]['pos2']) and (d['Pos'] + params['pcrLen'] > vrs[idname]['pos2']) and (d['Ori'] == "forward"):
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
                            #print(prleftjs[idname]['Name'], prrightjs[idname]['Name'])
                            vrs[idname]['Primer1Chrom'] = prleftjs[idname]['Chrom']
                            vrs[idname]['Primer1Pos'] = prleftjs[idname]['Pos']
                            vrs[idname]['Primer1Seq'] = prleftjs[idname]['Seq']
                            vrs[idname]['Primer1Tm'] = prleftjs[idname]['Tm']
                            vrs[idname]['Primer1Ori'] = prleftjs[idname]['Ori']
                            vrs[idname]['Primer1Hits'] = pricount[pleft]
                            vrs[idname]['Primer2Chrom'] = prrightjs[idname]['Chrom']
                            vrs[idname]['Primer2Pos'] = prrightjs[idname]['Pos']
                            vrs[idname]['Primer2Seq'] = prrightjs[idname]['Seq']
                            vrs[idname]['Primer2Tm'] = prrightjs[idname]['Tm']
                            vrs[idname]['Primer2Ori'] = prrightjs[idname]['Ori']
                            vrs[idname]['Primer2Hits'] = pricount[pright]
                            vrs[idname]['done'] = True

        # Update ToDo
        vartodo = 0
        for idname in range(idcounter):
            if not vrs[idname]['done']:
                vartodo += 1

        # Clean-up
        os.remove(params['fnSilicaIn'])
        os.remove(params['fnSilicaOut1'])
        os.remove(params['fnSilicaOut2'])
    return vartodo



def primerDesign(filename, genome, prefix):
    # Parameters
    params = {'fnPrimer3In': prefix + ".primer3.input", 'fnPrimer3Out': prefix + ".primer3.output", 'fnSilicaIn': prefix + ".silica.input", 'fnSilicaOut1': prefix + ".silica.primer.output", 'fnSilicaOut2': prefix + ".silica.amplicon.output", 'genome': genome, 'pcrLen': 500, 'spacer': 50, 'nprimer': 100, 'PRIMER_MAX_TM': 63, 'PRIMER_MIN_TM': 57, 'PRIMER_OPT_TM': 60, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 20, 'PRIMER_MAX_SIZE': 25, 'maxmatch': 5}

    # Load variants
    variants = dict()
    idcounter = 0
    with open(filename, 'rb') as csvfile:
        prireader = csv.DictReader(csvfile, delimiter='\t')
        for row in prireader:
            if not row['chr1'].startswith("#"):
                variants[idcounter] = {'chr1': row['chr1'], 'pos1': int(row['pos1']), 'chr2': row['chr2'], 'pos2': int(row['pos2']), 'type': row['type'], 'done': False}
                idcounter += 1

    # Iterate parameters
    vartodo = 0
    for (pcrlen, spacer, nprimer, maxmatch) in [(250, 50, 25, 5), (500, 50, 100, 5), (1000, 50, 500, 50)]:
    #for (pcrlen, spacer, nprimer) in [(500, 50, 100)]:
        params['pcrLen'] = pcrlen
        params['spacer'] = spacer
        params['nprimer'] = nprimer
        params['maxmatch'] = maxmatch
        params['PRIMER_OPT_TM'] = 60
        vartodo = primerDesignRun(params, variants)
        if not vartodo:
            break
        print("Parameter iteration done. Remaining variants", vartodo)
    if vartodo != 0:
        print("Leftover variants no primers were found", vartodo)
        for idname in range(idcounter):
            if not variants[idname]['done']:
                print(variants[idname])
                variants[idname]['Primer1Chrom'] = 'n/a'
                variants[idname]['Primer1Pos'] = 0
                variants[idname]['Primer1Seq'] = 'n/a'
                variants[idname]['Primer1Tm'] = 0
                variants[idname]['Primer1Ori'] = 'n/a'
                variants[idname]['Primer1Hits'] = 0
                variants[idname]['Primer2Chrom'] = 'n/a'
                variants[idname]['Primer2Pos'] = 0
                variants[idname]['Primer2Seq'] = 'n/a'
                variants[idname]['Primer2Tm'] = 0
                variants[idname]['Primer2Ori'] = 'n/a'
                variants[idname]['Primer2Hits'] = 0
    return list(variants.values())


    
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
