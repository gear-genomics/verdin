#! /usr/bin/env python

from __future__ import print_function
import sys
import argparse
import json

def jsonToTable(filename):
    print('chr1', 'pos1', 'chr2', 'pos2', 'type', 'Primer1Chrom', 'Primer1Pos', 'Primer1Seq', 'Primer1Ori', 'Primer1Hits', 'Primer1Tm', 'Primer2Chrom', 'Primer2Pos', 'Primer2Seq', 'Primer2Ori', 'Primer2Hits', 'Primer2Tm', sep="\t")
    with open(filename) as prijson:
        d = json.load(prijson)
        for r in d:
            print(r['chr1'], r['pos1'], r['chr2'], r['pos2'], r['type'], r['Primer1Chrom'], r['Primer1Pos'], r['Primer1Seq'], r['Primer1Ori'], r['Primer1Hits'], r['Primer1Tm'], r['Primer2Chrom'], r['Primer2Pos'], r['Primer2Seq'], r['Primer2Ori'], r['Primer2Hits'], r['Primer2Tm'], sep="\t")


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verdin')
    parser.add_argument('-v', '--verdin', required=True, metavar="verdin.json", dest='verdin', help='verdin json output')
    args = parser.parse_args()

    if args.verdin:
        jsonToTable(args.verdin)

