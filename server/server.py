#! /usr/bin/env python

from __future__ import print_function
from flask import Flask, jsonify, request
from flask_cors import CORS
from verdin import primerDesign
import os
import uuid
import datetime

VERDINWS = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
CORS(app)
app.config['VERDIN'] = os.path.join(VERDINWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['VERDIN'], "data")


def generate(inputjs):
    # Todo parse genome from UI
    genome = os.path.join(app.config['VERDIN'], "fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz")

    # Generate UUID
    uuidstr = str(uuid.uuid4())
    sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
    if not os.path.exists(sf):
        os.makedirs(sf)

    # Generate input file
    infile = os.path.join(sf, "verdin_" + uuidstr + ".variants")
    with open(infile, "w") as infi:
        print("chr1", "pos1", "chr2", "pos2", "type", sep='\t', file=infi)
        for v in inputjs:
            print(v['chr1'], v['pos1'], v['chr2'], v['pos2'], v['type'], sep='\t', file=infi)
        infi.close()

    # Run primer design
    prefix = os.path.join(sf, "verdin_" + uuidstr)
    primerlst = primerDesign(infile, genome, prefix)
    return jsonify(primerlst)
    

@app.route('/primers', methods=['POST'])
def generatePost():
    genome = os.path.join(app.config['VERDIN'], "fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz")
    variants = request.get_json()
    return generate(variants)
    
@app.route('/primers', methods=['GET'])
def generateGet():
    #genome = request.args.get("build")
    genome = os.path.join(app.config['VERDIN'], "fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz")
    variants = [{"chr1": request.args.get("chr1"), "pos1": request.args.get("pos1"), "chr2": request.args.get("chr2"), "pos2": request.args.get("pos2"), "type": request.args.get("type")}]
    return generate(variants)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=3300, debug=True, threaded=True)
