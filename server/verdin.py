#! /usr/bin/env python

from __future__ import print_function
from flask import Flask, jsonify, request
from flask_cors import CORS
import os
import uuid
import datetime

VERDINWS = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
CORS(app)
app.config['VERDIN'] = os.path.join(VERDINWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['VERDIN'], "data")

@app.route('/primers', methods=['GET'])
def generate():
    genome = request.args.get("build")
    chr1 = request.args.get("chr1")
    pos1 = request.args.get("pos1")
    chr2 = request.args.get("chr2")
    pos2 = request.args.get("pos2")
    svt = request.args.get("svtype")

    # Generate UUID
    uuidstr = str(uuid.uuid4())
    sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
    if not os.path.exists(sf):
        os.makedirs(sf)

    # Generate input file
    infile = os.path.join(sf, "verdin_" + uuidstr + ".variants")
    with open(infile, "w") as infi:
        print(chr1, pos1, chr2, pos2, svt, file=infi)
        infi.close()

    return jsonify(chrom1=chr1)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=3300, debug=True, threaded=True)
