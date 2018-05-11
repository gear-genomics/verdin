#! /usr/bin/env python

from __future__ import print_function
from flask import Flask, jsonify, request
from flask_cors import CORS
import os
import uuid
import datetime

BGENDIR = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
CORS(app)

@app.route('/primers', methods=['GET'])
def generate():
    chrom1 = request.args.get("chr1")
    return jsonify(chrom1=chrom1)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=3300, debug=True, threaded=True)
