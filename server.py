#! /usr/bin/env python

from __future__ import print_function
from flask import Flask, jsonify
import os
import uuid
import datetime

app = Flask(__name__, static_url_path='/static')

BGENDIR = os.path.dirname(os.path.abspath(__file__))

@app.route('/<chr1>/<int:pos1>/<chr2>/<int:pos2>/<svtype>', methods=['GET'])
def generate(chr1, pos1, chr2, pos2, svtype):
    return jsonify(chrom1=chr1, ps1=pos1, chrom2=chr2, ps2=pos2, svt=svtype)

@app.route('/')
def root():
    return app.send_static_file('index.html')


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=3300, debug=True, threaded=True)
