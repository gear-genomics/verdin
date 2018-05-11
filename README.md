Variant primer design that checks in-silico PCR primers across complete genomes.

Installing
----------

`sudo apt-get install -y build-essential g++ cmake git-all liblzma-dev zlib1g-dev libbz2-dev liblzma-dev`

`git clone --recursive https://github.com/gear-genomics/verdin.git`

`cd verdin`

`make all`


Install Genome Index
--------------------

Currently done by Indigo, place in fm folder.


Running
-------

`python server/verdin.py -v variants.tsv -g genome.fa.gz -p outprefix`
