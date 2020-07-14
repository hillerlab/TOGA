#!/usr/bin/env python3
"""Convert text chain file to bdb."""
import sys
import bsddb3

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def bed_bdb_index(inbed, outbdb):
    gene_data = {}
    f = open(inbed, "r")
    sys.stderr.write("Making gene_id: bed_data dict...\n")
    for line in f:
        gene_id = line.split("\t")[3]
        gene_data[gene_id] = line
    f.close()

    if len(gene_data.keys()) == 0:
        sys.stderr.write("(bed_bdb_index.py) Error! Input file {0} is empty! Aborted.\n".format(inbed))
        sys.exit(1)

    db = bsddb3.btopen(outbdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in gene_data.items():
        db[k.encode()] = v.encode()
    db.close()


if __name__ == "__main__":
    try:
        inbed = sys.argv[1]
        outbdb = sys.argv[2]
    except IndexError:
        sys.stderr.write("Usage: {0} [in_bed] [out_bdb]\n".format(sys.argv[0]))
        sys.exit(0)
    bed_bdb_index(inbed, outbdb)
