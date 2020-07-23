#!/usr/bin/env python3
"""Convert text bed file to bdb.

TOGA uses Berkeley DB to index chain and bed files.
These files could be quite long (especially chain files).
Search for a particular chain of gene in the text file
might take too long. Berkeley DB helps to speed it up.
"""
import sys
import bsddb3

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def bed_bdb_index(inbed, outbdb):
    # read the bed file
    gene_data = {}  # fill this dict with transcript_ID: corresponding line
    f = open(inbed, "r")  # assume each bed line has unique name (field 3)
    sys.stderr.write("Making gene_id: bed_data dict...\n")
    for line in f:
        gene_id = line.split("\t")[3]
        gene_data[gene_id] = line
    f.close()

    if len(gene_data.keys()) == 0:  # gene_data dict is empty -> meaning bed file was empty
        # this should not happen: halt TOGA
        sys.stderr.write(f"(bed_bdb_index.py) Error! Input file {inbed} is empty! Aborted.\n")
        sys.exit(1)

    # save data to Berkeley DB
    db = bsddb3.btopen(outbdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in gene_data.items():
        # need bytestring: so .encode() added to strings
        db[k.encode()] = v.encode()
    db.close()  # we're done


if __name__ == "__main__":
    try:  # read arguments
        inbed = sys.argv[1]
        outbdb = sys.argv[2]
    except IndexError:  # not enough arguments: show usage message and quit
        sys.stderr.write("Usage: {0} [in_bed] [out_bdb]\n".format(sys.argv[0]))
        sys.exit(0)
    bed_bdb_index(inbed, outbdb)
