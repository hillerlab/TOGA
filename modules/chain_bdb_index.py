#!/usr/bin/env python3
"""Convert text chain file to bdb.

TOGA uses Berkeley DB to index chain and bed files.
These files could be quite long (especially chain files).
Search for a particular chain of gene in the text file
might take too long. Berkeley DB helps to speed it up.
"""
import sys
from collections import defaultdict
import bsddb3

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def chain_bdb_index(inchain, outbdb):
    current_chain = 0  # initiate current_chain variable
    f = open(inchain, "r")  # read chain file
    # we will fill chain_data dictionary
    # structure is: chain_id: [corresponding lines]\
    # because one chain takes > 1 lines
    chain_data = defaultdict(list)
    sys.stderr.write("Making chain_id: chain_data dict...\n")
    for line in f:
        if line.startswith("#"):
            # this is a comment
            continue
        elif line.startswith("c"):
            # means that the first word is "chain"
            # this is a header line, contains chain_id
            chain_id = line[:-1].split()[-1]
            current_chain = chain_id
        elif current_chain == 0:
            # if so -> we didn't reach any chain header yet (init value)
            # so we read lines before the first chain in the file
            continue
        # add line to chains' list
        chain_data[current_chain].append(line)
    f.close()

    if len(chain_data.keys()) == 0:  # if the chain_data dict is empty -> the entire chain file was empty
        # need to stop TOGA execution
        sys.stderr.write("(chain_bdb_index.py) Error! Input file {0} is empty! Aborted.\n".format(inchain))
        sys.exit(1)

    # write chain data to bdb file
    db = bsddb3.btopen(outbdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in chain_data.items():
        db_key = k.encode()  # bdb requires bytestrings
        # now value is a list of strings
        # need to make list of str -> str -> bytes conversion:
        chain_data = "".join(v)[:-1].encode()
        db[db_key] = chain_data
    db.close()  # we are done


if __name__ == "__main__":
    try:  # read args
        inchain = sys.argv[1]
        outbdb = sys.argv[2]
    except IndexError:  # if args are wrong: show usage and quit
        sys.stderr.write("Usage: {0} [in_chain] [out_bdb]\n".format(sys.argv[0]))
        sys.exit(0)
    chain_bdb_index(inchain, outbdb)
