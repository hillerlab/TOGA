#!/usr/bin/env python3
"""Convert text chain file to bdb."""
import sys
from collections import defaultdict
import bsddb3

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def chain_bdb_index(inchain, outbdb):
    current_chain = 0
    chain_data = defaultdict(list)
    f = open(inchain, "r")
    sys.stderr.write("Making chain_id: chain_data dict...\n")
    for line in f:
        if line.startswith("#"):
            continue
        elif line.startswith("c"):
            chain_id = line[:-1].split()[-1]
            current_chain = chain_id
            # sys.stderr.write("Chain: {0}\r".format(current_chain))
        elif current_chain == 0:
            continue
        chain_data[current_chain].append(line)
    f.close()

    if len(chain_data.keys()) == 0:
        sys.stderr.write("(chain_bdb_index.py) Error! Input file {0} is empty! Aborted.\n".format(inchain))
        sys.exit(1)

    db = bsddb3.btopen(outbdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in chain_data.items():
        db_key = k.encode()
        chain_data = "".join(v)[:-1].encode()
        db[db_key] = chain_data
    db.close()


if __name__ == "__main__":
    try:
        inchain = sys.argv[1]
        outbdb = sys.argv[2]
    except IndexError:
        sys.stderr.write("Usage: {0} [in_chain] [out_bdb]\n".format(sys.argv[0]))
        sys.exit(0)
    chain_bdb_index(inchain, outbdb)
