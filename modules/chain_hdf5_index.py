#!/usr/bin/env python3
"""Convert text bed file to hdf5.

TOGA uses HDF5 to index chain and bed files.
These files could be quite big (especially chain files).
Search for a particular chain of gene in the text file
might take too long. HDF5 helps to speed it up.
"""
import sys
import os
from collections import defaultdict
import numpy as np
import h5py

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def chain_hdf5_index(inchain, outdb):
    current_chain = None  # initiate current_chain variable
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

    if len(chain_data.keys()) == 0: 
        # if the chain_data dict is empty -> the entire chain file was empty
        # need to stop TOGA execution
        sys.stderr.write(f"(chain_hdf5_index.py) Error! Input file {inchain} is empty! Aborted.\n")
        sys.exit(1)

    # write chain data to bdb file
    h = h5py.File(outdb, "w")
    sys.stderr.write("Writing to HDF5...\n")
    for k, v in chain_data.items():
        # now value is a list of strings
        # need to make list of str -> str -> bytes conversion:
        chain_data = "".join(v)[:-1]
        h.create_dataset(k, data=np.string_(chain_data))
    h.close()


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    try:  # read args
        inchain = sys.argv[1]
        outdb = sys.argv[2]
    except IndexError:  # if args are wrong: show usage and quit
        sys.stderr.write("Usage: {0} [in_chain] [out_hdf5]\n".format(sys.argv[0]))
        sys.exit(0)
    chain_hdf5_index(inchain, outdb)
