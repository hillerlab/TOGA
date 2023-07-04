#!/usr/bin/env python3
"""Convert text bed file to hdf5.

This allows TOGA to extract bed track
for a particular transcript ID immediately.
"""
import sys
import os
import numpy as np
import h5py
from version import __version__

__author__ = "Bogdan Kirilenko, 2020."
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def bed_hdf5_index(in_bed, out_db):
    # read the bed file
    f = open(in_bed, "r")  # assume each bed line has unique name (field 3)
    h = h5py.File(out_db, "w")
    lines_counter = 0
    sys.stderr.write("Making gene_id: bed_data dict...\n")
    for line in f:
        gene_id = line.split("\t")[3]
        h.create_dataset(gene_id, data=np.string_(line))
        lines_counter += 1
    f.close()
    h.close()

    if lines_counter == 0:  # meaning bed file was empty
        # this should not happen: halt TOGA
        sys.stderr.write(
            f"(bed_hdf5_index.py) Error! Input file {in_bed} is empty! Aborted.\n"
        )
        sys.exit(1)


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could fail
    try:  # read arguments
        in_bed_arg = sys.argv[1]
        out_db_arg = sys.argv[2]
    except IndexError:  # not enough arguments: show usage message and quit
        sys.stderr.write("Usage: {0} [in_bed] [out_hdf5]\n".format(sys.argv[0]))
        sys.exit(0)
    bed_hdf5_index(in_bed_arg, out_db_arg)
