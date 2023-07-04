#!/usr/bin/env python3
"""Make mutations index file."""
import sys
import os
from collections import defaultdict
import numpy as np
import h5py
from version import __version__

MUT_LINE_FIELDS = 8


def make_mut_index(in_txt, out_bdb):
    """Make mut index: transcript: mutations."""
    transcript_to_lines = defaultdict(list)
    f = open(in_txt, "r")
    for line in f:
        if not line.startswith("#"):
            continue
        line_data = line.rstrip().split("\t")
        trans = line_data[0].replace("# ", "")
        transcript_to_lines[trans].append(line)
    f.close()
    trans_to_text = {}
    for k, v in transcript_to_lines.items():
        lines_merge = "".join(v)
        trans_to_text[k] = lines_merge

    h = h5py.File(out_bdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in trans_to_text.items():
        h.create_dataset(k, data=np.string_(v))
    h.close()


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    try:
        in_txt_arg = sys.argv[1]
        out_db_arg = sys.argv[2]
    except IndexError:
        sys.exit(f"Usage: {sys.argv[0]} [in_muts] [out_bdb]")
    make_mut_index(in_txt_arg, out_db_arg)
