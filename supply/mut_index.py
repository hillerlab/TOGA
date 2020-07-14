#!/usr/bin/env python3
"""Make mutations index file."""
import sys
from collections import defaultdict
import bsddb3

MUT_LINE_FIELDS = 8


def make_mut_index(in_txt, out_bdb):
    """Make mut index: transcipt: mutations."""
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

    db = bsddb3.btopen(out_bdb, "w")
    sys.stderr.write("Writing to BDB...\n")
    for k, v in trans_to_text.items():
        db[k.encode()] = v.encode()
    db.close()


if __name__ == "__main__":
    try:
        in_txt = sys.argv[1]
        out_bdb = sys.argv[2]
    except IndexError:
        sys.exit(f"Usage: {sys.argv[0]} [in_muts] [out_bdb]")
    make_mut_index(in_txt, out_bdb)
