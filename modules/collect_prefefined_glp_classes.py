#!/usr/bin/env python3
"""Collect transcript/projection classifications that are already known.

For example, a transcript that was not intersected by any chain is classified as
Missing, or if projection locus is very short (<5%CDS) it could be missing or lost.
"""
import os
from version import __version__

PROJECTION = "PROJECTION"
TRANSCRIPT = "TRANSCRIPT"


def _add_transcripts_to_missing(*lists):
    """Generate a list of missing entries."""
    ret = []
    for lst in lists:
        for trans_ in lst:
            item = (TRANSCRIPT, trans_, "M")
            ret.append(item)
    return ret


def _collect_predefined_glp_cases(glp_class_file):
    """Read table produced by split exon realing jobs script."""
    ret = []
    if not os.path.isfile(glp_class_file):
        return ret
    f = open(glp_class_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        entry_class = line_data[1]
        entry = line_data[0]
        classification = line_data[2]
        item = (entry_class, entry, classification)
        ret.append(item)
    f.close()
    return ret
