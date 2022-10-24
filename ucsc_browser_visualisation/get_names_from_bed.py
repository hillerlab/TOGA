#!/usr/bin/env python3
"""Just extract names from toga output bed file.

Works like xenoRefGenelx.pl"""
import sys

if len(sys.argv) != 2:
    to_read = None
    sys.exit(f"Usage: {sys.argv[0]} [query_annotation.bed] | sort -u > ix.txt")
else:
    to_read = sys.argv[1]

f = open(to_read, "r")
for line in f:
    id_field = line.rstrip().split("\t")[3]
    no_chain_id = ".".join(id_field.split(".")[:-1])
    dot_split = no_chain_id.split(".")
    if len(dot_split) > 1:
        to_out = [id_field, no_chain_id] + dot_split
    else:
        to_out = [id_field, no_chain_id]
    line = "\t".join(to_out)
    print(line)
f.close()
