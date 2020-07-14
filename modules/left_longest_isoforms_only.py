#!/usr/bin/env python3
"""Save bed12 file with longest isoforms only."""
import argparse
import sys
from collections import defaultdict


def parse_args():
    """Read and parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("bed_file")
    app.add_argument("isoforms_file")
    app.add_argument("output")
    app.add_argument("--alt", default=False, action="store_true", dest="alt",
                     help="Output gene: isoform, not bed12")
    if len(sys.argv) < 4:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def bed_to_gene(line):
    """Extract gene name."""
    return line.split("\t")[3]


def read_isoforms(isoforms_file):
    """Read isiforms file."""
    f = open(isoforms_file, "r")
    gene_to_isoforms = defaultdict(list)
    for line in f:
        line_data = line[:-1].split("\t")
        gene = line_data[0]
        iform = line_data[1]
        gene_to_isoforms[gene].append(iform)
    f.close()
    return gene_to_isoforms


def read_bed_lines(bed_file):
    """Read bed file."""
    iform_to_bed = {}
    f = open(bed_file, "r")
    for line in f:
        iform = bed_to_gene(line)
        iform_to_bed[iform] = line.rstrip()
    f.close()
    return iform_to_bed


def connect_gene_to_bed(gene_to_iforms, iforms_to_bed):
    """Create gene: bed lines."""
    gene_to_lines = {}
    for gene, iforms in gene_to_iforms.items():
        bed_lines = [iforms_to_bed.get(i) for i in iforms if iforms_to_bed.get(i)]
        if len(bed_lines) == 0:
            continue
        gene_to_lines[gene] = bed_lines
    return gene_to_lines


def get_bed_len(bed_line):
    """Return summ exons length."""
    block_sizes = [int(x) for x in bed_line.split("\t")[10].split(",") if x]
    return sum(block_sizes)

def get_longest(gene_to_lines):
    """Extract longest tracks."""
    longest = {}
    for gene, bed_lines in gene_to_lines.items():
        if len(bed_lines) == 1:
            # nothing to do
            longest[gene] = bed_lines[0]
            continue
        lens = [get_bed_len(b) for b in bed_lines]
        bed_lens = sorted(zip(bed_lines, lens), key=lambda x: x[1])
        longest_bed = bed_lens[-1][0]
        longest[gene] = longest_bed
    return longest


def main():
    """Entry point."""
    args = parse_args()
    gene_to_isoforms = read_isoforms(args.isoforms_file)
    iform_to_bed = read_bed_lines(args.bed_file)
    gene_to_bed = connect_gene_to_bed(gene_to_isoforms, iform_to_bed)
    longest_lines = get_longest(gene_to_bed)
    if not args.alt:
        f = open(args.output, "w") if args.output != "stdout" else sys.stdout
        f.write("\n".join(longest_lines) + "\n")
        f.close() if args.output != "stdout" else None
    else:
        f = open(args.output, "w") if args.output != "stdout" else sys.stdout
        for k, v in longest_lines.items():
            f.write(f"{k}\t{bed_to_gene(v)}\n")
        f.close() if args.output != "stdout" else None

if __name__ == "__main__":
    main()
