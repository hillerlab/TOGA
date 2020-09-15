#!/usr/bin/env python3
"""Make togaPlot.tab."""
import argparse
import os
import sys
sys.path.append('../supply')
sys.path.append('supply')
from collections import defaultdict
try:
    from ucsc_browser_visualisation.generate_tab_files import split_proj_name
    from supply.plot_mutations import make_plot
except ImportError:
    from generate_tab_files import split_proj_name
    from plot_mutations import make_plot


def parse_args():
    """Read and parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("wd")
    app.add_argument("output_tab")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def get_projections(query_annot):
    """Extract tuples (transcript id, chain id) from query bed."""
    proj_list = []
    f = open(query_annot, "r")
    for line in f:
        projection_str = line.split("\t")[3]
        trans, chain = split_proj_name(projection_str)
        proj_tup = (trans, chain, )
        proj_list.append(proj_tup)
    f.close()
    return proj_list


def make_inact_dict(mut_file):
    """Make projection: mut lines dict."""
    proj_to_lines = defaultdict(list)
    f = open(mut_file, "r")
    for line in f:
        if not line.startswith("#"):
            continue
        line_data = line.rstrip().split("\t")
        trans = line_data[0].replace("# ", "")
        chain = line_data[1]
        projection = f"{trans}.{chain}"
        proj_to_lines[projection].append(line)
    f.close()
    return proj_to_lines


def main():
    """Entry point."""
    args = parse_args()
    query_annotation = os.path.join(args.wd, "query_annotation.bed")
    reference_annotation = os.path.join(args.wd, "toga_filt_ref_annot.bed")
    mut_file = os.path.join(args.wd, "inact_mut_data.txt")
    proj_to_lines = make_inact_dict(mut_file)
    projections = get_projections(query_annotation)
    projections_num = len(projections)
    print(f"Extracted {projections_num} projections")
    f = open(args.output_tab, "w")
    for num, (trans, chain) in enumerate(projections, 1):
        proj_key = f"{trans}.{chain}"
        inact_lines = proj_to_lines[proj_key]
        svg_line = make_plot(reference_annotation, mut_file, trans,
                             chain, None, None, None, inact_lines)
        print(f"{num}/{projections_num}", end="\r")
        svg_fmt = svg_line.replace("\n", "").replace("\t", " ")
        f.write(f"{proj_key}\t{svg_fmt}\n")
    f.close()


if __name__ == "__main__":
    main()
