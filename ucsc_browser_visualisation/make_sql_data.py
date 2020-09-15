#!/usr/bin/env python3
"""Perform all operations to create SQL tables for UCSC browser."""
import os
import sys
import argparse
import subprocess

GENERATE_TAB_FILES = "generate_tab_files.py"
GENERATE_PLOTS = "make_togaPlot.py"
LOCATION = os.path.dirname(__file__)


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("project_dir", help="Directory containing TOGA output")
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args

def call_gen_tab_files(project_dir):
    """Call Generate tab files subprocess."""
    print(f"Calling {GENERATE_TAB_FILES}")
    out_dir = os.path.join(project_dir, "tabs")
    os.mkdir(out_dir) if not os.path.isdir(out_dir) else None
    exe_ = os.path.join(LOCATION, GENERATE_TAB_FILES)
    cmd = f"{exe_} {project_dir} {out_dir}"
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        raise ValueError(f"Command {cmd} died - cannot generate tab files")
    print("5 of 6 tab files generated")


def call_gen_plot_files(project_dir):
    """Generate svg plots."""
    print(f"Calling {GENERATE_PLOTS}")
    out_file = os.path.join(project_dir, "tabs", "togaPlot.tab")
    exe_ = os.path.join(LOCATION, GENERATE_PLOTS)
    cmd = f"{exe_} {project_dir} {out_file}"
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        raise ValueError(f"Command {cmd} died - cannot generate plots")
    print("All tab files generated")


def read_tab_file(tab_file):
    """Make trans_id: data track."""
    trans_to_dat = {}
    f = open(tab_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        data = line_data[1:]
        trans_to_dat[trans] = data
    f.close()
    return trans_to_dat


def make_toga_data(project_dir):
    """Merge togaPlot, togaProt, togaInfo and togaInactFeat tables into togaData."""
    tabs_dir = os.path.join(project_dir, "tabs")
    toga_info_file = os.path.join(tabs_dir, "togaInfo.tab")
    toga_inact_feat_file = os.path.join(tabs_dir, "togaInactFeat.tab")
    toga_plot_file = os.path.join(tabs_dir, "togaPlot.tab")
    toga_prot_file = os.path.join(tabs_dir, "togaProt.tab")
    out_file = os.path.join(tabs_dir, "togaData.tab")

    print("Reading tab files to merge")
    trans_to_info = read_tab_file(toga_info_file)
    trans_to_ifeat = read_tab_file(toga_inact_feat_file)
    trans_to_plot = read_tab_file(toga_plot_file)
    trans_to_prot = read_tab_file(toga_prot_file)
    
    transcripts = trans_to_info.keys()
    print(f"Got data for {len(transcripts)} transcripts")

    f = open(out_file, "w")
    for t in transcripts:
        info = trans_to_info[t]
        ifeat = trans_to_ifeat[t]
        prot = trans_to_prot[t]
        plot = trans_to_plot[t]
        data_lst = [t] + info + ifeat + prot + plot
        tab_str = "\t".join(data_lst)
        f.write(tab_str)
        f.write("\n")
    f.close()
    print("Done")


def main():
    """Entry point."""
    args = parse_args()
    # call generate_tab_files.py and make_togaPlot.py
    # call_gen_tab_files(args.project_dir)
    # call_gen_plot_files(args.project_dir)
    # merge 4 transcriptID-related tables into a single one
    make_toga_data(args.project_dir)


if __name__ == "__main__":
    main()
