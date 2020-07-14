#!/usr/bin/env python3
"""Create a bed file + a fasta file."""
import sys
import argparse
import os
try:  # need to add this nasty construction
      # depending on from where this script was called
      # TODO: just move parse cesar bdb content here
    from modules.parse_cesar_bdb import parse_cesar_bdb
except ImportError:
    from parse_cesar_bdb import parse_cesar_bdb

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def die(msg, rc=1):
    """Show msg in stderr, exit with the rc given."""
    sys.stderr.write(msg + "\n")
    sys.stderr.write("Program finished with exit code {}\n".format(rc))
    sys.exit(rc)


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input_dir", help="Directory containing output BDB files")
    app.add_argument("output_bed", help="Save pre_final bed12 file to...")
    app.add_argument("output_fasta", help="Save fasta fasta to...")
    app.add_argument("meta_data", help="Save exons metadata to...")
    app.add_argument("prot_fasta", help="Save protein fasta to...")
    app.add_argument("skipped", help="Save skipped genes")
    app.add_argument("--output_trash", default=None)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    # print help if there are no args
    return args


def merge_cesar_output(input_dir, output_bed, output_fasta,
                       meta_data_arg, skipped_arg, prot_arg,
                       output_trash):
    """Merge multiple CESAR output files."""
    die(f"Error! {input_dir} is not a dir!") \
        if not os.path.isdir(input_dir) else None
    bdbs = [x for x in os.listdir(input_dir) if x.endswith(".bdb")]

    bed_summary = []
    fasta_summary = []
    trash_summary = []
    meta_summary = []
    prot_summary = []
    skipped = []

    task_size = len(bdbs)
    # extract data for all the files
    for num, bdb_file in enumerate(bdbs):
        bdb_path = os.path.join(input_dir, bdb_file)
        try:
            parsed_data = parse_cesar_bdb(bdb_path)
        except AssertionError:
            sys.exit(f"Error! Failed reading file {bdb_file}")

        bed_lines = parsed_data[0]
        trash_exons = parsed_data[1]
        fasta_lines = parsed_data[2]
        meta_data = parsed_data[3]
        prot_fasta = parsed_data[4]
        skip = parsed_data[5]

        if len(bed_lines) == 0:
            eprint(f"Warning! {bdb_file} is empty")
            continue  # it is empty
        bed_summary.append("\n".join(bed_lines) + "\n")
        fasta_summary.append(fasta_lines)

        # exons_left_summary.append(exons_left)
        trash_summary.append("".join(trash_exons))
        meta_summary.append(meta_data)
        skipped.append(skip)
        prot_summary.append(prot_fasta)
        eprint(f"Reading file {num + 1}/{task_size}", end="\r")

    # save output
    eprint("Saving the output")

    if len(bed_summary) == 0:
        # if so, no need to continue
        eprint("! merge_cesar_output.py:")
        die("No projections found! Abort.")

    with open(output_bed, "w") as f:
        f.write("".join(bed_summary))
    with open(output_fasta, "w") as f:
        f.write("".join(fasta_summary))
    with open(meta_data_arg, "w") as f:
        f.write("\n".join(meta_summary))
    with open(skipped_arg, "w") as f:
        f.write("\n".join(skipped))
    with open(prot_arg, "w") as f:
        f.write("\n".join(prot_summary))

    if output_trash:
        f = open(output_trash, "w")
        f.write("".join(trash_summary))
        f.close()


def main():
    """Entry point."""
    args = parse_args()
    merge_cesar_output(args.input_dir,
                       args.output_bed,
                       args.output_fasta,
                       args.meta_data,
                       args.skipped,
                       args.prot_fasta,
                       args.output_trash)

if __name__ == "__main__":
    main()
