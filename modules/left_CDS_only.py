#!/usr/bin/env python3
"""Remove UTRs from a bed-12 file.

Use of for pre-processing of your bed-12 file.
"""
import argparse
import sys
from collections import Counter

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def die(msg, rc=1):
    """Show msg in stderr, exit with the rc given."""
    sys.stderr.write(msg + "\n")
    sys.exit(rc)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input", help="Bed-12 formatted annotation track.")
    app.add_argument("output", default="stdout", help="Output destination, stdout as default")
    app.add_argument("--out_of_frame", "--ouf", action="store_true", dest="out_of_frame",
                     help="Do not skip out-of-frame genes.")
    app.add_argument("--save_rejected", "--sr", default=None)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def prepare_bed_file(bed_file, output, utred_out=False, ouf=False, save_rejected=None):
    """Return gene: exons dict."""
    new_lines = []
    utred_lines = []
    rejected = []
    names = Counter()  # we need to make sure that all names are unique

    f = open(bed_file, "r")
    for line in f:
        line_data = line[:-1].split("\t")

        if len(line_data) != 12:
            f.close()
            die(f"Error! Bed 12 file is required! Got a file with {len(line_data)} fields instead")

        chrom = line_data[0]
        chromStart = int(line_data[1])
        chromEnd = int(line_data[2])
        name = line_data[3]  # gene_name usually
        bed_score = int(line_data[4])  # never used
        strand = line_data[5]  # otherwise:
        # strand = True if line_data[5] == '+' else False
        thickStart = int(line_data[6])
        thickEnd = int(line_data[7])
        itemRgb = line_data[8]  # never used
        blockCount = int(line_data[9])
        blockSizes = [int(x) for x in line_data[10].split(',') if x != '']
        blockStarts = [int(x) for x in line_data[11].split(',') if x != '']
        blockEnds = [blockStarts[i] + blockSizes[i] for i in range(blockCount)]
        blockAbsStarts = [blockStarts[i] + chromStart for i in range(blockCount)]
        blockAbsEnds = [blockEnds[i] + chromStart for i in range(blockCount)]
        blockNewStarts, blockNewEnds = [], []
        names[name] += 1
        if names[name] > 1:
            name_upd = f"{name}_{names[name]}"
            rejected.append((name, f"Non uniq, renamed to {name_upd}"))
        else:
            name_upd = name

        if thickStart > thickEnd:
            f.close()
            die("Error! File is corrupted, thickEnd MUST be >= thickStart")

        for block_num in range(blockCount):
            blockStart = blockAbsStarts[block_num]
            blockEnd = blockAbsEnds[block_num]

            # skip the block if it is entirely UTR
            if blockEnd <= thickStart:
                continue
            elif blockStart >= thickEnd:
                continue

            # remove UTRs
            blockNewStart = blockStart if blockStart >= thickStart else thickStart
            blockNewEnd = blockEnd if blockEnd <= thickEnd else thickEnd
            blockNewStarts.append(blockNewStart - thickStart)
            blockNewEnds.append(blockNewEnd - thickStart)

        if len(blockNewStarts) == 0:
            rejected.append((name, "No CDS"))
            continue  # remove non-coding genes

        blockNewCount = len(blockNewStarts)
        blockNewSizes = [blockNewEnds[i] - blockNewStarts[i]
                         for i in range(blockNewCount)]

        if sum(blockNewSizes) % 3 != 0 and not ouf:
            rejected.append((name, "Out-of-frame gene"))
            continue  # out-of-frame gene

        new_track = [chrom, thickStart, thickEnd, name_upd, bed_score,
                     strand, thickStart, thickEnd, itemRgb, blockNewCount,
                     ",".join([str(x) for x in blockNewSizes]) + ",",
                     ",".join([str(x) for x in blockNewStarts]) + ","]
        new_line = "\t".join([str(x) for x in new_track])
        new_lines.append(new_line)

        # save UTRed track
        utr_track = line_data[:]
        utr_track[6] = utr_track[1]
        utr_track[7] = utr_track[2]
        utr_line = "\t".join(utr_track)
        utred_lines.append(utr_line)
    f.close()

    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write("\n".join(new_lines) + "\n")
    f.close() if output != "stdout" else None

    if save_rejected:
        f = open(save_rejected, "w") if save_rejected != "stdout" else sys.stdout
        for elem in rejected:
            f.write(f"{elem[0]}\t{elem[1]}\n")
        f.close() if save_rejected != "stdout" else None
    
    if utred_out:
        f = open(utred_out, "w")
        f.write("\n".join(utred_lines) + "\n")
        f.close()


def main():
    """Entry point."""
    args = parse_args()
    # maybe save rejected also here
    prepare_bed_file(args.input, args.output, ouf=args.out_of_frame, save_rejected=args.save_rejected)
    sys.exit(0)


if __name__ == "__main__":
    main()
