#!/usr/bin/env python3
"""Filter bed-12 file.

Remove:
- incomplete annotations
- genes without CDS
"""
import argparse
import sys
import re
from collections import Counter
from version import __version__

try:
    from modules.common import die
    from modules.common import eprint
    from modules.common import to_log
except ImportError:
    from common import die
    from common import eprint
    from common import to_log

__author__ = "Bogdan Kirilenko, 2023."

ALLOWED_CHARSET = "a-zA-Z0-9._-"
ALLOWED_CHARSET_RE = rf"[^{ALLOWED_CHARSET}]"
MODULE_NAME_FOR_LOG = "filter_bed"

def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input", help="Bed-12 formatted annotation track.")
    app.add_argument(
        "output", default="stdout", help="Output destination, stdout as default"
    )
    app.add_argument(
        "--out_of_frame",
        action="store_true",
        dest="out_of_frame",
        help="Do not skip out-of-frame genes.",
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def prepare_bed_file(bed_file, output, ouf=False, save_rejected=None, only_chrom=None):
    """Filter the bed file given and save the updated version."""
    new_lines = []  # keep updated lines
    rejected = []  # keep IDs of skipped transcripts + the reason why
    names = Counter()  # we need to make sure that all names are unique
    allowed_re = re.compile(ALLOWED_CHARSET_RE).search
    broken_names = []

    f = open(bed_file, "r")
    for num, line in enumerate(f, 1):
        # parse bed file according to specification
        line_data = line.rstrip().split("\t")

        if len(line_data) != 12:
            f.close()  # this is for sure an error
            # it is possible only if something except a bed12 was provided
            die(
                "Error! Bed 12 file is required! Got a file with {len(line_data)} fields instead"
            )

        chrom = line_data[0]
        if only_chrom and chrom != only_chrom:
            # TOGA allows to perform the analysis on a specific chromosome only
            # is so, we can skip all transcripts that located on other chromosomes
            continue
        chrom_start = int(line_data[1])
        chrom_end = int(line_data[2])
        name = line_data[3]  # gene_name usually
        corr_name = not bool(allowed_re(name))
        if corr_name is False:
            broken_names.append(name)
        # TODO: check weird characters in the transcript name
        # bed_score = int(line_data[4])  # never used
        # strand = line_data[5]  # otherwise:
        # strand = True if line_data[5] == '+' else False
        thick_start = int(line_data[6])
        thick_end = int(line_data[7])
        # itemRgb = line_data[8]  # never used
        block_count = int(line_data[9])
        block_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
        block_starts = [int(x) for x in line_data[11].split(",") if x != ""]
        block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
        block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
        block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
        block_new_starts, block_new_ends = [], []
        names[name] += 1

        if thick_start > thick_end:
            f.close()  # according to bed12 specification this should never happen
            err_msg = (
                f"{MODULE_NAME_FOR_LOG}: Error! Bed file is corrupted, thickEnd MUST be >= thickStart.\n"
                f"Problem occurred at line {num}, gene {name}"
            )
            to_log(err_msg)
            die(err_msg)
        elif thick_start == thick_end:
            # this means that this is a non-coding transcript
            # TOGA cannot process them: we can skip it
            rejected.append((name, "No CDS"))
            warn_msg = f"{MODULE_NAME_FOR_LOG}: transcript {name} skipped: no CDS"
            to_log(warn_msg)
            continue

        if thick_start < chrom_start or thick_end > chrom_end:
            # a very strange (but still possible) case
            f.close()  # for sure an error with input data
            err_msg = (
                f"{MODULE_NAME_FOR_LOG}: Error! Bed file is corrupted, thickRange is outside chromRange!\n"
                f"Problem occurred at line {num}, gene {name}"
            )
            to_log(err_msg)
            die(err_msg)

        # now select CDS only
        # we keep UTRs in the filtered file
        # however, we need CDS to check whether it's correct (% 3 == 0)
        for block_num in range(block_count):
            blockStart = block_abs_starts[block_num]
            blockEnd = block_abs_ends[block_num]

            # skip the block if it is entirely UTR
            if blockEnd <= thick_start:
                continue
            elif blockStart >= thick_end:
                continue

            # if we are here: this is not an entirely UTR exon
            # it might intersect the CDS border or to be in the CDS entirely
            # remove UTRs: block start must be >= CDS_start (thickStart)
            # block end must be <= CDS_end (thickEnd)
            block_new_start = blockStart if blockStart >= thick_start else thick_start
            block_new_end = blockEnd if blockEnd <= thick_end else thick_end
            block_new_starts.append(block_new_start - thick_start)
            block_new_ends.append(block_new_end - thick_start)

        if len(block_new_starts) == 0:
            # even it thickStart != thickEnd this transcript still can be non-coding
            # but if there are no blocks in the CDS -> we can catch this
            warn_msg = f"{MODULE_NAME_FOR_LOG}: transcript {name} skipped: no CDS"
            to_log(warn_msg)
            continue

        block_new_count = len(block_new_starts)
        block_new_sizes = [
            block_new_ends[i] - block_new_starts[i] for i in range(block_new_count)
        ]

        if sum(block_new_sizes) % 3 != 0 and not ouf:
            # this is an out-of-frame (or incomplete transcript)
            # ideally CDS length should be divisible by 3
            # not ouf means that we like to keep such transcripts for some reason
            rejected.append((name, "Out-of-frame gene"))
            warn_msg = f"{MODULE_NAME_FOR_LOG}: transcript {name} skipped: out-of-frame transcript"
            to_log(warn_msg)
            continue

        # we keep this transcript: add in to the list
        new_line = "\t".join([str(x) for x in line_data])
        new_lines.append(new_line)
    f.close()

    # if not allowed characters in transcript names: list them
    if len(broken_names) > 0:
        to_log(f"{MODULE_NAME_FOR_LOG}: Error! {len(broken_names)} transcript names contain not allowed characters")
        for t in broken_names:
            to_log(f"* {t}\n")
        err_msg = f"{MODULE_NAME_FOR_LOG}: Allowed characters are: {ALLOWED_CHARSET}"
        to_log(err_msg)
        die(err_msg)
    # if there are non-unique transcript IDs: die
    # I kill it there, not earlier to show them altogether
    if any(v > 1 for v in names.values()):
        err_msg = f"{MODULE_NAME_FOR_LOG}: Error! There are non-uniq transcript IDs:"
        to_log(err_msg)

        duplicates = [k for k, v in names.items() if v > 1]
        for d in duplicates:
            to_log(f"* {d}")
        die("")

    if len(new_lines) == 0:
        # no transcripts pass the filter: probably an input data mistake
        err_msg = f"Error! No reference annotation tracks left after filtering procedure! Abort"
        to_log(err_msg)
        die(err_msg)

    # write transcripts that passed the filter to the output file
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write("\n".join(new_lines) + "\n")
    f.close() if output != "stdout" else None

    if save_rejected:
        # save transcripts that didn't pass the filter + reason why
        f = open(save_rejected, "w")
        for elem in rejected:
            f.write(f"{elem[0]}\t{elem[1]}\n")
        f.close()


def main():
    """Entry point."""
    args = parse_args()
    prepare_bed_file(args.input, args.output, args.out_of_frame)
    sys.exit(0)


if __name__ == "__main__":
    main()
