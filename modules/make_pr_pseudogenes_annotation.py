#!/usr/bin/env python3
"""Create processed pseudogenes annotation track."""
import argparse
import sys
import os
from collections import defaultdict
from version import __version__

# TODO: figure this out
try:
    from modules.common import chain_extract_id
    from modules.common import bed_extract_id
    from modules.common import load_chain_dict
    from modules.common import setup_logger
    from modules.common import to_log
except ImportError:
    from common import chain_extract_id
    from common import bed_extract_id
    from common import load_chain_dict
    from common import setup_logger
    from common import to_log

__author__ = "Bogdan M. Kirilenko"

PINK_COLOR = "250,50,200"
DEF_SCORE = 100
PPGENE_SCORE_MARK = "-2.0"


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_classification", help="Chains classification file.")
    app.add_argument("chain_file", help="Chain file")
    app.add_argument("bed_db", help="Bed indexed file")
    app.add_argument("output", help="Bed-9 output")
    app.add_argument("--log_file", default=None, help="Log file")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args_ = app.parse_args()
    return args_


def verbose(msg):
    """Verbose message."""
    sys.stderr.write(str(msg))
    sys.stderr.write("\n")


def get_pp_gene_chains(chain_class_file):
    """Get gene: pp chains dict."""
    gene_to_pp_chains = defaultdict(list)  # init the dict
    f = open(chain_class_file, "r")  # open file with classifications
    f.__next__()  # skip header
    for line in f:
        line_data = line.rstrip().split("\t")
        # line contains the following fields"
        # gene orthologs paralogs trans p_pseudogenes
        trans = line_data[0]
        # proc_pseudogene chains are in the 4th field
        # pp_genes_field = line_data[4]
        # if pp_genes_field == "0":
        #     # it 0 -> no ppgene chains -> skip
        #     continue
        # parse comma-separated string and save to dict
        # pp_genes = [int(x) for x in pp_genes_field.split(",") if x != ""]
        # gene_to_pp_chains[trans] = pp_genes
        chain = int(line_data[1])
        stat = line_data[2]
        if stat == PPGENE_SCORE_MARK:
            gene_to_pp_chains[trans].append(chain)

    f.close()
    to_log(f"make_pr_pseudogenes anno: {len(gene_to_pp_chains)} transcripts have processed pseudogenes")
    return gene_to_pp_chains


def extract_chain(chain_file, chain_dict, chain):
    """Extract chain string.

    We have: chain file, chain_id, start byte and offset.
    """
    f = open(chain_file, "rb")
    start, offset = chain_dict.get(int(chain))
    f.seek(start)  # jump to start_byte_position
    chain = f.read(offset).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def get_corr_q_regions(gene_to_pp_chains, chain_file, chain_dict, bed_bdb):
    """Create projection: q region dict.

    We assume the following:
    1) processed pseudogene is a "single-exon" element
    2) ppgene chain covers the ppgene and nothing else.
    """
    proj_to_q_reg = {}  # save results here
    task_size = len(gene_to_pp_chains)
    to_log(f"make_pr_pseudogenes anno: identifying processed pseudogenes regions in the query")
    to_log(f"make_pr_pseudogenes anno: analysing data for {task_size} reference transcripts")
    chains_counter = 0
    # iterate over gene: [chain ids] elements
    for num, (gene, chains) in enumerate(gene_to_pp_chains.items()):
        # extract gene track
        gene_track = bed_extract_id(bed_bdb, gene).rstrip().split("\t")
        gene_strand = gene_track[5]  # we need the strand only
        chains_counter += len(chains)
        for chain_id in chains:
            # we have a list of chains
            projection = f"{gene}.{chain_id}"  # name this projection as usual
            # extract the chain and parse it's header
            # chain_body = chain_extract_id(chain_bdb, chain_id)
            chain_body = extract_chain(chain_file, chain_dict, chain_id)
            chain_header = chain_body.split("\n")[0].split()
            # we need chrom, start, end, strand and q_size
            q_chrom = chain_header[7]
            q_size = int(chain_header[8])
            q_strand = chain_header[9]
            q_start = int(chain_header[10])
            q_end = int(chain_header[11])
            if q_strand == "-":
                # if q_strand is - we need to invert coordinates
                # see chains documentation for details
                t_ = q_start
                q_start = q_size - q_end
                q_end = q_size - t_
            # get projection strand and sage the region
            proj_strand = "+" if q_strand == gene_strand else "-"
            proj_reg = (q_chrom, proj_strand, q_start, q_end)
            proj_to_q_reg[projection] = proj_reg
    result_size = len(proj_to_q_reg)
    to_log(
        f"make_pr_pseudogenes anno: for {task_size} reference transcripts "
        f"found {chains_counter} processed pseudogene chains in total"
    )
    to_log(f"make_pr_pseudogenes anno: resulted in {result_size} processed pseudogenes projections in query")
    return proj_to_q_reg


def merge_intervals(intervals):
    # Sort the intervals by start position
    to_log(f"make_pr_pseudogenes anno: collapsing intersecting bed lines (if any)")
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    to_log(f"make_pr_pseudogenes anno: {len(intervals)} records before collapse")

    # Initialize the merged list with the first interval
    merged = [intervals[0]]

    for current in intervals[1:]:
        # Get the last interval in the merged list
        last = merged[-1]

        # If the current interval overlaps with the last one and they are on the same chromosome and strand,
        # merge them by extending the end of the last interval to the end of the current one
        if current[0] == last[0] and current[5] == last[5] and current[1] <= last[2]:
            last = list(last)
            last[2] = max(last[2], current[2])
            merged[-1] = tuple(last)
        else:
            # Otherwise, just add the current interval to the list
            merged.append(current)
    to_log(f"make_pr_pseudogenes anno: {len(merged)} records after collapse")
    return merged


def make_bed_lines(proj_to_reg):
    """Create bed9 lines."""
    task_size = len(proj_to_reg)
    to_log(f"make_pr_pseudogenes_anno: generating bed lines for {task_size} projections")
    intervals = []
    for name, region in proj_to_reg.items():
        # just parse the region
        chrom = region[0]
        strand = region[1]
        chrom_start = region[2]
        chrom_end = region[3]
        # non coding region: then thick_start = thick_end
        thick_start = chrom_start
        thick_end = chrom_start
        # we don't need score field -> so there is some default value
        line_data = (
            chrom,
            chrom_start,
            chrom_end,
            f"{name}-like",
            DEF_SCORE,
            strand,
            thick_start,
            thick_end,
            PINK_COLOR,
        )
        intervals.append(line_data)

    # Merge overlapping intervals
    to_log(f"make_pr_pseudogenes_anno: generated {len(intervals)} intervals")
    intervals = merge_intervals(intervals)

    # Convert the merged intervals back into bed lines
    bed_lines = ["\t".join(str(x) for x in interval) for interval in intervals]
    to_log(f"make_pr_pseudogenes_anno: generated {len(bed_lines)} bed lines to be saved")
    return bed_lines


def save_bed(bed_lines, bed_file):
    """Save bed file."""
    f = open(bed_file, "w")
    f.write("\n".join(bed_lines))
    f.write("\n")
    f.close()


def create_ppgene_track(chain_class_file, chain_file, bed_bdb, output):
    """Create ppgene track."""
    # get processed pseudogene chains
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    to_log(f"make_pr_pseudogenes_anno: loading chain index...")
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    chain_dict = load_chain_dict(index_file)
    gene_to_pp_chains = get_pp_gene_chains(chain_class_file)
    if len(gene_to_pp_chains) == 0:
        # no proc pseudogenes
        to_log(f"make_pr_pseudogenes_anno: no processed pseudogenes found, skip")
        return
    # for each gene-chain pair get corresponding region in query
    projection_to_reg = get_corr_q_regions(
        gene_to_pp_chains, chain_file, chain_dict, bed_bdb
    )
    # convert the regions to bed-9 formatted-lines
    bed_lines = make_bed_lines(projection_to_reg)
    to_log(f"make_pr_pseudogenes_anno: {len(bed_lines)} processed pseudogene annotations in total")
    to_log(f"make_pr_pseudogenes_anno: saving results to {output}")
    save_bed(bed_lines, output)  # write lines to the output file


if __name__ == "__main__":
    args = parse_args()
    setup_logger(args.log_file)
    create_ppgene_track(
        args.chain_classification,
        args.chain_file,
        args.bed_db,
        args.output
    )
