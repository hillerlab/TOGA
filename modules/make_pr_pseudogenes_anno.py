#!/usr/bin/env python3
"""Create processed pseudogenes annotation track."""
import argparse
import sys
import os
from collections import defaultdict
try:
    from modules.common import chain_extract_id
    from modules.common import bed_extract_id
except ImportError:
    from common import chain_extract_id
    from common import bed_extract_id

PINK_COLOR = "250,50,200"
DEF_SCORE = 100


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_classification", help="Chains classification file.")
    app.add_argument("chain_file", help="Chain file")
    app.add_argument("bed_db", help="Bed indexed file")
    app.add_argument("output", help="Bed-9 output")
    app.add_argument("--verbose", "-v", action="store_true", dest="verbose",
                     help="Enable verbosity")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args_ = app.parse_args()
    return args_


def verbose(msg):
    """Verbose message."""
    sys.stderr.write(str(msg))
    sys.stderr.write("\n")


def get_pp_gene_chains(chain_class_file, v=False):
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
        pp_genes_field = line_data[4]
        if pp_genes_field == "0":
            # it 0 -> no ppgene chains -> skip
            continue
        # parse comma-separated string and save to dict
        pp_genes = [int(x) for x in pp_genes_field.split(",") if x != ""]
        gene_to_pp_chains[trans] = pp_genes
    f.close()
    if v:
        verbose(f"Extracted {len(gene_to_pp_chains)} genes with proc pseudogenes")
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


def get_corr_q_regions(gene_to_pp_chains, chain_file, chain_dict, bed_bdb, v=False):
    """Create projection: q region dict.
    
    We assume the following:
    1) processed pseudogene is a "single-exon" element
    2) ppgene chain covers the ppgene and nothing else.
    """
    proj_to_q_reg = {}  # save results here
    task_size = len(gene_to_pp_chains)
    # iterate over gene: [chain ids] elements
    for num, (gene, chains) in enumerate(gene_to_pp_chains.items()):
        if v:
            verbose(f"# Processing gene {gene} {num} / {task_size} with {len(chains)} chains")
        # extract gene track
        gene_track = bed_extract_id(bed_bdb, gene).rstrip().split("\t")
        gene_strand = gene_track[5]  # we need the strand only
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
    return proj_to_q_reg


def make_bed_lines(proj_to_reg):
    """Create bed9 lines."""
    bed_lines = []
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
        line_data = (chrom, chrom_start, chrom_end, name, DEF_SCORE, strand, thick_start, thick_end, PINK_COLOR)
        line = "\t".join(str(x) for x in line_data)
        bed_lines.append(line)
    return bed_lines


def save_bed(bed_lines, bed_file):
    """Save bed file."""
    f = open(bed_file, "w")
    f.write("\n".join(bed_lines))
    f.write("\n")
    f.close()


def load_chain_dict(chain_index_file):
    """Load dict chain ID: position in the file."""
    # TODO: put in common dir
    ans = {}
    if not os.path.isfile(chain_index_file):
        sys.exit(f"Error! File {chain_index_file} not found.")
    f = open(chain_index_file, "r")
    # tab-separated file
    # chain_ID, start_byte, offset
    for line in f:
        line_data = line.rstrip().split("\t")
        chain_id = int(line_data[0])
        start_byte = int(line_data[1])
        offset = int(line_data[2])
        val = (start_byte, offset)
        ans[chain_id] = val
    f.close()
    return ans


def create_ppgene_track(chain_class_file, chain_file, bed_bdb, output, v=None):
    """Create ppgene track."""
    # get processed pseudogene chains
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    verbose("Loading chains index...") if v else None
    index_file = chain_file.replace(".chain", ".chain_ID_position")
    chain_dict = load_chain_dict(index_file)
    verbose("Extracting pr pseudogene chains") if v else None
    gene_to_pp_chains = get_pp_gene_chains(chain_class_file, v)
    # for each gene-chain pair get corresponding region in query
    verbose("Extracting corresponding regions") if v else None
    projection_to_reg = get_corr_q_regions(gene_to_pp_chains, chain_file, chain_dict, bed_bdb, v)
    # convert the regions to bed-9 formatted-lines
    verbose("Saving results") if v else None
    bed_lines = make_bed_lines(projection_to_reg)
    verbose(f"There are {len(bed_lines)} bed lines") if v else None
    save_bed(bed_lines, output)  # write lines to the output file


if __name__ == "__main__":
    args = parse_args()
    create_ppgene_track(args.chain_classification,
                        args.chain_file,
                        args.bed_db,
                        args.output,
                        v=args.verbose)
