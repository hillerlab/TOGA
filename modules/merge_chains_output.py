#!/usr/bin/env python3
"""Parse raw chain runner output.

Chain features extraction steps results in numerous files.
This script merges these files and then
builds s table containing chain features.
"""
import argparse
import os
import sys
from datetime import datetime as dt
from collections import defaultdict
try:  # for robustness
    from modules.common import make_cds_track
    from modules.common import eprint
    from modules.common import die
    from modules.common import read_isoforms_file
except ImportError:
    from common import make_cds_track
    from common import eprint
    from common import die
    from common import read_isoforms_file


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

t0 = dt.now()


def verbose(msg, end="\n"):
    """Eprint for verbose messages."""
    eprint(msg + end)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("results_dir", type=str,
                     help="Directory containing the results.")
    app.add_argument("bed_file", type=str,
                     help="Bed file containing annotations for genes analyzed.")
    app.add_argument("output", type=str,
                     help="Save output here.")
    app.add_argument("--isoforms", "-i", type=str, default=None,
                     help="File with isoforms. Each line means "
                     "different ids of a gene, like gene_id<tab>[,-separated list of transcripts].")
    app.add_argument("--verbose", "-v", action="store_true",
                     dest="verbose", help="Verbose messages.")
    app.add_argument("--exon_cov_chains", "-e", action="store_true", dest="exon_cov_chains",
                     help="If yes, computing the gene_overs parameter we consider only "
                          "those chains that cover at least 1 base in exons. If no, "
                          "consider all the chains overlapping and introns only any.")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_bed_data(bed_file):
    """Get the necessary data from the bed file."""
    result = {}  # return this dictionary
    verbose(f"Reading {bed_file}")
    f = open(bed_file, "r")
    for line in f:
        # parse tab-separated bed file
        all_bed_info = line.rstrip().split("\t")
        cds_track = make_cds_track(line)  # we need CDS only
        cds_bed_info = cds_track.rstrip().split("\t")
    
        if len(all_bed_info) != 12 or len(cds_bed_info) != 12:
            # if there are not 12 fields - no guarantee that we parse what we want
            die(f"Error! Bed12 file {bed_file} is corrupted!")

        # extract fields that we need
        chromStart = int(all_bed_info[1])  # gene start
        chromEnd = int(all_bed_info[2])    # and end
        # blocks represent exons
        all_blockSizes = [int(x) for x in all_bed_info[10].split(',') if x != '']
        cds_blockSizes = [int(x) for x in cds_bed_info[10].split(',') if x != '']
        # data to save
        gene_len = abs(chromStart - chromEnd)

        # for decision tree I will need number of exons
        # and number of bases in exonic and intronic fractions
        exons_num = len(all_blockSizes)
        exon_fraction = sum(all_blockSizes)  # including UTR
        cds_fraction = sum(cds_blockSizes)  # CDS only
        intron_fraction = gene_len - exon_fraction
        gene_name = all_bed_info[3]

        # save the data
        result[gene_name] = {"gene_len": gene_len,
                             "exon_fraction": cds_fraction,
                             "intron_fraction": intron_fraction,
                             "exons_num": exons_num}
    f.close()
    verbose(f"Got data for {len(result.keys())} genes")
    return result


def process_gene_line(gene_line):
    """Parse gene line."""
    # contain [gene_id]=[intersected_chain_id] blocks
    # remove "gene" and "", get (gene. chain) tuples
    data_fields = [x for x in gene_line if "=" in x]
    data = [(x.split("=")[1], x.split("=")[0]) for x in data_fields]
    chain = data[0][0]  # chain is equal everywhere
    genes = [x[1] for x in data]
    return chain, genes


def parse_pairs(pairs):
    """Parse lines like X=5,Y=56, and returns a dict."""
    # ENST00000002501=0.1028238844578573,ENST00000006053=0.16846186988367085,
    # the last elem is always ""
    data = {x.split("=")[0]: x.split("=")[1] for x in pairs.split(",")[:-1]}
    return data


def process_chain_line(chain_line):
    """Parse chain line."""
    chain_id = chain_line[1]
    # collect the data in this dictionary
    data_fields = {chain_id: {}}
    # chain_id synteny global_score global_exo [local_exos] [coverage] [introns] [chain_len]

    # simple values:
    target = data_fields[chain_id]
    target["synteny"] = chain_line[2]
    target["global_score"] = chain_line[3]
    target["global_exo"] = chain_line[4]
    target["chain_Q_len"] = chain_line[5]
    target["chain_len"] = chain_line[10]

    # these values are more complex
    # there are comma-separated pairs of values (TRANSCRIPT, VALUE),
    target["local_exos"] = parse_pairs(chain_line[6])
    target["coverages"] = parse_pairs(chain_line[7])
    target["introns"] = parse_pairs(chain_line[8])
    target["flank_cov"] = parse_pairs(chain_line[9])
    return data_fields


def load_results(results_dir):
    """Load and sort the chain feature extractor results."""
    verbose("Loading the results...")
    results_files = os.listdir(results_dir)
    verbose(f"There are {len(results_files)} result files to combine")

    # to hold data from fields "genes":
    chain_genes_data = defaultdict(list)
    # to hold data from "chains" field:
    chain_raw_data = {}
    # read file-by-file, otherwise it takes too much place
    genes_counter, chain_counter = 0, 0  # count chain and genes lines

    for results_file in results_files:
        # there are N files: read them one-by-one
        path = os.path.join(results_dir, results_file)
        f = open(path, "r")
        for line in f:
            # read file line-by-line, all fields are tab-separated
            line_data = line.rstrip().split("\t")
            # define the class of this line
            # a line could be either gene or chain-related
            if line_data[0] == "genes":
                # process as a gene line
                chain, genes = process_gene_line(line_data)
                chain_genes_data[chain].extend(genes)
                genes_counter += 1
            elif line_data[0] == "chain":
                # chain related data
                the_chain_related = process_chain_line(line_data)
                # add this chain-related dict to the global one
                chain_raw_data.update(the_chain_related)
                chain_counter += 1
        # do not forget to close the file
        f.close()

    verbose(f"Got {len(chain_genes_data)} keys in chain_genes_data")
    verbose(f"Got {len(chain_raw_data)} keys in chain_raw_data")
    verbose(f"There were {genes_counter} genes lines and {chain_counter} chain lines")
    # actually, these values must be equal
    # just a sanity check
    if not genes_counter == chain_counter:
        eprint(f"WARNING! genes_counter and chain_counter hold different "
               f"values:\n{genes_counter} and {chain_counter} respectively")
        die("Some features extracting jobs died!")
    return chain_genes_data, chain_raw_data


def revert_dict(dct):
    """Revert {a: [list of b's]} to {b: [list of a's]}."""
    reverted = defaultdict(list)
    for k, value in dct.items():
        for v in value:
            reverted[v].append(k)
    verbose(f"chain_genes_data dict reverted, there are {len(reverted.keys())} keys now")
    return reverted


def make_synteny(genes, isoforms):
    """Return synteny for a list of genes and dict of isoforms."""
    return len(list(set([isoforms.get(gene) for gene in genes])))


def combine(bed_data, chain_data, genes_data, exon_cov, isoforms):
    """Combine chain and bed data and return gene-oriented table."""
    verbose("Combining the data...")
    combined = defaultdict(list)  # {gene: [related lines]}

    for chain, data in chain_data.items():
        # data dict contain the following keys:
        # 'synteny', 'global_score', 'global_exo', 'chain_len', 'local_exos', 'coverages', 'introns'
        # get the genes intersected this chain
        genes = list(data["local_exos"].keys())
        synteny = data["synteny"] if not isoforms else make_synteny(genes, isoforms)

        for gene in genes:  # iterate gene-by-gene
            # build a gene-oriented line
            gene_feat = bed_data.get(gene)
            if not gene_feat:
                continue  # it should not happen but...
            # extract gene-related features from this chain
            local_exo = data["local_exos"][gene]
            exon_coverage = str(data["coverages"][gene])
            intron_coverage = str(data["introns"][gene])
            flank_cov = str(data["flank_cov"][gene])
            # get a number of chains that covert this gene
            # else: if you need the chains that overlap EXONS
            gene_overs = len(genes_data[gene]) if not exon_cov else \
                len([x for x in genes_data[gene] if x != "None"])
            # fill the line
            line_data = [  # gene and chain overlap information
                         gene, gene_overs, chain, synteny,
                         # global chain features
                         data["global_score"], data["global_exo"],
                         data["chain_len"], data["chain_Q_len"],
                         # chain to gene local features
                         local_exo, exon_coverage, intron_coverage,
                         # gene features
                         gene_feat["gene_len"], gene_feat["exons_num"],
                         gene_feat["exon_fraction"], gene_feat["intron_fraction"],
                         flank_cov]

            line = "\t".join([str(x) for x in line_data]) + "\n"
            combined[gene].append(line)
    # got all the data --> return back
    verbose(f"Got combined dict with {len(combined.keys())} keys")
    return combined


def save(data, output):
    """Save the data into the file."""
    # make the header
    header_fields = "gene gene_overs chain synt gl_score gl_exo chain_len cds_qlen loc_exo exon_cover " \
                    "intr_cover gene_len ex_num ex_fract intr_fract flank_cov".split()
    header = "\t".join(header_fields) + "\n"
    # define the stream to write the data
    verbose("Writing output to {}".format(output))
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write(header)
    for lines in data.values():
        f.write("".join(lines))
    f.close()


def merge_chains_output(bed_file, isoforms_file, results_dir,
                        output, exon_cov_chains=False):
    """Chains output merger core function."""
    # read bed file, get gene features
    bed_data = read_bed_data(bed_file)
    # load isoforms data if provided
    # returns 3 values, we keep only isoform to gene
    if isoforms_file:
        _, isoforms, _ = read_isoforms_file(isoforms_file)
    else:
        isoforms = None
    # read result files from unit
    chain_genes_data, chain_raw_data = load_results(results_dir)
    # I need this dict reverted actually
    # not chain-genes-data but gene-chains-data
    genes_data = revert_dict(chain_genes_data)

    # combine all the data into one gene-oriented dictionary
    combined_data = combine(bed_data,
                            chain_raw_data,
                            genes_data,
                            exon_cov_chains,
                            isoforms)
    # save this data
    save(combined_data, output)
    # finish the program
    eprint(f"Estimated_time: {format(dt.now() - t0)}")


def main():
    """Entry point."""
    args = parse_args()
    merge_chains_output(args.bed_file,
                        args.isoforms,
                        args.results_dir,
                        args.output,
                        args.exon_cov_chains)


if __name__ == "__main__":
    main()
