#!/usr/bin/env python3
"""Join projection in genes for query annotation.

The rules are simple:
If two transcript have an intersecting exon
then they belong to the same gene.
Exons should be on the same strand!
"""
import argparse
import sys
import os
from collections import defaultdict
import networkx as nx
try:
    from modules.common import flatten
except ImportError:
    from common import flatten


def parse_args():
    """Read CMD args."""
    app = argparse.ArgumentParser()
    app.add_argument("query_bed", help="Query annotation bed file.")
    app.add_argument("output", help="Output containing query isoforms data")
    app.add_argument("--genes_track", help="Save gene borders track")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_query_bed(bed_file):
    """Read query bed, return exons and exon to gene dict."""
    f = open(bed_file, "r")
    exon_counter = 0  # each exon gets an unique ID
    exons_list = []  # save exons here
    exon_id_to_trans = {}  # exon_id to corresponding transcript
    trans_to_range = {}  # also save a region for each trascript

    for line in f:
        line_data = line.rstrip().split("\t")
        # parse bed 12 file according to specification
        chrom = line_data[0]
        transcript_name = line_data[3]
        strand = line_data[5]
        # we need CDS start, end == thickStart & thickEnd
        thickStart = int(line_data[6])
        thickEnd = int(line_data[7])
        blockCount = int(line_data[9])
        blockSizes = [int(x) for x in line_data[10].split(',') if x != '']
        blockStarts = [int(x) for x in line_data[11].split(',') if x != '']
        # we can now save transcript range:
        trans_range = (chrom, strand, thickStart, thickEnd)
        trans_to_range[transcript_name] = trans_range

        for i in range(blockCount):
            # get absolute coordinates for each exon
            exon_counter += 1  # update ID for new exon
            exon_abs_start = blockStarts[i] + thickStart
            exon_abs_end = exon_abs_start + blockSizes[i]
            # summarize exon data:
            exon_data = (exon_counter, chrom, strand, exon_abs_start, exon_abs_end)
            exons_list.append(exon_data)  # add it to the list
            # save exon_id to transcript:
            exon_id_to_trans[exon_counter] = transcript_name
    f.close()
    return exons_list, exon_id_to_trans, trans_to_range


def split_exons_in_chr_dir(exons_list):
    """Make a (chr, strand): [exons] dict."""
    chr_dir_exons_not_sorted = defaultdict(list)
    for exon in exons_list:
        # just rearrange a data a bit
        # go exon-by-exon
        exon_id = exon[0]
        chrom = exon[1]
        strand = exon[2]
        start = exon[3]
        end = exon[4]
        # split into two elements: key (chr & directon) and value (exon range + exon ID)
        chr_dir = (chrom, strand)
        exon_reduced = (exon_id, start, end)
        # add this to default dict
        chr_dir_exons_not_sorted[chr_dir].append(exon_reduced)
    # sort exons in each chr_dir track -> for efficiend algorithm
    # use start (field 1) as key
    chr_dir_exons = {k: sorted(v, key=lambda x: x[1])
                     for k, v in chr_dir_exons_not_sorted.items()}
    return chr_dir_exons


def intersect_ranges(region_1, region_2):
    """Return TRUE if ranges intersect."""
    inter_size = min(region_1[1], region_2[1]) - max(region_1[0], region_2[0])
    return True if inter_size > 0 else False


def intersect_exons(chr_dir_exons, exon_id_to_transcript):
    """Create graph of connected exons.
    
    We will use a graph where nodes are transcript IDs.
    Nodes are connected if they have a pair of intersected exons.
    """
    G = nx.Graph()  # init the graph
    for exons in chr_dir_exons.values():
        # this is the same chrom and direction now
        # add nodes now: to avoid missing transcripts
        # that dont intersect anything
        exon_ids = [e[0] for e in exons]
        transcripts = set(exon_id_to_transcript[x] for x in exon_ids)
        G.add_nodes_from(transcripts)

        if len(exons) == 1:
            # only one exon: nothing to intersect
            # gene id is on graph now: we can just continue
            continue
        exons_num = len(exons)
        for i in range(exons_num - 1):
            # nearly N^2 algorithm (a nested loop)
            # but actually faster
            # pick exon i
            exon_i = exons[i]
            exon_i_id = exon_i[0]
            exon_i_start = exon_i[1]
            exon_i_end = exon_i[2]
            # get range and corresponding transcript
            i_range = (exon_i_start, exon_i_end)
            i_trans = exon_id_to_transcript[exon_i_id]
            # then let's look at next exons
            for j in range(i + 1, exons_num):
                # exon j exists for sure (first range is up to exons_num - 1)
                # num of exons is > 1
                exon_j = exons[j]
                exon_j_id = exon_j[0]
                exon_j_start = exon_j[1]
                exon_j_end = exon_j[2]
                if exon_j_start > exon_i_end:
                    # if exon_j start is > exon_i end then we can break
                    # they are sorted: for all next exons start would be > exon_i end
                    # and if so, they do not intersect for sure
                    break
                # let's check that they intersect
                j_range = (exon_j_start, exon_j_end)
                j_trans = exon_id_to_transcript[exon_j_id]  # j exon transcript
                they_intersect = intersect_ranges(i_range, j_range)
                if they_intersect is False:
                    # don't intersect: nothing to do
                    continue
                # if we are in this branch: exons intersect
                # we can add an edge connectiong their corresponding transcripts
                G.add_edge(i_trans, j_trans)
    return G


def parse_components(components, trans_to_range):
    """Get genes data.

    Each gene has the following data:
    1) It's ID
    2) Included transcripts
    3) Genic range.
    """
    genes_data = []  # save gene objects here
    for num, component in enumerate(components, 1):
        gene_id = f"reg_{num}"  # need to name them somehow
        # get transcripts and their ranges
        transcripts = set(component.nodes())
        regions = [trans_to_range[t] for t in transcripts]
        # define gene range, chrom and strand are same everywhere
        # just get them from the 0'st region
        reg_0 = regions[0]
        gene_chrom = reg_0[0]
        gene_strand = reg_0[1]
        # start and end -> just min of starts and max of ends
        starts = [r[2] for r in regions]
        ends = [r[3] for r in regions]
        gene_start = min(starts)
        gene_end = max(ends)
        # pack everything
        gene_range = (gene_chrom, gene_strand, gene_start, gene_end)
        gene_data = {"ID": gene_id, "range": gene_range, "transcripts": transcripts}
        genes_data.append(gene_data)
    return genes_data


def save_isoforms(genes_data, output):
    """Save isoforms data."""
    # use the same logic as Ensembl Biomart
    # two tab-separated columns: gene_id and transcript_id
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write("Region_ID\tProjection_ID\n")
    for gene_datum in genes_data:
        gene_id = gene_datum["ID"]
        transcripts = gene_datum["transcripts"]
        for transcript in transcripts:
            f.write(f"{gene_id}\t{transcript}\n")
    f.close() if output != "stdout" else None


def save_regions(genes_data, output):
    """Save gene ranges into a bed-6 file, if required."""
    if output is None:
        # no, not required
        return
    f = open(output, "w") if output != "stdout" else sys.stdout
    for gene_datum in genes_data:
        gene_id = gene_datum["ID"]
        gene_range = gene_datum["range"]
        chrom, strand, start, end = gene_range
        f.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")
    f.close() if output != "stdout" else None


def get_query_isoforms_data(query_bed, query_isoforms, save_genes_track=None):
    """Create isoforms track for query."""
    # extract all exons
    # exon could be described as: <ID, chrom, strand, start, end>
    # + table exon_ID -> corresponding gene
    exons_list, exon_id_to_transcript, trans_to_range = read_query_bed(query_bed)
    # get {(chr, dir) -> [exons]} dict (sorted)
    # exons from different chrom/direction cannot actually intersect
    chr_dir_to_exons = split_exons_in_chr_dir(exons_list)
    # the main part -> get exon intersections
    # create a graph of transcripts with intersected exons
    conn_graph = intersect_exons(chr_dir_to_exons, exon_id_to_transcript)
    # split graph into connected components
    # if two transcripts are in the same component -> they belong to the same gene
    nx_v = float(nx.__version__)
    if nx_v < 2.4:  # TODO: keep it for ~2 months, then remove deprecated branch
        components = list(nx.connected_component_subgraphs(conn_graph))
        msg = f"Warning! You use networkx v{nx_v}; splitting components with " \
               "nx.connected_component_subgraphs(), which is deprecated. " \
               "Please upgrade networks to supress this warning\n"
        sys.stderr.write(msg)
    else:
        components = [conn_graph.subgraph(c) for c in nx.connected_components(conn_graph)]
    # covert components to isoforms table
    genes_data = parse_components(components, trans_to_range)
    # save the results
    save_isoforms(genes_data, query_isoforms)
    save_regions(genes_data, save_genes_track)

if __name__ == "__main__":
    args = parse_args()
    get_query_isoforms_data(args.query_bed, args.output, save_genes_track=args.genes_track)
