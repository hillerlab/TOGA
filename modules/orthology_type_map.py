#!/usr/bin/env python3
"""For a target and queue files given make orthology map.

Create a table that shows orthology relationships between
reference and query genes. Divide them as one2one, one2many, etc.
"""
import argparse
import sys
from collections import defaultdict
from collections import Counter
import networkx as nx
try:
    from modules.common import split_proj_name
    from modules.common import flatten
    from modules.common import eprint
    from modules.common import die
    from modules.common import get_graph_components
except ImportError:
    from common import split_proj_name
    from common import flatten
    from common import eprint
    from common import die
    from common import get_graph_components

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

INCLUDE_CLASSES = {"I", "PI", "G"}
ONE2ZERO = "one2zero"
ONE2ONE = "one2one"
ONE2MANY = "one2many"
MANY2ONE = "many2one"
MANY2MANY = "many2many"


def read_isoforms(isoforms_file, transcripts):
    """Read isoforms data."""
    gene_to_transcripts = defaultdict(list)
    transcript_to_gene = {}
    if isoforms_file is None:
        # special branch: one gene - one transcript
        transcript_to_gene = {x: x for x in transcripts}
        # gene to transcripts is dict key : list
        gene_to_transcripts = {x: [x, ] for x in transcripts}
        return gene_to_transcripts, transcript_to_gene
    # need to read a file
    f = open(isoforms_file, "r")
    f.__next__()  # first line is header
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        trans = line_data[1]
        if trans not in transcripts:
            # exclude unnecessary transcripts
            # such as filtered out or skipped due to memory reasons
            continue
        transcript_to_gene[trans] = gene
        gene_to_transcripts[gene].append(trans)
    return gene_to_transcripts, transcript_to_gene


def read_paralogs(paralogs_file):
    """Read list of paralogs."""
    if paralogs_file is None:
        # then nothing to return
        return set()
    f = open(paralogs_file, "r")
    paralogs = set(x.rstrip() for x in f.readlines())
    f.close()
    return paralogs


def extract_names_from_bed(bed_file):
    """Just extract a list of names from bed-12 file."""
    with open(bed_file, "r") as f:
        names = set(x.split("\t")[3] for x in f.readlines())
    assert len(names) > 0  # if not, there is a huge bug
    return names


def read_loss_data(loss_file):
    """Read loss classification file."""
    proj_to_status = {}
    f = open(loss_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        # GLP data: entry type (projection, trans or gene), name ans status
        entry_type = line_data[0]
        if entry_type != "PROJECTION":
            # we need only projections data
            continue
        # save projection status
        proj = line_data[1]
        status = line_data[2]
        proj_to_status[proj] = status
    f.close()
    return proj_to_status


def filter_query_transcripts(transcripts, paralogs, trans_to_status):
    """Keep intact orthologous transcripts."""
    filt_transcripts = []  # save transcripts we keep
    print(f"Extracted {len(transcripts)} query transcripts")
    for trans in transcripts:
        is_paral = trans in paralogs
        if is_paral:
            # paralogs do not participate in this game
            continue
        l_status = trans_to_status.get(trans, "N")
        if l_status not in INCLUDE_CLASSES:
            # only grey and intact participate
            continue
        filt_transcripts.append(trans)
    print(f"After filters {len(filt_transcripts)} transcripts left")
    return set(filt_transcripts)


def get_t_trans_to_projections(que_projections):
    """Make transcript: projections dict."""
    trans_to_proj = defaultdict(list)
    for proj in que_projections:
        # just need to split projection name to get transcript ID
        trans, _ = split_proj_name(proj)
        trans_to_proj[trans].append(proj)
    return trans_to_proj


def connect_genes(t_trans_to_gene, t_trans_to_q_proj, q_proj_to_q_gene):
    """Create orthology relationships graph.
    
    Nodes are reference and query genes.
    If nodes are connected -> they are orthologs.
    """
    o_graph = nx.Graph()  # init the graph
    genes = set(t_trans_to_gene.values())
    o_graph.add_nodes_from(genes)
    q_genes_all = []
    print(f"Added {len(genes)} reference genes on graph")
    conn_count = 0  # connections counter
    for t_trans, t_gene in t_trans_to_gene.items():
        # we know ref gene: ref transcripts relation by default
        # also we already know transcript-to-projections (trivial)
        projections = t_trans_to_q_proj[t_trans]
        # from projections (query transcripts) we know query gene ids
        q_genes = set(q_proj_to_q_gene[p] for p in projections)
        for q_gene in q_genes:
            # we can connect all of those q_genes with the reference gene
            conn_count += 1
            q_genes_all.append(q_gene)
            o_graph.add_edge(t_gene, q_gene)
    q_genes_all = set(q_genes_all)
    print(f"Added {len(q_genes_all)} query genes on graph")
    print(f"Graph contains {conn_count} connections")
    return o_graph


def extract_orth_connections(graph, r_genes_all, q_genes_all):
    """Split graph in orth connections."""
    orth_connections = []
    graph_components = get_graph_components(graph)

    for component in graph_components:
        # each component contains some reference and query gene ids
        # and they are orthologs
        nodes = component.nodes()
        # split all nodes into ref and query gene ids:
        r_genes = [n for n in nodes if n in r_genes_all]
        q_genes = [n for n in nodes if n in q_genes_all]
        # to define the class we need sizes of these groups:
        r_len = len(r_genes)
        q_len = len(q_genes)
        if q_len == 0:
            c_class = ONE2ZERO
        elif r_len == 1 and q_len == 1:
            c_class = ONE2ONE
        elif r_len == 1 and q_len > 1:
            c_class = ONE2MANY
        elif r_len > 1 and q_len == 1:
            c_class = MANY2ONE
        elif r_len > 1 and q_len > 1:
            c_class = MANY2MANY
        else:
            raise RuntimeError("impossible")
        # create connection object and save it
        conn = {"r_genes": r_genes, "q_genes": q_genes, "c_class": c_class}
        orth_connections.append(conn)
    # count different orthology classes
    class_list = [c["c_class"] for c in orth_connections]
    print(f"Detected {len(class_list)} orthology components")
    class_count = Counter(class_list)
    print(f"Orthology class sizes:")
    for k, v in class_count.most_common():
        print(f"{k}: {v}")
    return orth_connections


def save_data(orth_connections, r_gene_to_trans, q_trans_to_gene, t_trans_to_projections, out):
    """Save orthology data."""
    f = open(out, "w") if out != "stdout" else sys.stdout
    # write the header
    f.write("t_gene\tt_transcript\tq_gene\tq_transcript\torthology_class\n")
    non_orthologous_isoforms = []  # to save reference transcripts that don't have orthologs

    for conn in orth_connections:
        # connection: class + GENES
        conn_class = conn["c_class"]
        ref_genes = conn["r_genes"]
        # que_genes = conn["q_genes"]
        # one line per one isoforms line
        if conn_class == ONE2ZERO:
            # special case, nothing to show
            for ref_gene in ref_genes:
                # line per transcript
                ref_transcripts = r_gene_to_trans[ref_gene]
                for ref_transcript in ref_transcripts:
                    f.write(f"{ref_gene}\t{ref_transcript}\tNone\tNone\t{ONE2ZERO}\n")
            # goto the next component
            continue
        # not one2zero
        for ref_gene in ref_genes:
            # there could be > 1 reference gene
            ref_transcripts = r_gene_to_trans[ref_gene]
            for ref_transcript in ref_transcripts:
                # also, there might be > 1 transcript per gene
                projections = t_trans_to_projections[ref_transcript]
                if not projections:
                    # probably, due to some reasons some transcripts don't have any
                    # associated orthology projection
                    # maybe it's lost (and gene has an intact isoform), or was skipped due to
                    # technical reasons
                    non_orthologous_isoforms.append(ref_transcript)
                for proj in projections:
                    # and last, each transcript can be projected > once
                    proj_q_gene = q_trans_to_gene[proj]
                    f.write(f"{ref_gene}\t{ref_transcript}\t{proj_q_gene}\t{proj}\t{conn_class}\n")
    # close file and return non-orthologous reference isoforms
    f.close() if out != "stdout" else None
    return non_orthologous_isoforms


def orthology_type_map(ref_bed, que_bed, out, ref_iso=None, que_iso=None,
                       paralogs_arg=None, loss_data=None, save_skipped=None):
    """Make orthology classification track."""
    q_trans_paralogs = read_paralogs(paralogs_arg)  # do not include paralogs
    ref_transcripts = extract_names_from_bed(ref_bed)  # list of reference transcripts
    que_transcripts_all = extract_names_from_bed(que_bed)  # list of query transcripts
    if loss_data is None:
        # no loss data: assume everything is intact
        trans_to_L_status = {x: "I" for x in que_transcripts_all}
    else:  # we do add only I, PI and G projections
        trans_to_L_status = read_loss_data(loss_data)
    # remove I/PI/G or paralogous transcripts
    que_transcripts = filter_query_transcripts(que_transcripts_all,
                                               q_trans_paralogs,
                                               trans_to_L_status)
    # read reference and query isoform files; orthology is a story about genes
    r_gene_to_trans, r_trans_to_gene = read_isoforms(ref_iso, ref_transcripts)
    q_gene_to_trans, q_trans_to_gene = read_isoforms(que_iso, que_transcripts_all)
    r_genes_all = set(r_gene_to_trans.keys())
    q_genes_all = set(q_gene_to_trans.keys())
    # make transcript to projections dict:
    t_trans_to_projections = get_t_trans_to_projections(que_transcripts)
    # create graph to connect orthologous transcripts
    o_graph = connect_genes(r_trans_to_gene, t_trans_to_projections, q_trans_to_gene)
    # if a group of reference and query genes are in the same connected component
    # then they are orthologs
    orth_connections = extract_orth_connections(o_graph, r_genes_all, q_genes_all)
    # save data, get list of transcript that were not projected
    not_saved = save_data(orth_connections, r_gene_to_trans, q_trans_to_gene, t_trans_to_projections, out)
    if save_skipped:
        # save transcripts of orthologous genes that don't have any projection
        f = open(save_skipped, "w")
        f.write("\n".join(not_saved))
        f.write("\n")
        f.close()


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("target_bed", help="Annotation for the target genome.")
    app.add_argument("query_bed", help="Query annotation file.")
    app.add_argument("output", help="Output with orthology information")
    app.add_argument("--ref_isoforms", "--ri", type=str, default=None,
                     help="Isoforms data for reference (if exists)")
    app.add_argument("--que_isoforms", "--qi", type=str, default=None,
                     help="Isoforms data for query (if exists)")
    app.add_argument("--paralogs", "-p", default=None, help="List of paralogs")
    app.add_argument("--loss_data", "-l", default=None, help="GLP classification")
    app.add_argument("--save_skipped", "-s", default=None, help="Save orphan transcripts")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def main():
    """CLI entry point."""
    args = parse_args()
    orthology_type_map(args.target_bed,
                       args.query_bed,
                       args.output,
                       ref_iso=args.ref_isoforms,
                       que_iso=args.que_isoforms,
                       paralogs_arg=args.paralogs,
                       loss_data=args.loss_data,
                       save_skipped=args.save_skipped)


if __name__ == "__main__":
    main()
