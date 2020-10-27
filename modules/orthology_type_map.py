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
    from modules.common import read_isoforms_file
except ImportError:
    from common import split_proj_name
    from common import flatten
    from common import eprint
    from common import die
    from common import get_graph_components
    from common import read_isoforms_file


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

INCLUDE_CLASSES = {"I", "PI", "G"}
SCORE_THR = 0.9
ONE2ZERO = "one2zero"
ONE2ONE = "one2one"
ONE2MANY = "one2many"
MANY2ONE = "many2one"
MANY2MANY = "many2many"

R_GENES = "r_genes"
Q_GENES = "q_genes"
C_CLASS = "c_class"


def read_isoforms__otm(isoforms_file, transcripts):
    """Read isoforms data.

    Extended orthology type map version."""
    # gene_to_transcripts = defaultdict(list)
    # transcript_to_gene = {}
    if isoforms_file is None:
        # special branch: one gene - one transcript
        transcript_to_gene = {x: x for x in transcripts}
        # gene to transcripts is dict key : list
        gene_to_transcripts = {x: [x, ] for x in transcripts}
        return gene_to_transcripts, transcript_to_gene
    gene_to_transcripts, transcript_to_gene, _ = read_isoforms_file(isoforms_file,
                                                                    pre_def_trans_list=transcripts)
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


def read_proj_scores(orth_score_file):
    """Return projection: orthology score dict."""
    ret = {}
    if orth_score_file is None:
        return ret
    f = open(orth_score_file, "r")
    f.__next__()  # skip header
    for line in f:
        line_datum = line.rstrip().split("\t")
        trans = line_datum[0]
        chain = line_datum[1]
        projection = f"{trans}.{chain}"
        score = float(line_datum[2])
        if score == -1:
            # trans chain: without score
            # assign default 0.5 value
            score = 0.5
        ret[projection] = score
    f.close()
    return ret


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
        # is_low_qual = trans in low_score_trans
        # if is_low_qual:
        #     # chain orthology score is quite low
        #     continue
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


def connect_genes(t_trans_to_gene, t_trans_to_q_proj, q_proj_to_q_gene, proj_to_score):
    """Create orthology relationships graph.
    
    Nodes are reference and query genes.
    If nodes are connected -> they are orthologs.
    Also emit ref_gene -> que_gene -> scores dict.
    """
    o_graph = nx.Graph()  # init the graph
    ref_que_gene_scores = defaultdict(list)
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
        # also for each t_gene -> q_gene: add score
        q_genes = set()
        for proj_ in projections:
            q_gene = q_proj_to_q_gene[proj_]
            q_genes.add(q_gene)
            o_score = proj_to_score.get(proj_, 0.0)
            score_dict_key = (t_gene, q_gene)
            ref_que_gene_scores[score_dict_key].append(o_score)
        # add connections to graph
        for q_gene in q_genes:
            # we can connect all of those q_genes with the reference gene
            conn_count += 1
            q_genes_all.append(q_gene)
            o_graph.add_edge(t_gene, q_gene)
    q_genes_all = set(q_genes_all)
    print(f"Added {len(q_genes_all)} query genes on graph")
    print(f"Graph contains {conn_count} connections")
    return o_graph, ref_que_gene_scores


def get_c_class(r_num, q_num):
    """Given numbers of reference and query genes return orthology class."""
    # not many2many branches
    if q_num == 0:
        return ONE2ZERO
    elif r_num == 1 and q_num == 1:
        return ONE2ONE
    elif r_num == 1 and q_num > 1:
        return ONE2MANY
    elif r_num > 1 and q_num == 1:
        return MANY2ONE
    elif r_num > 1 and q_num > 1:
        return MANY2MANY
    else:  # something that must never happen
        err_msg = f"Orthology type map: reached corrupt orthology component\n"\
                  f"Got numbers of ref genes: {r_num}; query: {q_num}"\
                  f"Abort."
        raise RuntimeError(err_msg)


def is_complete_bipartite(graph, r_genes, q_genes):
    """Return True if a complete bipartite graph provided.

    Here we assume that the graph is bipartite.
    The graph is complete if:
    1) Each reference gene appears in #q_genes connections,
    2) Each query gene appears in #r_genes connections."""
    num_r_genes = len(r_genes)
    num_q_genes = len(q_genes)
    connections_flat = flatten(graph.edges())
    conn_count = Counter(connections_flat)
    # check criteria 1 and 2
    all_ref_ok = all(conn_count[r] == num_q_genes for r in r_genes)
    all_que_ok = all(conn_count[q] == num_r_genes for q in q_genes)
    if all_ref_ok and all_que_ok:
        # criteria 1 and 2 satisfied
        return True
    # something's not satisfied -> is not complete
    return False


def order_edges(edges_lst, ref_nodes):
    """Reorder edges such as reference node is on the first position."""
    ret = []
    for edge in edges_lst:
        f_, s_ = edge
        if f_ in ref_nodes:
            ref_node = f_
            que_node = s_
        else:
            ref_node = s_
            que_node = f_
        item = (ref_node, que_node)
        ret.append(item)
    return ret


def edges_to_dicts(graph, ref_nodes__set):
    """Convert list of tuples to dicts."""
    # in tuple ref and que may go in the different order
    # this is why we cannot convert this list to dict directly
    # need to check what on which position in each tuple first
    edges = order_edges(graph.edges(), ref_nodes__set)
    ref_nodes_edges = defaultdict(list)
    que_nodes_edges = defaultdict(list)
    for edge in edges:
        ref_nodes_edges[edge[0]].append(edge[1])
        que_nodes_edges[edge[1]].append(edge[0])
    return ref_nodes_edges, que_nodes_edges


def graph_is_empty(graph):
    """Check whether graph is empty."""
    return True if len(graph.nodes()) == 0 else False


def graph_has_isolated_points(graph):
    """Check whether graph has isolated nodes."""
    return True if len(list(nx.isolates(graph))) > 0 else False


def select_strongly_connected_ref_genes(ref_conns, que_conns):
    """Select strongly connected reference genes.

    Strongly connected genes have >1 common connections.
    """
    ret = {}
    for ref_gene, queries in ref_conns.items():
        ref_second_order = Counter(flatten(que_conns[q] for q in queries))
        del ref_second_order[ref_gene]
        strongly_connected = set(k for k, v in ref_second_order.items() if v > 1)
        if len(strongly_connected) == 0:
            continue
        ret[ref_gene] = strongly_connected
    return ret


def sep_strong_conn(graph, strongly_connected):
    """Check whether any of strongly connected ref genes are separated."""
    if len(strongly_connected) == 0:
        # nothing to check
        return False
    graph_nodes = set(graph.nodes())
    s_conn_appear_here = {k: v for k, v in strongly_connected.items() if k in graph_nodes}
    if len(s_conn_appear_here) == 0:
        # nothing to check: no strongly connected genes here
        return False
    for k, v in s_conn_appear_here.items():
        conn_num = len(v)
        if len(v.intersection(graph_nodes)) < conn_num:
            # need this set to be included into graph nodes
            # if not the case: strongly connected nodes are in a different graph
            return True
    # nothing bad happened
    return False


def check_low_score_edges_removed(scores_edges_left, scores_edges_removed):
    """Return True if scores of removed edges are lower in comparison to remaining edges."""
    if all(x == 0.0 for x in scores_edges_removed):
        # in this case no orthology scores provided
        return True
    # for now a primitive method would be used
    # TODO: check whether this require any improvement
    min_score_left = min(scores_edges_left)
    max_score_rem = max(scores_edges_removed)
    min_thr = 0.90 * min_score_left
    if max_score_rem < min_thr:
        return True
    elif all(x < 0.75 for x in scores_edges_removed):
        return True
    # we are not certain -> difference is not significant enough
    return False


def split_graph(graph, r_genes, q_genes, edge_to_score):
    """Split many2many graph if possible."""
    # 1: find lead edges: edges connected to leaf nodes
    ref_nodes_set = set(r_genes)
    # query_nodes_set = set(q_genes)
    ref_conns, que_conns = edges_to_dicts(graph, ref_nodes_set)
    ref_leaf_edges = [(k, v[0]) for k, v in ref_conns.items() if len(v) == 1]
    que_leaf_edges = [(v[0], k) for k, v in que_conns.items() if len(v) == 1]
    leaf_edges = set(ref_leaf_edges + que_leaf_edges)
    if len(leaf_edges) == 0:
        # if no leaves: too complicated case, don't try to resolve
        return [graph, ]
    # 2: get nodes connected to leaf edges, and set difference
    leaf_conn_nodes = set(flatten(leaf_edges))
    not_leaf_conn_nodes = set(graph.nodes()).difference(leaf_conn_nodes)
    # 3: get list of reference nodes that must be not separated
    ref_not_separate = select_strongly_connected_ref_genes(ref_conns, que_conns)
    # 4: get local graph copy for manipulations
    _loc_graph_copy = nx.Graph()
    for edge in graph.edges():
        # probably is faster than deepcopy()
        _loc_graph_copy.add_edge(edge[0], edge[1])
    # 5: get part of the graph consisting only leaf+ nodes
    trimmed_parts_united = _loc_graph_copy.edge_subgraph(leaf_edges).copy()
    trimmed_parts = get_graph_components(trimmed_parts_united)
    remainder = _loc_graph_copy.subgraph(not_leaf_conn_nodes)
    # 6: check that remainder has no isolated nodes
    # if so: we don't try to resolve this many2many
    # return original graph
    if graph_has_isolated_points(remainder):
        return [graph, ]
    # add remainder only if it's not empty
    # otherwise it makes no sense
    if not graph_is_empty(remainder):
        trimmed_parts.append(remainder)
    # check whether we separated strongly connected genes
    if any(sep_strong_conn(g, ref_not_separate) for g in trimmed_parts):
        # we separated genes that we didn't want to: return initial graph
        return [graph, ]
    # 7: primitive "statistics"
    # Just check that scores of deleted edges are significantly lower than
    # scores of remaining edges
    original_edges = set(order_edges(graph.edges(), ref_nodes_set))
    edges_left = set(order_edges(flatten(x.edges() for x in trimmed_parts), ref_nodes_set))
    edges_removed = original_edges.difference(edges_left)
    edges_left_scores = [edge_to_score.get(e, 0.0) for e in edges_left]
    edges_rem_scores = [edge_to_score.get(e, 0.0) for e in edges_removed]
    bool__removed_low_score_edges = check_low_score_edges_removed(edges_left_scores,
                                                                  edges_rem_scores)
    # 8: return graph parts or the original graph depending on edge scores
    if bool__removed_low_score_edges is True:
        return trimmed_parts
    else:
        # here we are not certain, better to return original graph
        return [graph, ]


def get_graph_conn(graph, over__ref_genes, over__que_genes):
    """Create connection for subgraph."""
    # TODO: remove code duplicate
    graph_nodes = graph.nodes()
    ref_genes = [x for x in over__ref_genes if x in graph_nodes]
    que_genes = [x for x in over__que_genes if x in graph_nodes]
    ref_len, que_len = len(ref_genes), len(que_genes)
    c_class = get_c_class(ref_len, que_len)
    conn = {R_GENES: ref_genes, Q_GENES: que_genes, C_CLASS: c_class}
    return conn


def resolve_many2many(graph, r_genes, q_genes, edge_to_score):
    """Resolve many2many graph."""
    is_b_complete = is_complete_bipartite(graph, r_genes, q_genes)
    if is_b_complete:
        # nothing to do actually
        conn = {R_GENES: r_genes, Q_GENES: q_genes, C_CLASS: MANY2MANY}
        return [conn, ]
    # not complete bipartite graph
    # first, select reference genes that have a single connection
    # edges = [(node, node), (node, node), ..] -> list of pairs
    ret = []
    graph_parts = split_graph(graph, r_genes, q_genes, edge_to_score)
    for elem in graph_parts:
        conn = get_graph_conn(elem, r_genes, q_genes)
        ret.append(conn)
    return ret


def extract_orth_connections(graph, r_genes_all, q_genes_all, edge_to_score):
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
        c_class = get_c_class(r_len, q_len)

        if c_class == MANY2MANY:
            # this is many2many, a different procedure
            # maybe this is not a complete bipartite graph
            # then we should split it into sub graphs
            connections = resolve_many2many(component, r_genes, q_genes, edge_to_score)
            orth_connections.extend(connections)
            continue
        # not many2many: just save it
        # create connection object and save it
        conn = {R_GENES: r_genes, Q_GENES: q_genes, C_CLASS: c_class}
        orth_connections.append(conn)

    # count different orthology classes
    class_list = [c[C_CLASS] for c in orth_connections]
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
        conn_class = conn[C_CLASS]
        ref_genes = conn[R_GENES]
        que_genes = set(conn[Q_GENES])
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
                proj_not_added = True
                if not projections:
                    # probably, due to some reasons some transcripts don't have any
                    # associated orthology projection
                    # maybe it's lost (and gene has an intact isoform), or was skipped due to
                    # technical reasons
                    non_orthologous_isoforms.append(ref_transcript)
                for proj in projections:
                    # and last, each transcript can be projected > once
                    proj_q_gene = q_trans_to_gene[proj]
                    if proj_q_gene not in que_genes:
                        # there are orthologous projections detected earlier but
                        # we removed them earlier
                        continue
                    proj_not_added = False  # if True -> consider this transcript skipped
                    f.write(f"{ref_gene}\t{ref_transcript}\t{proj_q_gene}\t{proj}\t{conn_class}\n")
                if proj_not_added is True:
                    # see comments for (if not projections)
                    non_orthologous_isoforms.append(ref_transcript)
    # close file and return non-orthologous reference isoforms
    f.close() if out != "stdout" else None
    return non_orthologous_isoforms


def get_edge_score(ref_que_conn_scores):
    """Return ref_gene: que_genene: score dict."""
    ret = {}
    for edge, scores in ref_que_conn_scores.items():
        max_score = max(scores)
        ret[edge] = max_score
    return ret


def orthology_type_map(ref_bed, que_bed, out, ref_iso=None, que_iso=None,
                       paralogs_arg=None, loss_data=None, save_skipped=None,
                       orth_scores_arg=None):
    """Make orthology classification track."""
    q_trans_paralogs = read_paralogs(paralogs_arg)  # do not include paralogs
    q_trans_l_score = read_proj_scores(orth_scores_arg)
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
    r_gene_to_trans, r_trans_to_gene = read_isoforms__otm(ref_iso, ref_transcripts)
    q_gene_to_trans, q_trans_to_gene = read_isoforms__otm(que_iso, que_transcripts_all)
    r_genes_all = set(r_gene_to_trans.keys())
    q_genes_all = set(q_gene_to_trans.keys())
    # make transcript to projections dict:
    t_trans_to_projections = get_t_trans_to_projections(que_transcripts)
    # create graph to connect orthologous transcripts
    o_graph, ref_que_conn_scores = connect_genes(r_trans_to_gene,
                                                 t_trans_to_projections,
                                                 q_trans_to_gene,
                                                 q_trans_l_score)
    edge_to_score = get_edge_score(ref_que_conn_scores)
    # if a group of reference and query genes are in the same connected component
    # then they are orthologs
    orth_connections = extract_orth_connections(o_graph, r_genes_all, q_genes_all, edge_to_score)
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
    app.add_argument("--orth_scores", "-o", default=None, help="Orthology scores")
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
                       save_skipped=args.save_skipped,
                       orth_scores_arg=args.orth_scores)


if __name__ == "__main__":
    main()
