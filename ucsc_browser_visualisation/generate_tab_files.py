#!/usr/bin/env python3
"""Generate tab files from project results.

Later they could be converted to MySQL tables
for track info visualisation.
"""
import argparse
import sys
import os
from collections import defaultdict
from version import __version__

SPACE = "&nbsp;"
PLACE_HOLDER_EXON_MID = "".join([SPACE for _ in range(5)])
INACT_FEATS = [
    "INTACT_PERC_IGNORE_M",
    "INTACT_PERC_INTACT_M",
    "INTACT_CODONS_PROP",
    "OUT_OF_CHAIN_PROP",
    "MIDDLE_IS_INTACT",
    "MIDDLE_IS_PRESENT",
]
# for assembled from fragments: we cannot get chain class features
FRAGM_FEATS = (0, 0.0, 0.0, 0.0, 0.0, 0.0)

REF_BED = "toga_filt_ref_annot.bed"
QUE_BED = "query_annotation.bed"
ORTH_SCORES = "orthology_scores.tsv"
CHAIN_RESULTS_DF = "chain_results_df.tsv"
LOSS_SUMM_DATA = "loss_summ_data.tsv"
INACT_MUT_DATA = "inact_mut_data.txt"
CESAR_RESULTS = "cesar_results.txt"

CODON_FASTA = "codon.fasta"
PROT_FASTA = "prot.fasta"
TEMP = "temp"

ZERO_S = "0"
ONE_S = "1"
NINE_S = "9"  # for NA values, 0 - False, 1 - True, 9 - N/A


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def parse_args():
    """Parse command line args."""
    app = argparse.ArgumentParser()
    app.add_argument("wd", help="TOGA project dir containing results")
    # app.add_argument("ref_annot", help="Reference bed annotation.")
    # app.add_argument("ref_name", help="Reference name, example: hg38")
    # app.add_argument("que_name", help="Query name, example: mm10")
    app.add_argument("output", help="Output dir")
    app.add_argument(
        "--no_raw_cesar_output",
        "--nrco",
        dest="no_raw_cesar_output",
        action="store_true",
        help=("If raw cesar output is not present, extract "
              "data from default output fasta files. Will be "
              "the default behaviour in the future"),
    )
    app.add_argument(
        "--force_read_prot_fa",
        dest="force_read_prot_fa",
        action="store_true",
        help=("For very special cases only: extract protein data "
              "from prot.fasta instead of raw casar output (default). "
              "This arg is implemented to fix TOGAs in Feb2023")
    )
    app.add_argument("--no_plots",
                    "--np",
                    dest="no_plots",
                    action="store_true",
                    help="If inactivating mutation plots are already generated, do not recreate them")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def split_proj_name(proj_name):
    """Split projection name.

    Projections named as follows: ${transcript_ID}.{$chain_id}.
    This function splits projection back into transcript and chain ids.
    We cannot just use split("."), because there migth be dots
    in the original transcript ID.
    """
    proj_name_split = proj_name.split(".")
    q_num_str = proj_name_split[-1]
    trans_name = ".".join(proj_name_split[:-1])
    return trans_name, q_num_str


def get_query_coordinates(wd):
    """Make a dict with query coordinates."""
    print("Reading query annotation bed data")
    query_bed = os.path.join(wd, QUE_BED)
    f = open(query_bed, "r")
    proj_to_q_coords = {}
    for line in f:
        line_data = line.rstrip().split("\t")
        chrom = line_data[0]
        start = line_data[1]
        end = line_data[2]
        name = line_data[3]
        coords = f"{chrom}:{start}-{end}"
        proj_to_q_coords[name] = coords
    f.close()
    print(f"There are {len(proj_to_q_coords.keys())} projections")
    return proj_to_q_coords


def get_ref_coordinates(wd):
    """Read reference trascript coordinates."""
    print("Reading reference annotation data")
    ref_coords_file = os.path.join(wd, TEMP, REF_BED)
    f = open(ref_coords_file, "r")
    ref_trans_to_region = {}
    for line in f:
        line_data = line.rstrip().split("\t")
        chrom = line_data[0]
        start = line_data[1]
        end = line_data[2]
        region = f"{chrom}:{start}-{end}"
        name = line_data[3]
        ref_trans_to_region[name] = region
    f.close()
    return ref_trans_to_region


def get_chain_scores(wd, projections_list):
    """Extract chain orthology scores."""
    print(f"Reading chain scores")
    proj_to_chain_score = {}
    chain_scores_file = os.path.join(wd, TEMP, ORTH_SCORES)
    f = open(chain_scores_file, "r")
    f.__next__()
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        chain = line_data[1]
        score = line_data[2]
        proj_ = f"{trans}.{chain}"
        if proj_ not in projections_list:
            continue
        proj_to_chain_score[proj_] = score
    f.close()
    return proj_to_chain_score


def get_chain_features(wd, projections_list):
    """Extract chain features used for classification."""
    print("Extracting chain features")
    # get projection -> chain features
    proj_to_chain_features = {}
    f = open(os.path.join(wd, TEMP, CHAIN_RESULTS_DF), "r")
    f.__next__()
    # remembed the order
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        # gene_overs_ = line_data[1]
        chain = line_data[2]
        projection = f"{trans}.{chain}"
        if projection not in projections_list:
            continue
        synt_ = line_data[3]
        # gl_score_ = line_data[4]
        gl_exo_ = line_data[5]
        # chain_len_ = line_data[6]
        # cds_qlen_ = line_data[7]
        loc_exon_ = line_data[8]
        exon_cover_ = line_data[9]
        intron_cover_ = line_data[10]
        exon_fract_ = line_data[13]
        intron_fract_ = line_data[14]
        flank_cov_ = line_data[15]

        exon_cov = (
            str(float(exon_cover_) / float(exon_fract_))
            if float(exon_fract_) != 0
            else "0"
        )
        intron_cov = (
            str(float(intron_cover_) / float(intron_fract_))
            if float(intron_fract_) != 0
            else "0"
        )
        tup = (synt_, flank_cov_, gl_exo_, loc_exon_, exon_cov, intron_cov)
        proj_to_chain_features[projection] = tup
    f.close()
    return proj_to_chain_features


def get_projection_class(wd):
    """Get projection class per projection."""
    f = open(os.path.join(wd, LOSS_SUMM_DATA), "r")
    print("Reading gene loss classification")
    proj_to_classification = {}
    for line in f:
        line_data = line.rstrip().split("\t")
        line_class = line_data[0]
        if line_class != "PROJECTION":
            continue
        proj = line_data[1]
        lclass_raw = line_data[2]
        if lclass_raw == "I":
            lclass = "Intact"
        elif lclass_raw == "UL":
            lclass = "Uncertain Loss"
        elif lclass_raw == "PI":
            lclass = "Partially intact"
        elif lclass_raw == "M":
            lclass = "Missing"
        elif lclass_raw == "L":
            lclass = "Lost"
        elif lclass_raw == "PG":
            lclass = "Paralogous projection"
        elif lclass_raw == "PM":
            lclass = "Partial missing"
        else:
            lclass = "Undefined"
        proj_to_classification[proj] = lclass
    f.close()
    return proj_to_classification


def save_toga_info_tab(
    out_dir,
    projections_list,
    proj_to_q_coords,
    ref_trans_to_region,
    proj_to_chain_score,
    proj_to_chain_features,
    projection_to_loss_class,
):
    """Save tab file for TOGAinfo table."""
    toga_info_tab_path = os.path.join(out_dir, "togaInfo.tab")
    print(f"Saving TOGAInfo tab file / {len(projections_list)} items")
    f = open(toga_info_tab_path, "w")
    for projection in projections_list:
        trans, chain = split_proj_name(projection)
        glp_class = projection_to_loss_class.get(projection, "Missing")
        # default 0.5 for fragmented assemblies
        chain_score = proj_to_chain_score.get(projection, 0.5)
        query_region = proj_to_q_coords[projection]
        ref_region = ref_trans_to_region[trans]
        chain_feats = proj_to_chain_features.get(projection, FRAGM_FEATS)
        # parse chain ml features
        synteny = chain_feats[0]
        flank = chain_feats[1]
        gl_exo = chain_feats[2]
        loc_exo = chain_feats[3]
        exon_cov = chain_feats[4]
        intr_cov = chain_feats[5]
        tab_row = (
            projection,
            trans,
            ref_region,
            query_region,
            chain_score,
            synteny,
            flank,
            gl_exo,
            loc_exo,
            exon_cov,
            intr_cov,
            glp_class,
        )
        tab_strs_ = map(str, tab_row)
        f.write("\t".join(tab_strs_))
        f.write("\n")
    print(f"Saved togaInfo tab at {toga_info_tab_path}")
    f.close()


def format_as_ali(seq_1, seq_2, w=80):
    """Format sequences as alignment."""
    lines = ["<TT>"]
    seq_zip = list(zip(seq_1, seq_2))
    zip_parts = parts(seq_zip, w)
    for part in zip_parts:
        upper_seq = "".join([x[0] for x in part])
        lower_seq = "".join([x[1] for x in part])
        seq_len = len(upper_seq)
        middle = "".join(
            [SPACE if upper_seq[i] != lower_seq[i] else "|" for i in range(seq_len)]
        )
        lines.append(f"ref:{SPACE}{upper_seq}<BR>")
        lines.append(f"{PLACE_HOLDER_EXON_MID}{middle}<BR>")
        lines.append(f"que:{SPACE}{lower_seq}<BR><BR>")
    lines.append("</TT>")
    return "".join(lines)


def ret_cesar_lines(f):
    """Return relevant CESAR output lines."""
    for line in f:
        if line.startswith("#"):
            continue
        line_s = line.rstrip()
        if len(line_s) == 0:
            continue
        yield line_s
    yield None


def read_simple_pairwise_fasta(wd, fname, all_projections):
    """Read protein or codon fasta file."""
    prot_fasta = os.path.join(wd, fname)
    # for each protein/codon ali: ref and que
    proj_to_reference_seq = {}
    proj_to_query_seq = {}
    curr_projection = None
    is_ref = None
    skip_next = False

    f = open(prot_fasta, "r")
    for line in f:
        if line.startswith(">"):
            # this is header
            header_data = line.rstrip().lstrip(">").split(" | ")
            if header_data[1] != "PROT" and fname == PROT_FASTA:
                skip_next = True
                continue
            elif header_data[1] != "CODON" and fname == CODON_FASTA:
                skip_next = True
                continue
            curr_projection = header_data[0]
            is_ref = header_data[2] == "REFERENCE"
        else:
            # this is a sequence, must be related to the previously selected header
            if skip_next is True:
                skip_next = False
                continue
            seq = line.rstrip()
            if is_ref:
                proj_to_reference_seq[curr_projection] = seq
            else:
                proj_to_query_seq[curr_projection] = seq
            # for safety: reset curr proj and isref
            curr_projection = None
            is_ref = None
    f.close()

    projection_to_ali = {}
    for proj, que_prot_seq in proj_to_query_seq.items():
        ref_prot_seq = proj_to_reference_seq.get(proj, None)
        if not ref_prot_seq:
            print(proj)
        prot_ali = format_as_ali(ref_prot_seq, que_prot_seq)
        projection_to_ali[proj] = prot_ali
    return projection_to_ali


def extract_exons_meta_data(meta_data_file):
    """Just read exons meta data file."""
    exon_id_to_meta_features = {}
    f = open(meta_data_file, "r")
    f.__next__()
    """Header:
    gene    exon_num        chain_id        act_region      exp_region      in_exp  pid
    blosum  gap     class paralog  q_mark
    """
    for line in f:
        line_data = line.rstrip().split("\t")
        transcript_id = line_data[0]
        exon_num = int(line_data[1])
        chain_id = int(line_data[2])
        projection_id = f"{transcript_id}.{chain_id}"
        coordinates = line_data[3]
        exp_range = line_data[4]
        in_expected = ONE_S if line_data[5] == "INC" else ZERO_S
        pid = float(line_data[6])
        blosum = float(line_data[7])
        is_gap = ONE_S if line_data[8] == "GAP" else ZERO_S
        ali_class = line_data[9]
        start_end_strs = coordinates.split(":")[1].split("-")
        q_start = int(start_end_strs[0])
        q_end = int(start_end_strs[1])
        len_in_q = abs(q_end - q_start)
        
        exon_id = (projection_id, exon_num)
        exon_data = {
            "location": coordinates,
            "pid": pid,
            "blosum": blosum,
            "ali_class": ali_class,
            "is_gap": is_gap,
            "exp_range": exp_range,
            "in_exp": in_expected,
            "len_in_Q": len_in_q
        }
        exon_id_to_meta_features[exon_id] = exon_data
    f.close()
    return exon_id_to_meta_features


def extract_full_nucl_sequences(nucl_fasta):
    """Extract full projection nucleotide sequences."""
    p_to_ref = {}
    p_to_que = {}
    f = open(nucl_fasta, "r")
    current_proj = None
    reading_ref = None

    for line in f:
        if line.startswith(">"):
            # reading header
            header = line.lstrip(">").rstrip()
            if header.startswith("ref_"):
                reading_ref = True
                current_proj = header[len("ref_"):]
            else:
                reading_ref = False
                current_proj = header
            continue
        # reading sequence
        seq = line.rstrip()
        if reading_ref:
            p_to_ref[current_proj] = seq
        else:
            p_to_que[current_proj] = seq
        current_proj = None
        reading_ref = None            
    f.close()
    return p_to_ref, p_to_que


def __join_proj_to_nucl(p_to_ref, p_to_que):
    """Join projection to nucleotide dicts."""
    ret = {}
    for p in p_to_ref.keys():
        ref = p_to_ref[p]
        que = p_to_que[p]
        ret[p] = {"ref": ref, "que": que}
    return ret


def get_nucl_data_from_fasta(wd, all_projections):
    """Extract nucleotide data."""
    meta_data = os.path.join(wd, "temp", "exons_meta_data.tsv")
    nucl_fasta = os.path.join(wd, "nucleotide.fasta")
    exon_to_meta_data = extract_exons_meta_data(meta_data)
    projection_to_ref, projection_to_q = extract_full_nucl_sequences(nucl_fasta)
    projection_to_nucl = __join_proj_to_nucl(projection_to_ref, projection_to_q)
    for k, v in projection_to_nucl.items():
        print(k, v)
    exit()
    return []


def get_seq_data_from_fasta(wd, all_projections):
    """Get seq data from fasta files."""
    projection_exon_data = get_nucl_data_from_fasta(wd, all_projections)
    projection_to_prot_ali = read_simple_pairwise_fasta(wd, PROT_FASTA, all_projections)
    projection_to_codon_ali = read_simple_pairwise_fasta(wd, CODON_FASTA, all_projections)
    return [], projection_to_prot_ali, projection_to_codon_ali


def read_prot_fasta(wd):
    prot_fasta = os.path.join(wd, PROT_FASTA)
    # for each protein/codon ali: ref and que
    proj_to_reference_seq = {}
    proj_to_query_seq = {}
    curr_projection = None
    is_ref = None
    skip_next = False

    f = open(prot_fasta, "r")
    for line in f:
        if line.startswith(">"):
            # this is header
            header_data = line.rstrip().lstrip(">").split(" | ")
            if header_data[1] != "PROT":
                skip_next = True
                continue
            curr_projection = header_data[0]
            is_ref = header_data[2] == "REFERENCE"
        else:
            # this is a sequence, must be related to the previously selected header
            if skip_next is True:
                skip_next = False
                continue
            seq = line.rstrip()
            if is_ref:
                proj_to_reference_seq[curr_projection] = seq
            else:
                proj_to_query_seq[curr_projection] = seq
            # for safety: reset curr proj and isref
            curr_projection = None
            is_ref = None
    f.close()

    projection_to_ali = {}
    for proj, que_prot_seq in proj_to_query_seq.items():
        ref_prot_seq = proj_to_reference_seq.get(proj, None)
        if not ref_prot_seq:
            continue
        prot_ali = format_as_ali(ref_prot_seq, que_prot_seq)
        projection_to_ali[proj] = prot_ali
    return projection_to_ali


def get_sequence_data(wd, all_projections, exon_to_stat):
    """Parse nucleotide and protein sequence data."""
    print("Reading sequence data")
    gigafasta = os.path.join(wd, TEMP, CESAR_RESULTS)
    # all_transcripts = set(split_proj_name(x)[0] for x in all_projections)
    f = open(gigafasta, "r")
    lines_gen = ret_cesar_lines(f)
    query_exon_to_nucl_data = {}
    ref_exon_to_nucl_seq = {}
    # protein sequences collector
    projection_to_prot_seq = {}
    reference_to_prot_seq = {}
    projection_to_codon_seq = {}
    reference_to_codon_seq = {}

    while True:
        # reading lines two-by-two
        header = lines_gen.__next__()
        if not header:
            break
        sequence = lines_gen.__next__()
        if not sequence:
            break
        if not header.startswith(">"):
            raise ValueError("Error! Broken order of lines")
        # parse header and check what is this
        header_data = header[1:].replace(" ", "").split("|")
        # four options: ref/query AND prot/nucl
        if "PROT" in header_data:
            projection_id = header_data[0]
            if projection_id not in all_projections:
                continue
            # transcript_id, _ = split_proj_name(projection_id)
            # this is prot seq: query or reference
            if header_data[-1] == "REFERENCE":
                # reference seq
                reference_to_prot_seq[projection_id] = sequence
            else:
                # query seq
                projection_to_prot_seq[projection_id] = sequence
            continue
        elif "CODON" in header_data:
            # codon alignment
            # see comments in the previous (PROT) branch
            projection_id = header_data[0]
            if projection_id not in all_projections:
                continue
            if header_data[-1] == "REFERENCE":
                # reference seq
                reference_to_codon_seq[projection_id] = sequence
            else:
                # query seq
                projection_to_codon_seq[projection_id] = sequence
            continue
        else:
            # nucleotide seq, again ref or query
            # projection -> the same fields
            transcript_id = header_data[0]
            # initially exon num is 0-based and string
            # convert it to 1-based via int and get string back
            exon_num = str(int(header_data[1]) + 1)
            chain_id = header_data[2]
            projection_id = f"{transcript_id}.{chain_id}"
            if projection_id not in all_projections:
                continue
            # need this int for sorting later
            exon_id = (projection_id, int(exon_num))
            if header_data[-1] == "query_exon":
                # obviously, this is a header of query exon
                location = header_data[3]
                pid = header_data[4]
                blosum = header_data[5]
                is_gap = ONE_S if header_data[6] == "GAP" else ZERO_S
                ali_class = header_data[7]
                exp_range = header_data[8]
                in_exp = ONE_S if header_data[9] == "INC" else ZERO_S
                is_del_or_no = exon_to_stat.get(exon_id, "I")
                exon_data = (
                    location,
                    pid,
                    blosum,
                    is_gap,
                    ali_class,
                    exp_range,
                    in_exp,
                    is_del_or_no,  # added 26 Aug 2022
                    sequence,
                )
                query_exon_to_nucl_data[exon_id] = exon_data
            else:
                ref_exon_to_nucl_seq[exon_id] = sequence
            continue
    f.close()
    # save protein and codon alignments
    projection_to_prot_ali = {}
    for proj, que_prot_seq in projection_to_prot_seq.items():
        ref_prot_seq = reference_to_prot_seq[proj]
        prot_ali = format_as_ali(ref_prot_seq, que_prot_seq)
        projection_to_prot_ali[proj] = prot_ali
    projection_to_codon_ali = {}
    for proj, que_codon_seq in projection_to_codon_seq.items():
        ref_codon_seq = reference_to_codon_seq[proj]
        prot_ali = format_as_ali(ref_codon_seq, que_codon_seq)
        projection_to_codon_ali[proj] = prot_ali

    projection_exon_data = []
    # save exons data
    # sort exons by 1) projection ID and 2) their number
    sorted_exons = sorted(query_exon_to_nucl_data.keys(), key=lambda x: (x[0], x[1]))
    # for exon_id, exon_data in query_exon_to_nucl_data.items():
    for exon_id in sorted_exons:
        exon_data = query_exon_to_nucl_data[exon_id]
        projection_id, exon_num_int = exon_id
        # and here we need str for "".join operation
        exon_num = str(exon_num_int)
        ref_sequence = ref_exon_to_nucl_seq[exon_id]
        que_sequence = exon_data[-1]
        if len(ref_sequence) > 0 and que_sequence == "-":
            # meaning that query exon is completely deleted
            # but we need the same sequence length
            que_sequence = "".join("-" for _ in range(len(ref_sequence)))
        ali = format_as_ali(ref_sequence, que_sequence)
        rest = exon_data[:-1]
        fields = (projection_id, exon_num, *rest, ali)
        projection_exon_data.append(fields)
    return projection_exon_data, projection_to_prot_ali, projection_to_codon_ali


def save_toga_nucl_tab(out_dir, exon_data):
    """Saving exon data tab."""
    print("Saving exon data tab.")
    toga_nucl_tab_path = os.path.join(out_dir, "togaNucl.tab")
    f = open(toga_nucl_tab_path, "w")
    for elem in exon_data:
        f.write("\t".join(elem))
        f.write("\n")
    f.close()
    print(f"Saved TOGANucl tab to {toga_nucl_tab_path} / {len(exon_data)} items")


def save_toga_prot_tab(out_dir, prot_data):
    """Save prot data."""
    print("Saving protein sequences")
    toga_prot_tab_path = os.path.join(out_dir, "togaProt.tab")
    f = open(toga_prot_tab_path, "w")
    for proj, seq in prot_data.items():
        f.write(f"{proj}\t{seq}\n")
    f.close()
    print(f"Saved TOGAProt tab to {toga_prot_tab_path} / {len(prot_data)} items")


def get_inact_data(wd, all_projections):
    """Get everything related to inactivating mutations."""
    inact_mut_file = os.path.join(wd, INACT_MUT_DATA)
    projection_to_inact_muts = []
    proj_to_inact_features = defaultdict(dict)
    proj_to_inact_features_rows = []
    exon_to_del_miss = {}
    f = open(inact_mut_file, "r")
    for line in f:
        if not line.startswith("#"):
            continue
        line_data = line.rstrip().replace("# ", "").split("\t")
        trans = line_data[0]
        chain = line_data[1]
        projection = f"{trans}.{chain}"
        if projection not in all_projections:
            continue

        if len(line_data) == 8:
            # inactivating mutation data

            # IF Deleted exon OR Missing exon -> also save the data
            # needed in exon ali visualizations
            mask_field = line_data[6]
            is_inact = ZERO_S if mask_field == "masked" else ONE_S
            line_data[6] = is_inact
            mut_trimmed = line_data[5][:20]
            line_data[5] = mut_trimmed
            mut_track = line_data[2:]
            sql_row = (projection, *mut_track)
            projection_to_inact_muts.append(sql_row)
            if line_data[4] == "Deleted exon":
                exon_to_del_miss[(projection, int(line_data[2]))] = "D"
            elif line_data[4] == "Missing exon":
                exon_to_del_miss[(projection, int(line_data[2]))] = "M"
            continue
        else:
            # features
            feat, val_ = line_data[2].split()
            # boolean features to num (for SQL table)
            if feat == "MIDDLE_IS_PRESENT":
                if val_ == "TRUE":
                    val = ONE_S
                else:
                    val = ZERO_S
            elif feat == "MIDDLE_IS_INTACT":
                if val_ == "TRUE":
                    val = ONE_S
                else:
                    val = ZERO_S
            else:
                # numeric value
                val = val_
            proj_to_inact_features[projection][feat] = val
            continue
    f.close()

    for proj, feat_val_d in proj_to_inact_features.items():
        vals = []
        for f in INACT_FEATS:
            val = feat_val_d.get(f, 0.0)
            # 0.0 - means not found, but not for booleans
            # where 0 is False
            if val == 0.0 and f == "MIDDLE_IS_INTACT":
                val = NINE_S
            elif val == 0.0 and f == "MIDDLE_IS_PRESENT":
                val = NINE_S
            vals.append(val)
        sql_row = (proj, *vals)
        proj_to_inact_features_rows.append(sql_row)
    return proj_to_inact_features_rows, projection_to_inact_muts, exon_to_del_miss


def save_toga_inact_mut_tab(out_dir, inact_mut_data):
    """Save togaInactMut table."""
    print("Saving toga inact mut data")
    toga_inact_mut_path = os.path.join(out_dir, "togaInactMut.tab")
    f = open(toga_inact_mut_path, "w")
    for elem in inact_mut_data:
        f.write("\t".join(elem))
        f.write("\n")
    f.close()
    print("Saved")


def save_toga_inact_feat_tab(out_dir, inact_feat_data):
    """Save togaInactMut table."""
    print("Saving toga feat mut data")
    toga_inact_feat_path = os.path.join(out_dir, "togaInactFeat.tab")
    f = open(toga_inact_feat_path, "w")
    for elem in inact_feat_data:
        f.write("\t".join(elem))
        f.write("\n")
    f.close()
    print("Saved")


def main():
    """Entry point."""
    args = parse_args()
    os.mkdir(args.output) if not os.path.isdir(args.output) else None
    proj_to_q_coords = get_query_coordinates(args.wd)
    all_projections = set(proj_to_q_coords.keys())
    trans_to_ref_coords = get_ref_coordinates(args.wd)
    proj_to_chain_score = get_chain_scores(args.wd, all_projections)
    proj_to_chain_features = get_chain_features(args.wd, all_projections)
    projection_to_loss_class = get_projection_class(args.wd)
    proj_to_inact_feat, proj_to_inact_mut, exon_to_stat = get_inact_data(args.wd, all_projections)

    if args.no_raw_cesar_output:
        raise NotImplementedError("Cancelled branch")
        seq_data = get_seq_data_from_fasta(args.wd, all_projections)
        exit()
    else:
        seq_data = get_sequence_data(args.wd, all_projections, exon_to_stat)
    projection_exon_data = seq_data[0]
    projection_to_prot_ali = seq_data[1]
    projection_to_codon_ali = seq_data[2]

    if args.force_read_prot_fa:
        # override default protein sequences source
        print("Force reading protein alignments from prot.fasta")
        projection_to_prot_ali = read_prot_fasta(args.wd)
    save_toga_info_tab(
        args.output,
        all_projections,
        proj_to_q_coords,
        trans_to_ref_coords,
        proj_to_chain_score,
        proj_to_chain_features,
        projection_to_loss_class,
    )
    save_toga_nucl_tab(args.output, projection_exon_data)
    save_toga_prot_tab(args.output, projection_to_prot_ali)
    save_toga_inact_feat_tab(args.output, proj_to_inact_feat)
    save_toga_inact_mut_tab(args.output, proj_to_inact_mut)


if __name__ == "__main__":
    main()
