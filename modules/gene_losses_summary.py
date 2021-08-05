#!/usr/bin/env python3
"""Gene losses summary.

Classify projections, transcripts and genes
as lost, uncertain, missing etc.
"""
import argparse
import sys
import os
from collections import defaultdict
from collections import Counter
from datetime import datetime as dt

try:
    from modules.common import make_cds_track
    from modules.common import split_proj_name
    from modules.common import die
    from modules.common import eprint
    from modules.common import read_isoforms_file
    from modules.GLP_values import *
except ImportError:
    from common import make_cds_track
    from common import split_proj_name
    from common import die
    from common import eprint
    from common import read_isoforms_file
    from GLP_values import *


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# GLP classes
# kind of enum
N_ = -1  # No data at all
PG = 0  # ParaloG
PM = 1  # Partial missing
# L = 2  # Lost
# M = 3  # Missing
# Updated order:
M = 2
L = 3
UL = 4  # Uncertain loss
PI = 5  # Partially intact
I = 6  # Intact
# N - skipped due to technical reasons
# NUM_TO_CLASS = {-1: "N", 0: "PG", 1: "PM", 2: "L", 3: "M", 4: "UL", 5: "PI", 6: "I"}
# updated order
NUM_TO_CLASS = {-1: "N", 0: "PG", 1: "PM", 2: "M", 3: "L", 4: "UL", 5: "PI", 6: "I"}
CLASS_TO_NUM = {v: k for k, v in NUM_TO_CLASS.items()}


# link GLP class to color
CLASS_TO_COL = {
    N_: BLACK,
    PG: BROWN,
    PM: GREY,
    L: LIGHT_RED,
    M: GREY,
    UL: SALMON,
    PI: LIGHT_BLUE,
    I: BLUE,
}


REM_T_L = 0.35  # less than REM_T_L of CDS left -> it's lost
REM_T_G = 0.49  # less than REM_T_G of CDS left -> it's UL
PART_THR = 0.5  # border between missing and partially intact

PROJECTION = "PROJECTION"
TRANSCRIPT = "TRANSCRIPT"


def parse_args():
    """Read and parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("loss_data", help="Directory containing loss checker output files")
    app.add_argument("ref_bed", help="Reference bed file.")
    app.add_argument("pre_final_bed", help="Pre-final query annotation track")
    app.add_argument("bed_out", help="Bed file output...")
    app.add_argument("summary", help="Save summary to...")
    app.add_argument("--isoforms", "-i", help="Provide isoforms data to classify genes")
    app.add_argument(
        "--trace", "-t", default=None, help="Trace a particular isoform fate"
    )
    app.add_argument(
        "--paral_projections",
        default=None,
        help="File containing paralogous projections",
    )
    app.add_argument("--exclude", default=None, help="List of transcripts to exclude")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_loss_data(loss_dir):
    """Read inact mutations data for each projection.

    Projection is a predicted transcript in the query.
    TOGA gets a projection when projects a transcript via a chain.
    We parse two sorts of information associated with each transcript.:
    1) There are 6 features such as %intact.
    2) A list of inactivating mutations (could be empty).
    """
    # initiate dictionaries for features / mutations:
    projection_to_p_intact_M_ignore = {}
    projection_to_p_intact_M_intact = {}
    projection_to_i_codon_prop = {}
    proj_to_prop_oub = {}
    proj_to_80_p_intact = {}
    proj_to_80_p_present = {}
    projection_to_mutations = defaultdict(list)
    loss_files = os.listdir(loss_dir)

    for l_file in loss_files:
        # go file-by-file; because CESAR jobs produce a number of files
        path = os.path.join(loss_dir, l_file)
        f = open(path, "r")
        for line in f:
            # then line-by-line
            if not line.startswith("#"):
                # mutations-related lines should start with #
                continue
            # parse inact mutation data
            # [2:] to cut "# "
            line_data = line[2:].rstrip().split("\t")
            transcript_id = line_data[0]
            query_name = line_data[1]  # synonym for chain_id
            projection_id = f"{transcript_id}.{query_name}"

            # a section of %intact-related features
            if line_data[2].startswith("INTACT_PERC_IGNORE_M"):
                # intact percent branch; ignore missing sequence mode
                perc = float(line_data[2].split()[1])
                projection_to_p_intact_M_ignore[projection_id] = perc
                continue
            elif line_data[2].startswith("INTACT_PERC_INTACT_M"):
                # intact percent branch; consider missing part as intact
                perc = float(line_data[2].split()[1])
                projection_to_p_intact_M_intact[projection_id] = perc
                continue
            elif line_data[2].startswith("MIDDLE_80%_INTACT"):
                # flag: are there inact mutations in the middle 80% of CDS?
                raw_val = line_data[2].split()[1]
                val = True if raw_val == "TRUE" else False
                proj_to_80_p_intact[projection_id] = val
                continue
            elif line_data[2].startswith("MIDDLE_80%_PRESENT"):
                # flag: any missing fragment in the middle 80% if CDS?
                raw_val = line_data[2].split()[1]
                val = True if raw_val == "TRUE" else False
                proj_to_80_p_present[projection_id] = val
                continue
            elif line_data[2].startswith("INTACT_CODONS_PROP"):
                # proportion of intact codons
                # codons that are not deleted, missing and have no inact mutations
                perc = float(line_data[2].split()[1])
                projection_to_i_codon_prop[projection_id] = perc
                continue
            elif line_data[2].startswith("OUT_OF_CHAIN_PROP"):
                # proportion of transcript that lies beyond the chain
                perc = float(line_data[2].split()[1])
                proj_to_prop_oub[projection_id] = perc
                continue

            # a section of inactivating mutations
            exon_num = int(line_data[2])
            codon_num = line_data[3]
            mut_class = line_data[4]
            mut_itself = line_data[5]
            masked = True if line_data[6] == "masked" else False
            mut_ = (exon_num, codon_num, mut_class, mut_itself, masked)
            projection_to_mutations[projection_id].append(mut_)
        f.close()
    # to avoid returning a bunch of variables I packed them into a tuple
    output = (
        projection_to_mutations,
        projection_to_p_intact_M_ignore,
        projection_to_p_intact_M_intact,
        projection_to_i_codon_prop,
        proj_to_prop_oub,
        proj_to_80_p_intact,
        proj_to_80_p_present,
    )
    return output


def read_bed(bed_file):
    """Read bed file."""
    gene_line = {}
    f = open(bed_file, "r")
    for line in f:
        gene_name = line.split("\t")[3]
        gene_line[gene_name] = line
    f.close()
    return gene_line


def get_l_exon_num(exon_num):
    """20% of exons must be affected to lost the gene."""
    if exon_num == 1:
        return 1
    elif exon_num <= 10:
        return 2
    else:
        twenty_perc = exon_num / 5
        return twenty_perc


def get_projection_classes(
    all_projections,
    trans_exon_sizes,
    p_to_pint_m_ign,
    p_to_pint_m_int,
    projection_to_mutations,
    p_to_i_codon_prop,
    p_to_p_out_of_bord,
    p_80_int,
    p_80_pre,
    trace=None,
    paral_=None,
):
    """Classify projections as intact, lost, uncertain, etc."""
    projection_class = {}  # our answer: projection -> class
    # deal with paral_ argument
    if paral_ is None:
        # means that ids of paralogous projections are not provided
        # = set() because we check each projection whether it appears
        # in the paral set
        paral = set()
    else:
        # if is defined: it's already a set
        paral = paral_

    for num, projection in enumerate(all_projections):
        if projection in paral:
            # paralogous projection -> separate class automatically
            projection_class[projection] = PG
            continue
        # unpack all features for this projection
        p_intact_M_ign = p_to_pint_m_ign.get(projection, -1)
        p_intact_M_int = p_to_pint_m_int.get(projection, -1)
        p_i_codons = p_to_i_codon_prop.get(projection, -1)
        no_loss_in_80_p = p_80_int.get(projection, None)
        m_80_present = p_80_pre.get(projection, None)
        frame_oub = p_to_p_out_of_bord.get(projection, 0.0)
        transcript, _ = split_proj_name(projection)
        # tracing_: works if called as a standalone script
        # if set, script writes additional information about
        # the traced transcript -> like what decision made and why
        tracing_ = trace and transcript == trace

        exon_sizes = trans_exon_sizes.get(transcript)
        if exon_sizes is None:
            # sanity check, this must never happen
            print(
                f"Cannot find transcript {transcript}; probably an error"
            ) if tracing_ else None
            projection_class[projection] = N_
            continue

        # parse inactivating mutations
        all_mutations = projection_to_mutations.get(projection, [])
        # get only inactivating mutations, which are not compensations and not masked
        mutations = [
            m for m in all_mutations if m[4] is False or m[2] == MISS_EXON
        ]  # m[4]: bool MASKED
        mutations = [m for m in mutations if m[2] != COMPENSATION]

        if tracing_:
            # verbosity
            print(f"Projection: {projection}")
            print(f"%intact_Mign: {p_intact_M_ign} | %intact_Mint: {p_intact_M_int}")
            print(f"No int in m80%: {no_loss_in_80_p} | no miss in m80: {m_80_present}")
            print(f"Overall events: {len(all_mutations)}")
            print(f"Inact mutations: {len(mutations)}")
            print(f"Prop intact codons: {p_i_codons}")
            print(f"Exon sizes:\n{exon_sizes}")

        if len(mutations) == 0 and p_intact_M_ign > 0.6:
            # intact, >60% CDS is presented
            if tracing_:
                print("No mutations, p_M_ign > 0.6: Intact")
            # nothing to talk about
            projection_class[projection] = I
            continue
        elif p_i_codons < REM_T_L:
            # too small fraction of CDS is presented
            if tracing_:
                print(f"In this projection only {p_i_codons} codons remain intact")
                print(f"Threshold is {REM_T_L} codons")
                print(f"-> class L")
            projection_class[projection] = L
            continue

        # initiate exon -> status dict, exon could be intact, deleted, missing or have inact mut
        exon_status = {k: "I" for k in exon_sizes}
        exon_num = len(exon_sizes.keys())
        overall_seq_len = sum(exon_sizes.values())
        # get list of exons that hold > 40% of CDS
        exon_40_p = {k: v / overall_seq_len > 0.4 for k, v in exon_sizes.items()}
        exon_40_p_nums = [k for k, v in exon_40_p.items() if v is True]

        # select missing and deleted exons
        missing_exons = [m[0] for m in mutations if m[2] == MISS_EXON]
        deleted_exons = [m[0] for m in mutations if m[2] == DEL_EXON]
        # update status for deleted and missing exons:
        for me in missing_exons:
            exon_status[me] = "M"
        for de in deleted_exons:
            exon_status[de] = "D"

        print(f"Exon statuses initial\n{exon_status}") if tracing_ else None
        # deal with smaller mutations (such as frameshifts, stop codons, and so on)
        other_muts = [
            m
            for m in mutations
            if m[2] != MISS_EXON and m[2] != DEL_EXON and m[2] != COMPENSATION
        ]
        # don't need events happened in the deleted/missing exons
        other_muts = [
            m
            for m in other_muts
            if m[0] not in deleted_exons and m[0] not in missing_exons
        ]
        # also compute % of missing sequence
        missed_seq_len = sum(exon_sizes[k] for k, v in exon_status.items() if v == "M")
        missing_prop = missed_seq_len / overall_seq_len
        if tracing_:
            print(f"% Missing: {missing_prop}")
            print(
                f"Missing seq len: {missed_seq_len}; overall seq len: {overall_seq_len}"
            )
        # compute number of exons that we require to have inact mutations to call the gene Lost:
        affected_thr = get_l_exon_num(exon_num)

        # special cases: all exons are missed or deleted
        if all(v == "M" for v in exon_status.values()):
            print("All exons are Missing: transcript is missing") if tracing_ else None
            print(f"Out of borders prop is {frame_oub}") if tracing_ else None
            if frame_oub > PART_THR:
                print(
                    "-> class PM, too big fraction ouf of chain borders"
                ) if tracing_ else None
                projection_class[projection] = PM
            else:
                projection_class[projection] = M
            continue
        elif all(v == "D" for v in exon_status.values()):
            print("All exons are Deleted -> class L") if tracing_ else None
            projection_class[projection] = L
            continue

        if p_intact_M_int < 0.2:
            # < 20% intact, this gene is clearly lost
            print("P_intact_M_int < 20%: class L") if tracing_ else None
            projection_class[projection] = L
            continue

        # the main classification process
        if no_loss_in_80_p is True:
            print("GO TO BRANCH 1: No inact mut in m80%") if tracing_ else None
            # first branch -> no inact mutations in the middle 80% of CDS
            # possible: Intact, partially intact, missed
            # uncertain if too small fraction is left
            # if there are any inact mutations -> it doesn't really matter
            if p_i_codons < REM_T_G:
                # %intact codons is less than UL threshold: it's UL
                # if it is also less that lost threshold: it was classified above
                print(f"Prop of intact codons is {p_i_codons}") if tracing_ else None
                print(f"Need > {REM_T_G} to be I/PI") if tracing_ else None
                print(f"-> Class Uncertain Loss") if tracing_ else None
                projection_class[projection] = UL
                continue
            if len(missing_exons) == 0:
                # everything is fine, no missing sequence at all, middle 80% intact
                print(
                    "No missing exons, no inact mut in m80%: class I"
                ) if tracing_ else None
                projection_class[projection] = I
                continue

            # there are missing sequence, need to check carefully
            if m_80_present:
                # means that there are no missing sequence in the middle 80%
                print("No M exons in m80%: class I") if tracing_ else None
                projection_class[projection] = I
                continue
            else:
                # there are missing sequence in the middle 80%
                # need a %missing value to decide, M or PI class (or even PM)
                print("There are M in m80%...") if tracing_ else None
                # partially intact or missed
                if missing_prop < 0.5:
                    print("Missing prop < 50%: class PI") if tracing_ else None
                    # > 50% CDS is presented
                    projection_class[projection] = PI
                    continue
                else:
                    print("Missing prop > 50%: PM or M branch") if tracing_ else None
                    # > 50% of CDS missed, just call it missed
                    if frame_oub > PART_THR:
                        print(
                            "-> class PM, too big fraction out of chain borders"
                        ) if tracing_ else None
                        projection_class[projection] = PM
                    else:
                        projection_class[projection] = M
                    continue

        else:  # second major branch
            print("GO TO BRANCH 2: there are inact mut in m80%") if tracing_ else None
            # second major branch: there ARE inact mutations in the middle 80% of CDS
            # possible classes: uncertain loss or lost
            for m in other_muts:  # update exon status
                exon_status[m[0]] = "L"
            print(f"Exon status is:\n{exon_status}") if tracing_ else None

            # different decision tree for single and multi-exon genes
            if exon_num == 1:
                # special case, for a single exon gene we require > 2 mutations
                print("Single exon branch") if tracing_ else None
                if exon_status[1] == "D":
                    print("Single exon deleted: class L") if tracing_ else None
                    # the only exon is deleted -> lost
                    projection_class[projection] = L
                elif exon_status[1] == "M":
                    print("Single exon missing: class M") if tracing_ else None
                    # the only exon is missed -> M (redundant branch highly likely)
                    if frame_oub > PART_THR:
                        print(
                            "-> class PM, too big fraction out of chain borders"
                        ) if tracing_ else None
                        projection_class[projection] = PM
                    else:
                        projection_class[projection] = M
                elif p_intact_M_ign < 0.6 and len(other_muts) >= 2:
                    # less than 60% intact and >= 2 inact mut: Lost according to our methodology
                    print("%intact < 60% && 2 incat mut: class L") if tracing_ else None
                    projection_class[projection] = L
                else:
                    # there are inact mutations in the middle 80% of CDS
                    # not enough evidence to say it's lost: Uncertain
                    print("Not enough evidence for L -> UL") if tracing_ else None
                    projection_class[projection] = UL
                continue
            # multi-exon gene branch
            num_exons_affected = len(
                [k for k, v in exon_status.items() if v == "D" or v == "L"]
            )
            # the first question: is %intact < 60%?
            if p_intact_M_int < 0.6:
                print(f"% intact M int < 60 branch") if tracing_ else None
                # well, %intact < 60, maybe a loss!
                if tracing_:
                    print(
                        f"Affected exons: {num_exons_affected}; required: {affected_thr}"
                    )
                if num_exons_affected >= affected_thr:
                    # check number of affected exons (with inactivating mutations)
                    print(f"Enough affected exons -> L") if tracing_ else None
                    projection_class[projection] = L
                    continue

                # also check whether there is an exon covering > 40% CDS that has TWO inact mutations
                muts_occur = Counter(m[0] for m in other_muts)
                muts_in_40_exons = [muts_occur[x] for x in exon_40_p_nums]
                if any(x >= 2 for x in muts_in_40_exons):
                    # if there is an exon that takes >40% of CDS and has >= 2 inact muts: it's lost
                    print(
                        f"There are exons > 40% with 2+ mutations -> L"
                    ) if tracing_ else None
                    projection_class[projection] = L
                    continue
                if any(exon_status[x] == "D" for x in exon_40_p_nums):
                    # also: if any of significant exons is deleted
                    print(
                        f"Some of exons > 40% are Deleted -> Lost"
                    ) if tracing_ else None
                    projection_class[projection] = L
                    continue
                print(
                    f"Not enough evidence for Lost -> Uncertain"
                ) if tracing_ else None
                projection_class[projection] = UL
            else:
                # if %intact > 60: cannot be intact
                # but not enough evidence to call it lost
                # class uncertain loss OR missing
                print(f"% intact M int > 60 branch") if tracing_ else None
                if frame_oub > PART_THR:
                    print(
                        "-> class PM, too big fraction out of chain borders"
                    ) if tracing_ else None
                    projection_class[projection] = PM
                else:
                    print(
                        f"not enough evidence for L -> Uncertain"
                    ) if tracing_ else None
                    projection_class[projection] = UL
                continue
    return projection_class


def get_exon_sizes(ref_bed):
    """Get exon: size dict for reference bed."""
    trans_exon_sizes = {}
    f = open(ref_bed, "r")
    for line in f:
        # first of all, we need only coding exons
        # so I apply make_cds_track to remove UTRs
        cds_line = make_cds_track(line)
        # parse the line according to bed12 specification
        line_data = cds_line.split("\t")
        strand = line_data[5]
        exon_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
        if strand == "-":
            # if strand is -: exon sizes ordered from last to the first
            # need to reverse the list
            exon_sizes = exon_sizes[::-1]
        # create exon_num: size dict and save
        exon_to_size = {n: v for n, v in enumerate(exon_sizes, 1)}
        # make_cds_track adds _CDS substring to the gene name
        trans_id = line_data[3].replace("_CDS", "")
        trans_exon_sizes[trans_id] = exon_to_size
    f.close()
    return trans_exon_sizes


def get_paralogs_data(paral_file):
    """Extract paralogous projections."""
    # basically just read a file and save to a set
    if paral_file is None:
        return set()
    with open(paral_file, "r") as f:
        paral_proj = set(x.rstrip() for x in f.readlines())
    return paral_proj


def color_bed_file(bed_in, bed_out, proj_to_class):
    """Assing colors to bed tracks according to projection status."""
    in_ = open(bed_in, "r")
    out_ = open(bed_out, "w")
    for line in in_:
        if line == "\n":
            # for some reason sometimes we get empty lines, get rid of them
            continue
        line_data = line.rstrip().split("\t")
        if line_data[0] == "None":
            # force skip None chrom bed tracks
            # for some reason, they can reach this point
            continue
        # take projection ID from the line and then it's class
        projection_id = line_data[3]
        projection_class = proj_to_class.get(projection_id, N_)
        if projection_class == N_:
            eprint(f"Warning! Class of {projection_id} not found")
        color = CLASS_TO_COL[projection_class]
        # write corresponding color to file
        line_data[8] = color
        line_upd = "\t".join(line_data)
        out_.write(line_upd)
        out_.write("\n")
    # close files, that's it
    in_.close()
    out_.close()


def remove_unused_trans(gene_to_trans, iforms_to_save):
    """Remove unnecessary transcripts from the isoforms file."""
    ret = defaultdict(list)
    for g, trans in gene_to_trans.items():
        for t in trans:
            if t not in iforms_to_save:
                continue
            ret[g].append(t)
    return ret


def read_excluded(exc_arg):
    """Read the list of excluded genes."""
    if exc_arg is None:
        return set()
    f = open(exc_arg, "r")
    ret = set(x.rstrip() for x in f)
    f.close()
    return ret


def read_predefined_glp_data(predef_class):
    """Read predefined classfications."""
    trans_class = {}
    proj_class = {}
    if predef_class is None:
        return proj_class, trans_class
    for item in predef_class:
        entry_class = item[0]
        entry = item[1]
        entry_status = item[2]
        if entry_class == TRANSCRIPT:
            trans_class[entry] = CLASS_TO_NUM[entry_status]
        else:
            proj_class[entry] = CLASS_TO_NUM[entry_status]
    return proj_class, trans_class


def gene_losses_summary(
    loss_data_arg,
    ref_bed,
    pre_final_bed_arg,
    bed_out,
    summary_arg,
    trace_arg=None,
    iforms_file=None,
    paral=None,
    exclude_arg=None,
    predefined_class=None,
):
    """Gene losses summary core function."""
    t0 = dt.now()
    # TOGA don't make any conclusions about projections via paralogous chains
    paralogs_set = get_paralogs_data(paral)
    excluded_genes = read_excluded(exclude_arg)
    # we need transcript exons sizes for decision tree
    trans_exon_sizes = get_exon_sizes(ref_bed)
    # parse inactivating mutations data
    loss_data_all = read_loss_data(loss_data_arg)
    # unpack the returned tuple:
    projection_to_mutations = loss_data_all[0]
    p_to_pintact_M_ign = loss_data_all[1]
    p_to_pintact_M_int = loss_data_all[2]
    p_to_i_codon_prop = loss_data_all[3]
    p_to_oub_prop = loss_data_all[4]
    p_80_int = loss_data_all[5]
    p_80_pre = loss_data_all[6]
    # read predefined glp classes if they are:
    predef_proj_class, predef_trans_class = read_predefined_glp_data(predefined_class)

    # for consistency: get a set of all possible projections
    all_projections = set(p_to_pintact_M_ign.keys()).union(
        set(projection_to_mutations.keys())
    )
    # call this function to classify projections
    projection_class = get_projection_classes(
        all_projections,
        trans_exon_sizes,
        p_to_pintact_M_ign,
        p_to_pintact_M_int,
        projection_to_mutations,
        p_to_i_codon_prop,
        p_to_oub_prop,
        p_80_int,
        p_80_pre,
        trace=trace_arg,
        paral_=paralogs_set,
    )
    projection_class.update(predef_proj_class)
    # add projections that we added from predefined list
    all_projections.update(projection_class.keys())

    # projections are classified, we can color the bed file now:
    color_bed_file(pre_final_bed_arg, bed_out, projection_class)

    # get transcript -> [projections] dict
    # one transcript might have > 1 orthologous chain:
    # so it might have > 1 projection
    transcript_to_projections = defaultdict(list)
    for proj in all_projections:
        # projection is: $transcript DOT $chain_id
        # split it to trans and chain (we don't need chain here)
        trans, _ = split_proj_name(proj)
        if trans in excluded_genes:
            continue
        transcript_to_projections[trans].append(proj)

    # classify transcripts
    transcript_class = {}
    for trans, projections in transcript_to_projections.items():
        # idea is the following: we get classes of all projections
        # associated with this transcript
        # then we get the "best" class and assign it to this transcript
        # it a transcript has M, L and UL projections, we will
        # classify this transcript as Uncertain Loss
        # classes are "enum" so we can just apply max() function
        p_classes = set(projection_class.get(p) for p in projections)
        status = max(p_classes)  # just use Enum I > UL > L > M > N
        transcript_class[trans] = status
    transcript_class.update(predef_trans_class)
    all_transcripts = set(transcript_class.keys())

    # classify genes part
    if iforms_file:  # if isoforms provided, get gene: [transcripts] dict
        gene_to_trans, _, _ = read_isoforms_file(
            iforms_file, pre_def_trans_list=all_transcripts
        )
    else:  # no isoforms provided: nothing we can do next
        gene_to_trans = {}

    gene_class = {}
    for gene, transcripts in gene_to_trans.items():
        # classify genes applying the same logic we used to classify transcripts
        trans_statuses = [transcript_class.get(t, -1) for t in transcripts]
        status = max(trans_statuses)
        gene_class[gene] = status

    # save results, create 3 column table
    # 1st column: what is classified: projection, transcript or a gene
    # 2nd column: item identifier (like gene name)
    # 3rd: the class itself, like I (intact)
    f = open(summary_arg, "w")
    # first of all save projections classification
    for k, v in projection_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"PROJECTION\t{k}\t{v_ch}\n")
    # then save transcripts classification
    for k, v in transcript_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"TRANSCRIPT\t{k}\t{v_ch}\n")
    # at last genes classification
    for k, v in gene_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"GENE\t{k}\t{v_ch}\n")
    f.close()
    print(f"Elapsed: {dt.now() - t0}")


def main():
    """Entry point for CLI."""
    args = parse_args()
    gene_losses_summary(
        args.loss_data,
        args.ref_bed,
        args.pre_final_bed,
        args.bed_out,
        args.summary,
        trace_arg=args.trace,
        iforms_file=args.isoforms,
        paral=args.paral_projections,
        exclude_arg=args.exclude,
    )


if __name__ == "__main__":
    main()
