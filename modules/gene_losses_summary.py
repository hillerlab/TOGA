#!/usr/bin/env python3
"""Gene losses summary."""
import argparse
import sys
import os
from collections import defaultdict
from collections import Counter
from datetime import datetime as dt
try:
    from modules.common import make_cds_track
    from modules.common import split_proj_name
except ImportError:
    from common import make_cds_track
    from common import split_proj_name

# classes
# kinda enum
N_ = -1
PG = 0  # ParaloG
PM = 1  # Partial missing
L = 2  # Lost
M = 3  # Missing
G = 4  # Grey
PI = 5  # Partially intact
I = 6  # Intact
# N - skipped due to technical reasons
NUM_TO_CLASS = {-1: "N", 0: "PG", 1: "PM", 2: "L", 3: "M", 4: "G", 5: "PI", 6: "I"}

# colors
BLUE = "0,0,200"
LIGHT_BLUE = "0,200,255"
LIGHT_RED = "255,50,50"
SALMON = "255,160,120"
GREY = "130,130,130"
BROWN = "159,129,112"
BLACK = "10,10,10"

CLASS_TO_COL = {N_: BLACK,  PG: BROWN, PM: GREY, L: LIGHT_RED,
                M: GREY, G: SALMON, PI: LIGHT_BLUE, I: BLUE}

# mut classes:
MISS_EXON = "Missing exon"
DEL_EXON = "Deleted exon"
DEL_MISS = {MISS_EXON, DEL_EXON}
COMPENSATION = "COMPENSATION"
SSM = "SSM"
START_MISSING = "START_MISSING"
FS_DEL = "FS_DEL"
FS_INS = "FS_INS"
BIG_DEL = "BIG_DEL"
BIG_INS = "BIG_INS"
STOP = "STOP"

REM_T_L = 0.35
REM_T_G = 0.49
PART_THR = 0.5


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(str(msg) + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_args():
    """Read and parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("loss_data", help="Directory containing loss checker output files")
    app.add_argument("ref_bed", help="Reference bed file.")
    app.add_argument("pre_final_bed", help="Pre-final query annotation track")
    app.add_argument("bed_out", help="Bed file output...")
    app.add_argument("summary", help="Save summary to...")
    app.add_argument("--isoforms", "-i", help="Provide isoforms data to classify genes")
    app.add_argument("--trace", "-t", default=None, help="Trace a particular isoform fate")
    app.add_argument("--paral_projections", default=None, help="File containing paralogous projections")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_loss_data(loss_dir):
    """Read loss data for each transcript."""
    projection_to_p_intact_M_ignore = {}
    projection_to_p_intact_M_intact = {}
    projection_to_i_codon_prop = {}
    proj_to_prop_oub = {}
    proj_to_80_p_intact = {}
    proj_to_80_p_present = {}
    projection_to_mutations = defaultdict(list)
    loss_files = os.listdir(loss_dir)

    for lfile in loss_files:
        # go file-by-file
        path = os.path.join(loss_dir, lfile)
        f = open(path, "r")
        for line in f:
            if not line.startswith("#"):
                # how I encode mutation-related lines
                continue
            line_data = line[2:].rstrip().split("\t")
            transcript_id = line_data[0]
            query_name = line_data[1]
            projection_id = f"{transcript_id}.{query_name}"

            if line_data[2].startswith("INTACT_PERC_IGNORE_M"):
                # intact percent branch
                perc = float(line_data[2].split()[1])
                projection_to_p_intact_M_ignore[projection_id] = perc
                continue
            elif line_data[2].startswith("INTACT_PERC_INTACT_M"):
                # intact percent branch
                perc = float(line_data[2].split()[1])
                projection_to_p_intact_M_intact[projection_id] = perc
                continue
            elif line_data[2].startswith("MIDDLE_80%_INTACT"):
                raw_val = line_data[2].split()[1]
                val = True if raw_val == "TRUE" else False
                proj_to_80_p_intact[projection_id] = val
                continue
            elif line_data[2].startswith("MIDDLE_80%_PRESENT"):
                raw_val = line_data[2].split()[1]
                val = True if raw_val == "TRUE" else False
                proj_to_80_p_present[projection_id] = val
                continue
            elif line_data[2].startswith("INTACT_CODONS_PROP"):
                perc = float(line_data[2].split()[1])
                projection_to_i_codon_prop[projection_id] = perc
                continue
            elif line_data[2].startswith("OUT_OF_CHAIN_PROP"):
                perc = float(line_data[2].split()[1])
                proj_to_prop_oub[projection_id] = perc
                continue
            # ordinary mutation line
            exon_num = int(line_data[2])
            codon_num = line_data[3]
            mut_class = line_data[4]
            mut_itself = line_data[5]
            masked = True if line_data[6] == "masked" else False
            mut_ = (exon_num, codon_num, mut_class, mut_itself, masked)
            projection_to_mutations[projection_id].append(mut_)
        f.close()
    output = (projection_to_mutations,
              projection_to_p_intact_M_ignore,
              projection_to_p_intact_M_intact,
              projection_to_i_codon_prop,
              proj_to_prop_oub,
              proj_to_80_p_intact,
              proj_to_80_p_present)
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


def get_L_exon_num(exon_num):
    """20% of exons must be affected to lost the gene."""
    if exon_num == 1:
        return 1
    elif exon_num <= 10:
        return 2
    else:
        twenty_perc = exon_num / 5
        return twenty_perc


def get_projection_classes(all_projections, trans_exon_sizes, p_to_pint_M_ign,
                           p_to_pint_M_int, projection_to_mutations, p_to_i_codon_prop,
                           p_to_p_out_of_bord, p_80_int, p_80_pre,
                           trace=None, paral_=None):
    projection_class = {}
    if paral_ is None:
        # to be able to call if X in paral
        paral = set()
    else:  # if is defined: it's a set
        paral = paral_

    for num, projection in enumerate(all_projections):
        if projection in paral:
            projection_class[projection] = PG
            continue
        p_intact_M_ign = p_to_pint_M_ign.get(projection, -1)
        p_intact_M_int = p_to_pint_M_int.get(projection, -1)
        p_i_codons = p_to_i_codon_prop.get(projection, -1)
        no_loss_in_80_p = p_80_int.get(projection, None)
        m_80_present = p_80_pre.get(projection, None)
        frame_oub = p_to_p_out_of_bord.get(projection, 0.0)
        transcript, _ = split_proj_name(projection)
        tracing_ = trace and transcript == trace

        exon_sizes = trans_exon_sizes.get(transcript)
        if exon_sizes is None:
            print(f"Cannot find transcript {transcript}; probably an error") if tracing_ else None
            projection_class[projection] = N_  # must never happen
            # not APPRIS or so
            continue

        all_mutations = projection_to_mutations.get(projection, [])
        mutations = [m for m in all_mutations if m[4] is False or m[2] == MISS_EXON]  # m[4] = MASKED
        mutations = [m for m in mutations if m[2] != COMPENSATION]
        # transcript_to_projections[transcript].append(projection)

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
            if tracing_:
                print("No mutations, p_M_ign > 0.6: Intact")
            # nothing to talk about
            projection_class[projection] = I
            continue
        elif p_i_codons < REM_T_L:
            if tracing_:
                print(f"In this projection only {p_i_codons} codons remain intact")
                print(f"Threshold is {REM_T_L} codons")
                print(f"-> class L")
            projection_class[projection] = L
            continue

        exon_status = {k: "I" for k in exon_sizes}
        exon_num = len(exon_sizes.keys())
        overall_seq_len = sum(exon_sizes.values())
        exon_40_p = {k: v / overall_seq_len > 0.4 for k, v in exon_sizes.items()}
        exon_40_p_nums = [k for k, v in exon_40_p.items() if v is True]
        
        # select missing and deleted exons
        missing_exons = [m[0] for m in mutations if m[2] == MISS_EXON]
        deleted_exons = [m[0] for m in mutations if m[2] == DEL_EXON]
        for me in missing_exons:
            exon_status[me] = "M"
        for de in deleted_exons:
            exon_status[de] = "D"
        print(f"Exon statuses initial\n{exon_status}") if tracing_ else None
        other_muts = [m for m in mutations if m[2] != MISS_EXON and m[2] != DEL_EXON and m[2] != COMPENSATION]
        other_muts = [m for m in other_muts if m[0] not in deleted_exons and m[0] not in missing_exons]
        missed_seq_len = sum(exon_sizes[k] for k, v in exon_status.items() if v == "M")
        missing_prop = missed_seq_len / overall_seq_len
        if tracing_:
            print(f"% Missing: {missing_prop}")
            print(f"Missing seq len: {missed_seq_len}; overall seq len: {overall_seq_len}")
        affected_thr = get_L_exon_num(exon_num)

        # special cases: all exons are missed or deleted
        if all(v == "M" for v in exon_status.values()):
            print("All exons are Missing: transcript is missing") if tracing_ else None
            print(f"Out of borders prop is {frame_oub}") if tracing_ else None
            if frame_oub > PART_THR:
                print("-> class PM, too big fraction ouf of chain borders") if tracing_ else None
                projection_class[projection] = PM
            else:
                projection_class[projection] = M
            continue
        elif all(v == "D" for v in exon_status.values()):
            print("All exons are Deleted -> class L") if tracing_ else None
            projection_class[projection] = L
            continue

        if p_intact_M_int < 0.2:
            print("P_intact_M_int < 20%: class L") if tracing_ else None
            projection_class[projection] = L
            continue
        
        # main classification process
        if no_loss_in_80_p is True:
            print("GO TO BRANCH 1: No inact mut in m80%") if tracing_ else None
            # first branch -> no inact mutations in the middle 80% of CDS
            # possible: Intact, partially intact, missed
            # grey if too small fraction is left
            # if there are any inact mutations -> it doesn't really matter
            if p_i_codons < REM_T_G:
                print(f"Prop of intact codons is {p_i_codons}") if tracing_ else None
                print(f"Need > {REM_T_G} to be I/PI") if tracing_ else None
                print(f"-> Class Grey") if tracing_ else None
                projection_class[projection] = G
                continue
            if len(missing_exons) == 0:
                print("No missing exons, no inact mut in m80%: class I") if tracing_ else None
                projection_class[projection] = I
                continue
            # m_80_p = m_80_present(exon_sizes, exon_status)
            if m_80_present:
                print("No M exons in m80%: class I") if tracing_ else None
                # means that there are no missing sequence in the middle 80%
                projection_class[projection] = I
                continue
            
            else:
                print("There are M in m80%...") if tracing_ else None
                # partially intact or missed
                if missing_prop < 0.5:
                    print("Missed prop < 50%: class PI") if tracing_ else None
                    # > 50% CDS is presented
                    projection_class[projection] = PI
                    continue
                else:
                    print("Missed prop > 50%: PM or M branch") if tracing_ else None
                    # > 50% of CDS missed, just call it missed
                    if frame_oub > PART_THR:
                        print("-> class PM, too big fraction ouf of chain borders") if tracing_ else None
                        projection_class[projection] = PM
                    else:
                        projection_class[projection] = M
                    continue

        else:
            print("GO TO BRANCH 2: there are inact mut in m80%") if tracing_ else None
            # second major branch: there ARE inact mutations in the middle 80% of CDS
            # possible classes: grey or lost
            for m in other_muts:
                exon_status[m[0]] = "L"
            print(f"Exon status is:\n{exon_status}") if tracing_ else None

            if exon_num == 1:
                print("Single exon branch") if tracing_ else None
                # special case, for a single exon gene we require > 2 mutations
                if exon_status[1] == "D":
                    print("Single exon deleted: class L") if tracing_ else None
                    # the only exon is deleted -> lost
                    projection_class[projection] = L
                elif exon_status[1] == "M":
                    print("Single exon missing: class M") if tracing_ else None
                    # the only exon is missed -> M (reduntant branch highly likely)
                    if frame_oub > PART_THR:
                        print("-> class PM, too big fraction ouf of chain borders") if tracing_ else None
                        projection_class[projection] = PM
                    else:
                        projection_class[projection] = M
                elif p_intact_M_ign < 0.6 and len(other_muts) >= 2:
                    print("%intact < 60% && 2 incat mut: class L") if tracing_ else None
                    projection_class[projection] = L
                else:
                    print("Not enough evidence for L -> G") if tracing_ else None
                    projection_class[projection] = G
                continue
            num_exons_affected = len([k for k, v in exon_status.items() if v == "D" or v == "L"])
            if p_intact_M_int < 0.6:
                print(f"% intact M int < 60 branch") if tracing_ else None
                # well, %intact < 60, maybe a loss!
                if tracing_:
                    print(f"Affected exons: {num_exons_affected}; required: {affected_thr}")
                if num_exons_affected >= affected_thr:
                    # standard way: number of affected exons
                    print(f"Enough affected exons -> L") if tracing_ else None
                    projection_class[projection] = L
                    continue
                # also check whether there is an exon covering > 40% CDS that has TWO inact mutations
                muts_occur = Counter(m[0] for m in other_muts)
                muts_in_40_exons = [muts_occur[x] for x in exon_40_p_nums]
                if any(x >= 2 for x in muts_in_40_exons):
                    print(f"There are exons > 40% with 2+ mutations -> L") if tracing_ else None
                    projection_class[projection] = L
                    continue
                if any(exon_status[x] == "D" for x in exon_40_p_nums):
                    print(f"Some of exons > 40% are Deleted -> Lost") if tracing_ else None
                    projection_class[projection] = L
                    continue
                print(f"Not enough evidence for Lost -> Grey") if tracing_ else None
                projection_class[projection] = G
            else:
                print(f"% intact M int > 60 branch") if tracing_ else None
                print(f"not enough evidence for L -> Grey") if tracing_ else None
                # if %intact > 60: cannot be intact
                projection_class[projection] = G
                continue
    return projection_class


def get_exon_sizes(ref_bed):
    """Get exon: size dict for reference bed."""
    trans_exon_sizes = {}
    f = open(ref_bed, "r")
    for line in f:
        cds_line = make_cds_track(line)
        line_data = cds_line.split("\t")
        strand = line_data[5]
        exon_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
        if strand == "-":
            exon_sizes = exon_sizes[::-1]
        exon_to_size = {n: v for n, v in enumerate(exon_sizes, 1)}
        trans_id = line_data[3].replace("_CDS", "")
        trans_exon_sizes[trans_id] = exon_to_size
    f.close()
    return trans_exon_sizes


def read_isoforms(isoforms_file, all_transcripts):
    """Read isoforms file."""
    gene_to_trans = defaultdict(list)
    f = open(isoforms_file, "r")
    f.__next__()
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        trans = line_data[1]
        if trans not in all_transcripts:
            continue
        gene_to_trans[gene].append(trans)
    f.close()
    return gene_to_trans


def get_paralogs_data(paral_file):
    """Extract paralogous projections."""
    if paral_file is None:
        return set()
    paral_proj = []
    with open(paral_file, "r") as f:
        paral_proj = set(x.rstrip() for x in f.readlines())
    return paral_proj


def color_bed_file(bed_in, bed_out, proj_to_class):
    """Assing colors to bed tracks."""
    in_ = open(bed_in, "r")
    out_ = open(bed_out, "w")
    for line in in_:
        line_data = line.rstrip().split("\t")
        projection_id = line_data[3]
        projection_class = proj_to_class.get(projection_id, N_)
        # N_ must never occur, TODO: raise an error
        color = CLASS_TO_COL[projection_class]
        line_data[8] = color
        line_upd = "\t".join(line_data)
        out_.write(line_upd)
        out_.write("\n")
    in_.close()
    out_.close()


def gene_losses_summary(loss_data_arg, ref_bed, pre_final_bed_arg,
                        bed_out, summary_arg, trace_arg=None,
                        iforms=None, paral=None):
    """Entry point."""
    t0 = dt.now()
    # classify projections
    paralogs_set = get_paralogs_data(paral)
    trans_exon_sizes = get_exon_sizes(ref_bed)
    loss_data_all = read_loss_data(loss_data_arg)
    projection_to_mutations = loss_data_all[0]
    p_to_pintact_M_ign = loss_data_all[1]
    p_to_pintact_M_int = loss_data_all[2]
    p_to_i_codon_prop = loss_data_all[3]
    p_to_oub_prop = loss_data_all[4]
    p_80_int = loss_data_all[5]
    p_80_pre = loss_data_all[6]

    all_projections = set(p_to_pintact_M_ign.keys()).union(set(projection_to_mutations.keys()))
    projection_class = get_projection_classes(all_projections,
                                              trans_exon_sizes,
                                              p_to_pintact_M_ign,
                                              p_to_pintact_M_int,
                                              projection_to_mutations,
                                              p_to_i_codon_prop,
                                              p_to_oub_prop,
                                              p_80_int,
                                              p_80_pre,
                                              trace=trace_arg,
                                              paral_=paralogs_set)

    # color bed file
    color_bed_file(pre_final_bed_arg, bed_out, projection_class)    

    # get transcript to projections dict
    transcript_to_projections = defaultdict(list)
    for proj in all_projections:
        trans, _ = split_proj_name(proj)
        transcript_to_projections[trans].append(proj)
    transcript_class = {}

    # classify transcripts
    for trans, projections in transcript_to_projections.items():
        p_classes = set(projection_class.get(p) for p in projections)
        status = max(p_classes)  # just use Enum I > G > L > M > N
        transcript_class[trans] = status
    all_transcripts = set(transcript_class.keys())

    # classify genes
    if iforms:
        gene_to_trans = read_isoforms(iforms, all_transcripts)
    else:
        gene_to_trans = {}

    gene_class = {}
    for gene, transcripts in gene_to_trans.items():
        trans_statuses = [transcript_class.get(t, -1) for t in transcripts]
        status = max(trans_statuses)
        gene_class[gene] = status
    
    # save summary
    f = open(summary_arg, "w")
    for k, v in projection_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"PROJECTION\t{k}\t{v_ch}\n")
    for k, v in transcript_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"TRANSCRIPT\t{k}\t{v_ch}\n")
    for k, v in gene_class.items():
        v_ch = NUM_TO_CLASS.get(v, "N")
        f.write(f"GENE\t{k}\t{v_ch}\n")
    f.close()
    print(f"Elapsed: {dt.now() - t0}")


def main():
    """Entry point for CLI."""
    args = parse_args()
    gene_losses_summary(args.loss_data,
                        args.ref_bed,
                        args.pre_final_bed,
                        args.bed_out,
                        args.summary,
                        trace_arg=args.trace,
                        iforms=args.isoforms,
                        paral=args.paral_projections)


if __name__ == "__main__":
    main()
