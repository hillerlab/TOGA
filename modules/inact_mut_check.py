#!/usr/bin/env python3
"""Scan reading frame for inactivating mutations."""
import sys
import argparse
from re import finditer, IGNORECASE
from collections import defaultdict
from collections import Counter
from collections import namedtuple
try:
    from modules.parse_cesar_output import parse_cesar_out
    from modules.parse_cesar_output import classify_exon
except ImportError:
    from parse_cesar_output import parse_cesar_out
    from parse_cesar_output import classify_exon

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


STOPS = {"TAG", "TAA", "TGA"}
D_M = {"D", "M"}
left_splice_corr = ("ag", )  # acceptor
right_splice_corr = ("gt", "gc", )  # donor
pattern = r"-{1,}"
LEFT_SSID = 0
RIGHT_SSID = 1
ACCEPTOR = 0
DONOR = 1

BIG_INDEL_SIZE = 50
SAFE_EXON_DEL_SIZE = 40  # actually 39
FIRST_LAST_DEL_SIZE = 20
BIG_EXON_THR = BIG_INDEL_SIZE * 5
# mutation "struc"
# gene, chain, exon -> exon identifier
# position -> num of codon, 0 if entire exon missed/deleted
# class -> SSM, FS_DEL, FS_INS, STOP, DELETED, MISSED, BIG_DEL, BIG_INS
# mut:
#    if stop -> stop codon
#    if INS/DEL -> ins/del size
#    if SSM -> show what exactly
# masked -> if in first/last 10% -> ignore
Mutation = namedtuple("Mutation", "gene chain exon position mclass mut masked mut_id")
# mut classes
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

# if False: frame-preserving indels might be inactivating
# if true: no more big indels
# frame-preserving exon deletions are not inact mutations
IGNORE_FP_INDELS = True


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(str(msg) + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("cesar_output", help="File containing raw CESAR output")
    app.add_argument("--verbose", "-v", dest="verbose", action="store_true")
    app.add_argument("--gene", default="None", help="Gene name")
    app.add_argument("--u12", default=None)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_cesar_out(cesar_line):
    """Return ref and query sequence."""
    cesar_content = cesar_line.split("\n")
    # del cesar_content[0]
    fractions = parts(cesar_content, 4)
    cesar_fractions = []
    for fraction in fractions:
        if len(fraction) == 1:
            continue
        ref_seq = fraction[1]
        query_name = fraction[2][1:]
        query_seq = fraction[3]
        if len(ref_seq) != len(query_seq):
            die("Error! Ref and query sequences must have the same length!")
        elif len(ref_seq) == 0:
            die("Error! The input is empty!")
        fraction = (query_name, ref_seq, query_seq)
        cesar_fractions.append(fraction)
    return cesar_fractions


def parse_u12_opt(gene, u12_data):
    """Parse U12 introns data."""
    if u12_data is None:
        return set()
    if gene == "None":
        return set()

    ans = set()
    f = open(u12_data, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        if trans != gene:
            continue
        exon_num = int(line_data[1])
        side = 0 if line_data[2] == "A" else 1
        u12_site = (exon_num, side)
        ans.add(u12_site)
    f.close()
    return ans


def mask_mut(mut):
    """Mutation is immutable, need a func to return a masked version."""
    # namedtuple("Mutation", "gene chain exon position mclass mut masked mut_id")
    gene = mut.gene
    chain = mut.chain
    exon = mut.exon
    position = mut.position
    mclass = mut.mclass
    mut_ = mut.mut
    masked = True
    mut_id = mut.mut_id
    upd_mut = Mutation(gene, chain, exon, position, mclass, mut_, masked, mut_id)
    return upd_mut


def analyse_splice_sites(ref, query, gene, chain, u12_data=None, v=None):
    """Check correctness of the splice sites."""
    u12_data = set() if not u12_data else u12_data
    mut_counter = 1
    eprint(f"U12 introns: {u12_data}") if v else None
    sps_report = []
    cds_indexes = [i for i, c in enumerate(ref) if c != " " and c != ">"]
    exon_num = 0
    exon_num_indexes = defaultdict(list)

    for i in range(1, len(cds_indexes)):
        prev = cds_indexes[i - 1]
        current = cds_indexes[i]
        delta = current - prev
        exon_num_indexes[exon_num].append(prev)
        if delta > 1:
            exon_num += 1

    switch = False
    exon_num = 0
    indexes_exon_num = {}
    cds_index = 0
    for i in range(len(ref)):
        i_ref = ref[i]
        if i_ref == " ":
            switch = True
        elif i_ref != " " and switch:
            switch = False
            exon_num += 1
            indexes_exon_num[cds_index] = exon_num
            cds_index += 1
        elif i_ref != " " and not switch:
            indexes_exon_num[cds_index] = exon_num
            cds_index += 1

    exon_num_start_end = {k: (min(v), max(v)) for k, v in exon_num_indexes.items()}
    exons_num = len(exon_num_start_end)

    for exon_num_, start_end in exon_num_start_end.items():
        exon_num = exon_num_ + 1
        # u12_sps = u12_data.get(true_exon_num)
        start, end = start_end
        left_splice_site = query[start - 2: start]
        left_splice_site_N = "n" in left_splice_site or "N" in left_splice_site
        right_splice_site = query[end + 1: end + 3]
        right_splice_size_N = "n" in right_splice_site or "N" in right_splice_site

        left_site_wrong = left_splice_site not in left_splice_corr
        right_site_wrong = right_splice_site not in right_splice_corr
        eprint(f"Exon {exon_num}; L_SPS: {left_splice_site}; R_SPS: {right_splice_site}") if v else None

        # if u12_sps:
        #     right_site_wrong = False if RIGHT_SSID in u12_sps else right_site_wrong
        #     left_site_wrong = False if LEFT_SSID in u12_sps else left_site_wrong

        # something wrong with NN and --
        if exon_num != 1 and left_site_wrong:
            mask = True if (exon_num, 0) in u12_data else False
            mask = True if left_splice_site_N else mask
            mut_ = f"{left_splice_corr}->{left_splice_site}"
            mut_id = f"SSM_{mut_counter}"
            mut_counter += 1
            mut = Mutation(gene=gene, chain=chain, exon=exon_num, position=0,
                           mclass=SSM, mut=mut_, masked=mask, mut_id=mut_id)
            sps_report.append(mut)
            
        if exon_num != exons_num and right_site_wrong:
            mask = True if (exon_num, 1) in u12_data else False
            mask = True if right_splice_size_N else mask
            mut_ = f"{right_splice_corr}->{right_splice_site}"
            mut_id = f"SSM_{mut_counter}"
            mut_counter += 1
            mut = Mutation(gene=gene, chain=chain, exon=exon_num, position=1,
                           mclass=SSM, mut=mut_, masked=mask, mut_id=mut_id)
            sps_report.append(mut)
    if len(sps_report) <= 1:
        return sps_report
    # if more than 1 -> exclude intron deletions
    sps_report = sorted(sps_report, key=lambda x: (x.exon, x.position))
    for i in range(1, len(sps_report)):
        j = i - 1
        prev = sps_report[j]
        curr = sps_report[i]
        curr_exon = curr.exon
        prev_exon = prev.exon
        if curr_exon != prev_exon + 1:
            continue
        curr_pos = curr.position
        prev_pos = prev.position
        if not (curr_pos == 0 and prev_pos == 1):
            continue
        prev_to_what = prev.mut.split("->")[1]
        curr_to_what = curr.mut.split("->")[1]
        if prev_to_what == curr_to_what == "--":
            # intron deletion
            sps_report[j] = mask_mut(prev)
            sps_report[i] = mask_mut(curr)
        else:
            continue
    return sps_report


def split_to_codons(ref, que):
    """Create codon_num->(ref_seq, que_seq) dict."""
    ref_not_space = [i for i, c in enumerate(ref) if c != " "]
    ref_merge, que_merge = "", ""
    for i in ref_not_space:
        ref_merge += ref[i]
        que_merge += que[i]
    ref_codon_to_seq = defaultdict(str)
    que_codon_to_seq = defaultdict(str)
    curr_codon = 1
    for i in range(len(ref_merge)):
        ref_codon_to_seq[curr_codon] += ref_merge[i].upper()
        que_codon_to_seq[curr_codon] += que_merge[i].upper()
        ref_codon_len = len([c for c in ref_codon_to_seq[curr_codon] if c.isalpha()])
        if ref_codon_len == 3:
            curr_codon += 1
    return ref_codon_to_seq, que_codon_to_seq


def get_codon_to_exon_nums(ref):
    """Get codon_num -> exon_num correspondence."""
    ref_parts = [x for x in ref.split(" ") if x]
    prev_rem = 0
    curr_codon = 1
    codon_to_exon = {}

    for num, exon_seq in enumerate(ref_parts, 1):
        lett_only = [c for c in exon_seq if c.isalpha()]
        ex_len = len(lett_only) - prev_rem
        full_codons = ex_len // 3
        rem = ex_len % 3
        prev_rem = 0 if rem == 0 else 3 - rem
        codons_num = full_codons if rem == 0 else full_codons + 1
        for _ in range(codons_num):
            codon_to_exon[curr_codon] = num
            curr_codon += 1
    return codon_to_exon


def corr_exon_num_or_no_fs(codon, exon_num):
    """Need to decide which codon num to assign."""
    # get left size
    cut_at = codon["split_"]
    left_side_ref = codon["ref_codon"][:cut_at]
    left_side_que = codon["que_codon"][:cut_at]
    ls_ref_gaps = left_side_ref.count("-")
    ls_que_gaps = left_side_que.count("-")
    fs_left = abs(ls_ref_gaps - ls_que_gaps) % 3 != 0
    if fs_left:
        return exon_num - 1
    else:
        return exon_num


def scan_rf(codon_table, gene, chain, exon_stat=None, v=False,
            big_indel_thrs=None, sec_codons=None, no_fpi=False):
    """Scan codon table for inact mutations."""
    sec_codons_set = sec_codons if sec_codons else set()
    codons_num = len(codon_table)
    in_mut_report = []
    perc_10 = codons_num // 10
    left_t = perc_10
    right_t = codons_num - perc_10
    mut_number_stop = 1
    mut_number_fs = 1
    mut_number_big_indel = 1
    q_dels_in_a_row = 0

    for num, codon in enumerate(codon_table, 1):
        mask = True if num <= left_t or num >= right_t else False
        last_codon = True if num == codons_num else False
        first_codon = True if num == 1 else False

        ex_num = codon["t_exon_num"] + 1
        ref_codon = codon["ref_codon"]
        que_codon = codon["que_codon"]
        codon_is_split = codon["split_"] > 0
        # it depends on the exon length:
        big_indel_thr = big_indel_thrs[ex_num] if big_indel_thrs else BIG_INDEL_SIZE

        if que_codon == "---" and codon_is_split is False:
            q_dels_in_a_row += 1
        else:
            big_ins_cond = q_dels_in_a_row * 3 > big_indel_thr
            if big_ins_cond and no_fpi is True:
                # no fpi is False -> fpi is True -> frame pres indels ignored
                # no fpi is True -> fpi is False -> frame pres indels not ignored
                position = num - q_dels_in_a_row  # -1 +1
                # if FS and big ins at the same time???
                mut_ = f"-{q_dels_in_a_row * 3}"
                mclass = BIG_DEL
                mut_id = f"BI_{mut_number_big_indel}"
                mut_number_big_indel += 1
                mut = Mutation(gene=gene, chain=chain, exon=ex_num, position=position,
                               mclass=mclass, mut=mut_, masked=mask, mut_id=mut_id)
                in_mut_report.append(mut)
            q_dels_in_a_row = 0

        if first_codon and que_codon != "ATG":
            mut_id = "START_1"
            mut = Mutation(gene=gene, chain=chain, exon=ex_num, position=1,
                           mclass=START_MISSING, mut=que_codon, masked=mask,
                           mut_id=mut_id)
            in_mut_report.append(mut)

        # check for FS
        ref_gaps = ref_codon.count("-")
        que_gaps = que_codon.count("-")
        delta = ref_gaps - que_gaps
        fs = abs(delta) % 3 != 0

        if fs: 
            # we have a frameshifing indel!
            mut_ = f"+{delta}" if delta > 0 else f"{delta}"
            mclass = FS_INS if delta > 0 else FS_DEL
            mut_id = f"FS_{mut_number_fs}"
            mut_number_fs += 1
            if not codon_is_split:
                fs_ex_num = ex_num
            else:
                # we need to decide which exon is it
                fs_ex_num = corr_exon_num_or_no_fs(codon, ex_num)
            mut = Mutation(gene=gene, chain=chain, exon=fs_ex_num, position=num,
                           mclass=mclass, mut=mut_, masked=mask, mut_id=mut_id)
            in_mut_report.append(mut)
            if v:  # verbose
                eprint("Detected FS")
                eprint(codon)
            mask = True  # to must all next inact mutations in this codon
        if delta > big_indel_thr and no_fpi is True:
            # big insertion!
            # if FS and big ins at the same time???
            mut_ = f"+{delta}"
            mclass = BIG_INS
            mut_id = f"BI_{mut_number_big_indel}"
            mut_number_big_indel += 1
            if not codon_is_split:
                bi_ex_num = ex_num
            else:
                # we need to decide which exon is it
                bi_ex_num = corr_exon_num_or_no_fs(codon, ex_num)
            mut = Mutation(gene=gene, chain=chain, exon=bi_ex_num, position=num,
                           mclass=mclass, mut=mut_, masked=mask, mut_id=mut_id)
            in_mut_report.append(mut)
            mask = True

        # check for inframe stop codons
        que_codon_no_gap = que_codon.replace("-", "")
        if not que_codon_no_gap:
            continue
        triplets = parts(que_codon_no_gap, n=3)
        stop_triplets = [x for x in triplets if x in STOPS]
        if len(stop_triplets) > 0 and not last_codon:
            # we have stops!
            mut_ = f"{ref_codon.replace('-', '')}->{stop_triplets[0]}"
            mclass = STOP
            mut_id = f"STOP_{mut_number_stop}"
            if not codon_is_split:
                # not split, no correction
                st_ex_num = ex_num
            else:
                prev_exon = ex_num - 1
                prev_exon_stat = exon_stat[prev_exon] if exon_stat else "I"
                this_exon_stat = exon_stat[ex_num] if exon_stat else "I"
                if this_exon_stat in D_M or prev_exon_stat in D_M:
                    # split codon, one of exons is D/M: skip this
                    continue
                st_ex_num = ex_num if codon["split_"] == 1 else ex_num - 1
            mut_number_stop += 1
            is_sec_pos = num - 1 in sec_codons_set  # 0-based in that set
            if is_sec_pos and stop_triplets[0] == "TGA":
                # mask if U-coding codon
                stop_mask = True
            else:  # apply usual rules
                stop_mask = mask
            mut = Mutation(gene=gene, chain=chain, exon=st_ex_num, position=num,
                           mclass=mclass, mut=mut_, masked=stop_mask, mut_id=mut_id)
            in_mut_report.append(mut)

            if v:
                eprint("Detected STOP")
                eprint(codon)
    # I can inder number of codons in each exon directly from codon table
    return in_mut_report


def detect_compensations(inact_mut, codon_table):
    """Detect compensation events."""
    m_counter = 1
    answer = []
    fs = [mut for mut in inact_mut if mut.mclass in {FS_DEL, FS_INS}]
    fs_num = len(fs)
    if fs_num <= 1:
        # need at least 2 frameshifs
        return answer
    potent_compensations = []  # based only on FS values
    for i_num in range(fs_num - 1):
        # skip the last mutation, it cannot be compensated
        init_mut = fs[i_num]
        init_mut_val = int(init_mut.mut)
        init_mut_id = init_mut.mut_id
        fs_values = [init_mut_val, ]
        fs_ids = [init_mut_id, ]
        # iter over next mutations
        for j_num in range(i_num + 1, fs_num):
            j_mut = fs[j_num]
            j_mut_val = int(j_mut.mut)
            j_mut_id = j_mut.mut_id
            fs_values.append(j_mut_val)
            fs_ids.append(j_mut_id)
            # check whether it's potential compensation
            if sum(fs_values) % 3 == 0:
                potent_compensations.append(fs_ids)
                break
    if len(potent_compensations) == 0:
        return []  # no potential compensations, skip this

    # verify potential compensations
    what_is_compensated = set()  # to avoid twice> compensated FS
    for comp in potent_compensations:
        if len(what_is_compensated.intersection(comp)) > 0:
            # it means that some of these FS are already compensated
            continue
        comp_muts = sorted([m for m in fs if m.mut_id in comp], key=lambda x: x.position)
        start_pos = comp_muts[0].position
        end_pos = comp_muts[-1].position
        # positions are 1-based, need to correct
        codons_seq = codon_table[start_pos - 1: end_pos]
        Q_seq = "".join([c["que_codon"] for c in codons_seq]).replace("-", "")
        upd_codons = parts(Q_seq, n=3)
        if len(STOPS.intersection(upd_codons)) > 0:
            # there are stops in the compensated sequence
            continue
        # add compensation track
        ethalon_mut = comp_muts[0]
        gene = ethalon_mut.gene
        chain = ethalon_mut.chain
        exon = ethalon_mut.exon  # not applicable mostly, but
        position = ethalon_mut.position
        mclass = COMPENSATION
        mut_id = f"C_{m_counter}"
        fs_ids = [x.split("_")[1] for x in comp]
        mut = "FS_{}".format(",".join(fs_ids))
        comp_mut = Mutation(gene=gene, chain=chain, exon=exon, position=position,
                            mclass=mclass, mut=mut, masked=False, mut_id=mut_id)
        answer.append(comp_mut)
        for c in comp:
            what_is_compensated.add(c)
        m_counter += 1
    return answer


def mask_compensated_fs(mut_list):
    """Mask compensated mutations."""
    comp_muts = [m for m in mut_list if m.mclass == COMPENSATION]
    if len(comp_muts) == 0:
        # no compensations, no worries
        return mut_list
    comp_fs_ids = []
    for cmp in comp_muts:
        fs_numbs = cmp.mut.split("_")[1].split(",")
        fs_ids = [f"FS_{x}" for x in fs_numbs]
        comp_fs_ids.extend(fs_ids)
    comp_fs_ids = set(comp_fs_ids)
    not_touch = [m for m in mut_list if m.mut_id not in comp_fs_ids]
    to_mask = [m for m in mut_list if m.mut_id in comp_fs_ids]
    masked = [mask_mut(m) for m in to_mask]
    filtered = not_touch + masked
    return filtered


def compute_percent_id(seq_1, seq_2):
    """Return % identity for two sequences."""
    assert len(seq_1) == len(seq_2)  # otherwise it is a bug
    matches = sum([1 for i in range(len(seq_1))
                  if seq_1[i] == seq_2[i]
                  and seq_1[i] != "N"
                  and seq_2[i] != "N"
                  and seq_1[i] != "-"
                  and seq_2[i] != "-"])
    length = len(seq_1.replace("N", "").replace("-", ""))  # we ignore N's
    pid = matches * 100 / length if length != 0 else 0
    return pid



def classify_exons(gene, que, codon_table, exon_class, exon_gap,
                   exon_pid, exon_blosum, missing_exons, ex_inc):
    """Classify exons as intact, deleted and missing."""
    del_num, miss_num = 1, 1
    exon_nums = list(range(codon_table[-1]["t_exon_num"] + 1))
    exons_report = []
    exon_stat = ["X", ]
    for exon_num in exon_nums:
        ex_num_ = exon_num + 1
        if exon_num in missing_exons:
            # a priori missing
            # missing exons list is 0-based
            # exon_num - 0 based, exon_num_ (with underscore) is 1-based
            mut = Mutation(gene=gene, chain=que, exon=ex_num_, position=0, mut="-",
                           mclass=MISS_EXON, masked=False, mut_id=f"MIS_{miss_num}")
            miss_num += 1
            exons_report.append(mut)
            exon_stat.append("M")
            continue
        ex_class = exon_class.get(exon_num, None)
        ex_gap = exon_gap.get(exon_num, None)
        ex_pid = exon_pid.get(exon_num, 0)
        ex_blosum = exon_blosum.get(exon_num, 0)
        exon_excl = ex_inc.get(exon_num, None)
        del_, q = classify_exon(ex_class, exon_excl, ex_pid, ex_blosum)

        if ex_class == "M" or ex_gap:
            mut = Mutation(gene=gene, chain=que, exon=ex_num_, position=0, mut="-",
                           mclass=MISS_EXON, masked=False, mut_id=f"MIS_{miss_num}")
            miss_num += 1
            exons_report.append(mut)
            exon_stat.append("M")
            continue
        elif del_ is False:
            mut = Mutation(gene=gene, chain=que, exon=ex_num_, position=0, mut="-",
                           mclass=DEL_EXON, masked=False, mut_id=f"DEL_{del_num}")
            del_num += 1
            exons_report.append(mut)
            exon_stat.append("D")
            continue
        else:
            exon_stat.append("I")
            pass
    return exons_report, exon_stat


def muts_to_text(mutations, perc_intact_1, perc_intact_2, i_prop, oub, m_80_i, m_80_p, gene):
    """Convert mutations array to text."""
    ordered = sorted(mutations, key=lambda x: (x.chain, x.exon, x.position))
    strings = []
    for elem in ordered:
        elem_values = list(elem._asdict().values())
        # last field masked values False and True migth be confusing
        if elem_values[6] is False:
            elem_values[6] = "unmasked"
        else:
            elem_values[6] = "masked"
        elem_values = [str(x) for x in elem_values]
        # hash-tag to distinguish CESAR output from GeneLossScanner
        elem_string = "# " + "\t".join(elem_values)
        strings.append(elem_string)
    for q, val in perc_intact_1.items():
        p_intact_line = f"# {gene}\t{q}\tINTACT_PERC_IGNORE_M {val}"
        strings.append(p_intact_line)
    for q, val in perc_intact_2.items():
        p_intact_line = f"# {gene}\t{q}\tINTACT_PERC_INTACT_M {val}"
        strings.append(p_intact_line)
    for q, val in i_prop.items():
        i_prop_line = f"# {gene}\t{q}\tINTACT_CODONS_PROP {val}"
        strings.append(i_prop_line)
    for q, val in oub.items():
        out_of_ch_line = f"# {gene}\t{q}\tOUT_OF_CHAIN_PROP {val}"
        strings.append(out_of_ch_line)
    for q, val in m_80_i.items():
        val_str = "TRUE" if val is True else "FALSE"
        m_80_line = f"# {gene}\t{q}\tMIDDLE_80%_INTACT {val_str}"
        strings.append(m_80_line)
    for q, val in m_80_p.items():
        val_str = "TRUE" if val is True else "FALSE"
        m_80_line = f"# {gene}\t{q}\tMIDDLE_80%_PRESENT {val_str}"
        strings.append(m_80_line)
    return "\n".join(strings) + "\n"


def get_exon_num_corr(codons_data):
    """Get correspondence between exon numbers in Q and T."""
    ans = defaultdict(set)
    for codon in codons_data:
        ans[codon["q_exon_num"]].add(codon["t_exon_num"])
    if not ans[0]:
        ans[0] = {0}
    return ans


def compute_intact_perc(codon_table, mutations, q_name, v=False):
    """Compute intact %ID."""
    query_muts = [m for m in mutations if m.chain == q_name]
    gene_len = len(codon_table)
    # codon_valid = [True for _ in codon_table]
    # initiate codon_status, mark deleted codons with D, the rest with I
    codon_status = ["I" if c["que_codon"] != "---" else "D" for c in codon_table]
    # del_miss_mut_exons = {m.exon - 1 for m in query_muts if m.mclass in DEL_MISS}
    del_exons = {m.exon - 1 for m in query_muts if m.mclass == DEL_EXON and m.masked is False}
    safe_del_exons = {m.exon - 1 for m in query_muts if m.mclass == DEL_EXON and m.masked is True}
    miss_exons = {m.exon - 1 for m in query_muts if m.mclass == MISS_EXON}
    # inval_nums = [n for n, c in enumerate(codon_table) if c["t_exon_num"] in del_miss_mut_exons]
    # for inval_num in inval_nums:
    #     codon_valid[inval_num] = False
    del_codon_nums = [n for n, c in enumerate(codon_table) if c["t_exon_num"] in del_exons]
    safe_del_codon_nums = [n for n, c in enumerate(codon_table) if c["t_exon_num"] in safe_del_exons]
    miss_codon_nums = [n for n, c in enumerate(codon_table) if c["t_exon_num"] in miss_exons]

    for del_codon in del_codon_nums:
        codon_status[del_codon] = "L"
    for del_codon in safe_del_codon_nums:
        codon_status[del_codon] = "D"
    for miss_codon in miss_codon_nums:
        codon_status[miss_codon] = "M"
    
    compensations = [m.mut for m in query_muts if m.mclass == COMPENSATION]
    comp_nums = ",".join([c.split("_")[1] for c in compensations]).split(",")
    comp_fs = {f"FS_{c}" for c in comp_nums}
    
    for m in query_muts:
        if m.mut_id in comp_fs:
            continue
        elif m.mclass in DEL_MISS:
            continue
        elif m.mclass == COMPENSATION:
            continue
        elif m.mclass == SSM:
            if m.masked is True:
                # U12
                continue
            to_what = m.mut.split("->")[1]
            if to_what == "??" or to_what.upper() == "NN":
                # we don't know actually
                continue
            ssm_exon = m.exon - 1
            codon_pos_at_exon = [n for n, c in enumerate(codon_table) if c["t_exon_num"] == ssm_exon]
            if len(codon_pos_at_exon) == 0 and ssm_exon == 0:
                affected_num = 0
            elif m.position == 0:
                affected_num = codon_pos_at_exon[0]
            else:
                affected_num = codon_pos_at_exon[-1]
            if codon_status[affected_num] != "M":
                codon_status[affected_num] = "L"
            continue
        elif m.mclass == "START_MISSED":
            continue

        affected_num = m.position - 1
        if codon_status[affected_num] != "M":
            codon_status[affected_num] = "L"

    if all(x is "I" for x in codon_status):
        # nearly impossible
        return 1.0, 1.0, 1.0, True, True

    if v:
        eprint(f"Num of codons: {len(codon_status)}")
        eprint(f"Num of I: {codon_status.count('I')}")
        eprint(f"Num of D: {codon_status.count('D')}")
        eprint(f"Num of M: {codon_status.count('M')}")
        eprint(f"Num of L: {codon_status.count('L')}")
    codon_status_string = "".join(codon_status)
    repl_d = codon_status_string.replace("D", "")
    ignore_m = repl_d.replace("M", "")
    m_intact = repl_d.replace("M", "I")

    ignore_m_spans = ignore_m.split("L")
    m_intact_spans = m_intact.split("L")

    ignore_m_max_span = max(len(x) for x in ignore_m_spans)
    m_intact_max_span = max(len(x) for x in m_intact_spans)

    p_intact_ignore_m = ignore_m_max_span / gene_len
    p_intact_intact_m = m_intact_max_span / gene_len

    ten_perc = gene_len // 10
    middle = codon_status_string[ten_perc : - ten_perc]
    middle_80_intact = False if "L" in middle else True
    middle_80_present = False if "M" in middle else True
    not_m = len(codon_status) - codon_status.count("M")
    if not_m > 0:
        num_of_I_codons = codon_status.count("I") / not_m
    else:  # meaning all codons are missing
        num_of_I_codons = 1.0

    return p_intact_ignore_m, p_intact_intact_m, num_of_I_codons, middle_80_intact, middle_80_present


def filter_mutations(mut_list):
    """Remove mutations of del or missed exons."""
    del_exon_muts = [m for m in mut_list if m.mclass == DEL_EXON]
    mis_exon_muts = [m for m in mut_list if m.mclass == MISS_EXON]
    del_exons = [m.exon for m in del_exon_muts]
    missing_exons = [m.exon for m in mis_exon_muts]
    del_and_missed = set(missing_exons + del_exons)
    mut_to_keep = [m for m in mut_list if m.exon not in del_and_missed]
    filtered = mut_to_keep + del_exon_muts + mis_exon_muts
    return filtered


def get_D_runs(ex_stat):
    """Get D runs."""
    d_inds = [n for n, v in enumerate(ex_stat) if v == "D" or v == "mD"]
    if len(d_inds) <= 1:
        # nothing we can do
        return []
    result = []
    curr_list = [d_inds[0], ]
    for elem in d_inds[1:]:
        prev_elem = curr_list[-1]
        if prev_elem + 1 == elem:
            # append to the current list
            curr_list.append(elem)
        else:
            # start new list
            curr_copy = curr_list.copy()
            result.append(curr_copy)
            curr_list = [elem, ]
    result.append(curr_list)
    final_res = [x for x in result if len(x) > 1]
    return result


def find_safe_ex_dels(mut_list, ex_stat_, ex_lens, no_fpi=False):
    """Select safe exon deletions."""
    if ex_lens is None:
        return mut_list, ex_stat_
    upd_mut_list = []
    mdel_num = 1
    compensated_ex = set()
    # get compensated exdels
    D_runs = get_D_runs(ex_stat_)
    for D_run in D_runs:
        if len(D_run) < 2:
            # single ex del cannot compensate itself
            continue
        D_lens = [ex_lens[d] for d in D_run]
        D_run_len = len(D_run)
        for i_num in range(D_run_len - 1):
            d_elem = D_run[i_num]
            if d_elem in compensated_ex:
                continue
            d_len = D_lens[i_num]
            d_elems = [d_elem, ]
            d_lens = [d_len, ]
            # maybe there is a bit more elegant solution?
            for j_num in range(i_num + 1, D_run_len):
                j_elem = D_run[j_num]
                j_len = D_lens[j_num]
                d_elems.append(j_elem)
                d_lens.append(j_len)
                if sum(d_lens) % 3 != 0:
                    continue
                for x in d_elems:
                    compensated_ex.add(x)
                break

    # filter D/M mutations
    for m in mut_list:
        ex_len = ex_lens[m.exon]
        is_first = m.exon == 1
        is_last = m.exon == len(ex_lens)
        last_or_first = is_first or is_last
        frame_pres = ex_len % 3 == 0 or m.exon in compensated_ex
        # if we ignore FS indels at all: then any exon size fits, always true
        # no fpi is True -> frame-pres indels not ignored -> ignore deletions of short exons only
        # no fpi is False -> ignore all frame-preserving exon deletions
        fp_ex_len_cond = ex_len < SAFE_EXON_DEL_SIZE if no_fpi is True else True

        if last_or_first and ex_len < FIRST_LAST_DEL_SIZE:
            # then we say it's masked
            new_m = Mutation(gene=m.gene, chain=m.chain, exon=m.exon, position=m.position,
                             mclass=MISS_EXON, mut_id=f"MDEL_{mdel_num}", mut=m.mut, masked=m.masked)
            mdel_num += 1
            upd_mut_list.append(new_m)
            ex_stat_[m.exon] = "M"
            continue
        elif frame_pres and fp_ex_len_cond:
            masked_m = mask_mut(m)
            upd_mut_list.append(masked_m)
            ex_stat_[m.exon] = "mD"
        else:
            upd_mut_list.append(m)
    return upd_mut_list, ex_stat_


def infer_big_indel_thresholds(ex_lens):
    """For each exon define big indel threshold."""
    ex_T = {}
    if ex_lens is None:
        return {}
    for ex_num, ex_len in ex_lens.items():
        if ex_len <= BIG_EXON_THR:
            ex_T[ex_num] = BIG_INDEL_SIZE
            continue
        # then take 20%
        thr = int(ex_len / 5)
        ex_T[ex_num] = thr
    return ex_T


def get_exon_pairs(exon_stat):
    """Get pairs of I exons to check for split stop codons."""
    pairs = []
    pair_init = None
    for num, stat in enumerate(exon_stat):
        if num == 0:  # X - placeholder
            continue
        prev_num = num - 1
        prev_stat = exon_stat[prev_num]
        
        if stat == "I" and pair_init is None:
            continue
        elif stat == "I" and pair_init:
            # pair is initiated, let's close it
            pair = (pair_init, num)
            pairs.append(pair)
            pair_init = None
        
        if stat == "M":
            # something like I-D-D-M-D-I -> then we skip it
            pair_init = None
            continue
        
        if stat == "D" or stat == "mD":
            # the most interesting case
            # deleted or masked deleted
            if pair_init is None and prev_stat == "I":
                # initiate new pair
                # prev elem was I
                pair_init = prev_num
                continue
            else:
                continue
    return pairs


def detect_split_stops(codon_table, gene, q_name, exon_stat):
    """Considering all exon deletions find all split stop codons."""
    i_exon_pairs = get_exon_pairs(exon_stat)
    if len(i_exon_pairs) == 0:
        return []
    mut_num = 1
    muts = []
    for pair in i_exon_pairs:
        first_exon = pair[0]
        second_exon = pair[1]
        # in codon table numbers are 0-based, so correct
        c_f_exon = first_exon - 1
        c_s_exon = second_exon - 1
        # get split codon for first exon
        # there are two exons, N and M, betweend them - deleted guys
        # exon_NX -------- Xexon_M
        # X marks split codons that potentially contain stop
        # to take last codon of N, if it's split, I should take 0'st codon for exon N+1
        # if it's split -> it goes to N + 1 exon with "split" != 0 field
        # for exon M -> just take the 0'st codon
        try:
            f_ex_split = [c for c in codon_table if c["t_exon_num"] == c_f_exon + 1][0]
            s_ex_split = [c for c in codon_table if c["t_exon_num"] == c_s_exon][0]
        except IndexError:
            # case of ultrashort exons (1-2bp)
            # better to avoid any conclusions
            continue
        if f_ex_split["split_"] == 0 or s_ex_split["split_"] == 0:
            # one of those codons is complete -> split stop is impossible
            continue
        # cut corresponding seq
        f_ex_seq = f_ex_split["que_codon"][:f_ex_split["split_"]]
        s_ex_seq = s_ex_split["que_codon"][s_ex_split["split_"]:]
        # f_part_len = len(f_ex_seq)
        # s_part_len = len(s_ex_seq)
        split_codon_seq = f_ex_seq + s_ex_seq
        split_triplets = parts(split_codon_seq, 3)
        stops_in = [x for x in split_triplets if x in STOPS]
        if len(stops_in) == 0:
            # no stops on split
            continue
        # there are stop codons!
        mut_ = stops_in[0]
        # ex_num = second_exon if s_part_len > f_part_len else first_exon
        mut_id = f"SP_STOP_{mut_num}"
        mut_num += 1
        mut = Mutation(gene=gene, chain=q_name, exon=second_exon, position=0,
                        mclass=STOP, mut=mut_, masked=False, mut_id=mut_id)
        muts.append(mut)
    return muts


def get_out_of_borders_prop(codon_table, miss_exons):
    """Compute a proportion of out-of-chain-borders bases."""
    gene_len = len(codon_table)
    m_codons_len = 0
    for m_exon in miss_exons:
        m_codons_len += len([c for c in codon_table if c["t_exon_num"] == m_exon])
    if gene_len == 0:
        return 0.0
    m_prop = m_codons_len / gene_len
    return m_prop


def inact_mut_check(cesar_data, u12_introns=None, v=False, gene="None",
                    ex_prop=None, ref_ss=None, sec_codons=None, no_fpi=False):
    """Importable function."""
    cesar_fractions = read_cesar_out(cesar_data)\
    # TODO: optimise this part
    u12_introns_data = parse_u12_opt(gene, u12_introns)
    if ref_ss:  # fake U12 introns, still the same result
        u12_introns_data = u12_introns_data.union(ref_ss)
    # intact_percent = {}
    p_intact_ignore_M = {}
    p_intact_intact_M = {}
    middle_80_intact = {}
    middle_80_present = {}
    i_codons_prop = {}
    out_of_b_vals = {}
    mutations = []

    for cesar_fraction in cesar_fractions:
        # get data for SSM
        fraction_mutations = []
        q_name = cesar_fraction[0]
        q_name_d_key = int(q_name) if q_name.isnumeric() else q_name
        ref = cesar_fraction[1]
        query = cesar_fraction[2]
        if v:
            eprint(f"Detecting inactivating mutations for query: {q_name}")
        # try to get properties dicts
        # chain_to_exon_to_properties = (chain_exon_class, chain_exon_gap, pIDs, pBl, chain_missed)
        if ex_prop is None:
            exon_class = {}
            exon_gap = {}
            exon_pid = {}
            exon_blosum = {}
            missing_exons = {}
            ex_inc = {}
            ex_lens = {}
        else:
            exon_class = ex_prop[0].get(q_name_d_key, {})
            exon_gap = ex_prop[1].get(q_name_d_key, {})
            exon_pid = ex_prop[2].get(q_name_d_key, {})
            exon_blosum = ex_prop[3].get(q_name_d_key, {})
            missing_exons = ex_prop[4].get(q_name_d_key, set())
            ex_inc = ex_prop[5].get(q_name_d_key, {})
            ex_lens = ex_prop[6]

        # report.append(f"###ANALYSING QUERY {q_name}")
        # one loop -> for splice site mutations
        sps_mutations = analyse_splice_sites(ref, query, gene, q_name, u12_introns_data, v=v)
        fraction_mutations.extend(sps_mutations)
        # another -> for stops and FS
        codon_table = parse_cesar_out(ref, query)

        # next loop -> for deleted/missed exons
        if ex_prop:
            exon_del_miss_, exon_stat_ = classify_exons(gene, q_name, codon_table, exon_class,
                                                        exon_gap, exon_pid, exon_blosum, missing_exons,
                                                        ex_inc)
            exon_del_miss, exon_stat = find_safe_ex_dels(exon_del_miss_, exon_stat_, ex_lens, no_fpi=no_fpi)
            fraction_mutations.extend(exon_del_miss)
        else:
            # we don't have a lot of data
            # will do it outselves
            exon_stat = None
            pass
        big_indel_thrs = infer_big_indel_thresholds(ex_lens)
        inact_muts = scan_rf(codon_table,
                             gene,
                             q_name,
                             exon_stat=exon_stat,
                             v=v,
                             big_indel_thrs=big_indel_thrs,
                             sec_codons=sec_codons,
                             no_fpi=no_fpi)
        fraction_mutations.extend(inact_muts)
        split_stop_codons = detect_split_stops(codon_table, gene, q_name, exon_stat)
        fraction_mutations.extend(split_stop_codons)

        compensations = detect_compensations(inact_muts, codon_table)
        fraction_mutations.extend(compensations)
        fraction_mutations = mask_compensated_fs(fraction_mutations)
        fraction_mutations = filter_mutations(fraction_mutations)

        pintact_features = compute_intact_perc(codon_table, fraction_mutations, q_name, v=v)
        p_intact_ignore_M[q_name] = pintact_features[0]
        p_intact_intact_M[q_name] = pintact_features[1]
        i_codons_prop[q_name] = pintact_features[2]
        middle_80_intact[q_name] = pintact_features[3]
        middle_80_present[q_name] = pintact_features[4]

        out_of_borders_prop = get_out_of_borders_prop(codon_table, missing_exons)
        out_of_b_vals[q_name] = out_of_borders_prop

        mutations.extend(fraction_mutations)

    report = muts_to_text(mutations, p_intact_ignore_M, p_intact_intact_M, i_codons_prop,
                          out_of_b_vals, middle_80_intact, middle_80_present, gene)
    return report


def main():
    """Entry point of a standalone script."""
    args = parse_args()
    f = open(args.cesar_output, "r")
    cesar_line = f.read()
    f.close()
    report = inact_mut_check(cesar_line, u12_introns=args.u12, gene=args.gene, v=args.verbose)
    print(report)
    sys.exit(0)


if __name__ == "__main__":
    main()
