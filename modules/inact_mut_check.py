#!/usr/bin/env python3
"""Scan reading frame for inactivating mutations."""
import sys
import argparse
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import asdict
from constants import Constants
from constants import InactMutClassesConst as MutClasses
from modules.parse_cesar_output import (parse_cesar_out)
from modules.parse_cesar_output import classify_exon
from modules.common import die
from modules.common import eprint
from modules.common import parts
# please look for all named constants in the
# modules/GLP_values.py  # TODO: revise it later
from modules.GLP_values import *

__author__ = "Bogdan Kirilenko, 2020."
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


@dataclass
class Mutation:
    gene: str  # gene, chain, exon -> exon identifier
    chain: str
    exon: int
    position: int  # num of codon, 0 if entire exon missed/deleted
    mclass: str  # MutClasses.SSM, FS_DEL, etc.
    mut: str  # if stop -> stop codon, if INS/DEL -> ins/del size, etc.
    masked: bool  # if in first/last 10% -> ignore
    mut_id: str  # mutation ID


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
    # each four lines stand for one CESAR output unit
    # one CESAR unit -> one pairwise alignment
    cesar_content = cesar_line.split("\n")
    fractions = parts(cesar_content, 4)
    cesar_fractions = []

    for fraction in fractions:
        # parse fractions
        # we need sequences and query sequence names
        # to distinguish then
        if len(fraction) < 4:
            # the last line of the file
            continue
        # fraction[0] -> reference sequence name, we don't need it
        ref_seq = fraction[1]
        query_name = fraction[2][1:]  # remove > - fasta file
        query_seq = fraction[3]

        if len(ref_seq) != len(query_seq):
            # ref and query seq must have the same length in the pairwise alignment
            die("Error! Ref and query sequences must have the same length!")
        elif len(ref_seq) == 0:
            # also must never happen -> there is an error
            die("Error! The input is empty!")
        # save fraction data
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


def create_masked_mut(mut):
    """To refactor later.

    Previous implementation used named tuple which is immutable.
    Dataclass is mutable, however, the logic should be still adjusted.
    """
    return Mutation(
        gene=mut.gene,
        chain=mut.chain,
        exon=mut.exon,
        position=mut.position,
        mclass=mut.mclass,
        mut=mut.mut,
        masked=True,
        mut_id=mut.mut_id,
    )


def analyse_splice_sites(
    ref,
    query,
    gene,
    chain,
    codon_table,
    atg_codons_data,
    mask_all_first_10p=False,
    u12_data=None,
    v=None,
):
    """Check correctness of the splice sites."""
    u12_data = set() if not u12_data else u12_data  # load U12 data if provided
    mut_counter = 1  # mut counter -> for mutations IDs
    eprint(f"U12 introns: {u12_data}") if v else None

    sps_report = []  # save mut data here
    # CESAR alignment looks like this:
    # >reference
    #       ATGGCAa            aaGTC>>>CTGGGGAtt          cCCC...
    # >query
    # atcagcATGGGAAagtacgtagcgtAAGTCCCCCTACCGATAaggatcgtgtCCCC...
    # CESAR aligns coding exons only
    # spaces in reference sequence: nothing, not exon at all
    # > in the reference mean intron deletion
    # in ref: uppercase letters - full codons, lowercase - split codons
    # in query: uppercase letters - CDS, lowercase - something else
    # so we can extract indexes of CDS in the reference - not spaces and >
    cds_indexes = [i for i, c in enumerate(ref) if c != " " and c != ">"]
    exon_num = 0  # initiate exons counter
    # if letter indexes follow each other such as 1,2,3,4 -> they belong to same exon
    # it there is a gap, as 1,2,3,10,11,12 -> there are different exons
    exon_num_indexes = defaultdict(list)

    for i in range(1, len(cds_indexes)):
        # go letter-by-letter, keep current and previous letter
        prev = cds_indexes[i - 1]
        current = cds_indexes[i]
        delta = current - prev
        exon_num_indexes[exon_num].append(prev)
        # if delta curr - prev > 1: curr and prev are on different exons
        # increment exon number then
        if delta > 1:
            exon_num += 1

    # for each exon get min and max index
    exon_num_start_end = {k: (min(v), max(v)) for k, v in exon_num_indexes.items()}
    exons_num = len(exon_num_start_end)

    for exon_num_, start_end in exon_num_start_end.items():
        # exon_num_ is 0-based
        # but we used to use 1-based exon numbers
        exon_num = exon_num_ + 1
        start, end = start_end
        # we know start and end of the exon
        # so we can get splice site coordinates
        # if there is N in the splice site -> we don't know what's there
        # -> we are not sure -> mask this mutation
        acceptor_splice_site = query[start - 2: start]  # donor
        acceptor_splice_site_N = (
            "n" in acceptor_splice_site or "N" in acceptor_splice_site
        )
        donor_splice_site = query[end + 1: end + 3]  # acceptor
        donor_splice_size_N = "n" in donor_splice_site or "N" in donor_splice_site

        # assign to the last / first codon of the exon (depends on what splice site is affected)
        codon_pos_at_exon = [
            n for n, c in enumerate(codon_table) if c["t_exon_num"] == exon_num_
        ]
        acceptor_codon_num = (
            codon_pos_at_exon[0] if len(codon_pos_at_exon) > 0 else None
        )
        donor_codon_num = codon_pos_at_exon[-1] if len(codon_pos_at_exon) > 0 else None

        # check that splice sites are canonical
        acceptor_site_wrong = acceptor_splice_site not in LEFT_SPLICE_CORR
        donor_site_wrong = donor_splice_site not in RIGHT_SPLICE_CORR

        eprint(
            f"Exon {exon_num}; L_SPS: {acceptor_splice_site}; R_SPS: {donor_splice_site}"
        ) if v else None

        # print(exon_num)
        # print(acceptor_site_wrong)
        # print(donor_site_wrong)

        if exon_num != 1 and acceptor_site_wrong and acceptor_codon_num:
            # add mutation for left (acceptor) splice site
            # doesn't apply to the first exon obviously
            # mask this mutation if it's suspected to be U12 splice site
            mask = True if (exon_num, 0) in u12_data else False
            # if N in the splice site -> also mask it
            mask = True if acceptor_splice_site_N else mask

            # if splice site mutation in first 10% but followed by ATG in first 10% then
            atg_mask = _define_whether_mask(
                acceptor_codon_num,
                atg_codons_data["left_t"],
                atg_codons_data["right_t"],
                atg_codons_data["atg_codon_nums"],
                mask_all_first_10p=mask_all_first_10p,
            )

            mask = True if atg_mask is True else mask

            # create mutation object, describe what happened
            mut_ = f"{LEFT_SPLICE_CORR}->{acceptor_splice_site}"
            mut_id = f"{MutClasses.SSM_A}_{mut_counter}"
            # print(mut_id)
            mut_counter += 1
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=exon_num,
                position=acceptor_codon_num,
                mclass=MutClasses.SSM_A,
                mut=mut_,
                masked=mask,
                mut_id=mut_id,
            )
            sps_report.append(mut)  # add mutation to the list

        if exon_num != exons_num and donor_site_wrong and donor_codon_num:
            # add mutation for right (donor) splice site
            # doesn't apply to the last exon
            # mask this mutation if it's suspected to be U12 splice site
            mask = True if (exon_num, 1) in u12_data else False
            mask = True if donor_splice_size_N else mask  # if N in mutation -> mask it

            # if splice site mutation in first 10% but followed by ATG in first 10% then
            atg_mask = _define_whether_mask(
                donor_codon_num,
                atg_codons_data["left_t"],
                atg_codons_data["right_t"],
                atg_codons_data["atg_codon_nums"],
                mask_all_first_10p=mask_all_first_10p,
            )

            mask = True if atg_mask is True else mask

            # create mutation object
            mut_ = f"{RIGHT_SPLICE_CORR}->{donor_splice_site}"
            mut_id = f"{MutClasses.SSM_D}_{mut_counter}"
            # print(mut_id)
            mut_counter += 1
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=exon_num,
                position=donor_codon_num,
                mclass=MutClasses.SSM_D,
                mut=mut_,
                masked=mask,
                mut_id=mut_id,
            )
            sps_report.append(mut)  # add it to the list
    if len(sps_report) <= 1:
        # 0 or 1 mutations: return them
        return sps_report

    # if more than 1 -> exclude intron deletions
    sps_report = sorted(sps_report, key=lambda x: (x.exon, x.position))

    # sort mutations from the beginning to the end
    for i in range(1, len(sps_report)):
        # iterate over pairs of mutations: current and previous
        j = i - 1
        prev = sps_report[j]
        curr = sps_report[i]
        curr_exon = curr.exon
        prev_exon = prev.exon

        if curr_exon != prev_exon + 1:
            # if current exon doesn't follow the previous immediately -> not the case
            # like prev mut exon is 3 and current is 6
            continue

        # if exons follow each other (like 3 and 4) then continue
        curr_pos = curr.position
        prev_pos = prev.position
        curr_type = curr.mclass
        prev_type = prev.mclass

        # they must belong to the same intron, check this
        if not (prev_type == MutClasses.SSM_D and curr_type == MutClasses.SSM_A):
            continue

        prev_to_what = prev.mut.split("->")[1]
        curr_to_what = curr.mut.split("->")[1]
        # if it was -- in both cases: this is intron deletion
        # mask these mutations
        if prev_to_what == curr_to_what == "--":
            # intron deletion
            sps_report[j].masked = True
            sps_report[i].masked = True
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
        letters_only = [c for c in exon_seq if c.isalpha()]
        ex_len = len(letters_only) - prev_rem
        full_codons = ex_len // 3
        rem = ex_len % 3
        prev_rem = 0 if rem == 0 else 3 - rem
        codons_num = full_codons if rem == 0 else full_codons + 1
        for _ in range(codons_num):
            codon_to_exon[curr_codon] = num
            curr_codon += 1
    return codon_to_exon


def corr_exon_num_or_no_fs(codon, exon_num):
    """Need to decide which exon num to assign."""
    # count letters in the left and rigth exon
    cut_at = codon["split_"]
    left_side_ref = codon["ref_codon"][:cut_at]
    left_side_que = codon["que_codon"][:cut_at]
    # count gaps on each side
    ls_ref_gaps = left_side_ref.count("-")
    ls_que_gaps = left_side_que.count("-")
    # if number of gaps on the left exon % 3 != 0:
    # we say the left exon if mutated (exon_num - 1)
    # otherwise the rigth one (exon_num)
    fs_left = abs(ls_ref_gaps - ls_que_gaps) % 3 != 0
    if fs_left:
        return exon_num - 1
    else:
        return exon_num


def _find_atg_codons(codon_table):
    """Find reference codons aligned to ATG."""
    atg_codon_nums = []
    for num, codon in enumerate(codon_table, 1):
        que_codon = codon["que_codon"]
        que_codon_no_gap = que_codon.replace("-", "")
        if not que_codon_no_gap:
            # there are only gaps -> nothing to catch
            continue
        triplets = parts(que_codon_no_gap, n=3)
        start_triplets = [x for x in triplets if x == "ATG"]
        if len(start_triplets) > 0:
            atg_codon_nums.append(num)
    return atg_codon_nums


def _get_next_bigger_num(num, lst):
    for elem in lst:
        if elem >= num:
            return elem
    return 999999999  # TODO: ideally, last codon position


def _define_whether_mask(
    num, left_t, right_t, atg_codon_nums, mask_all_first_10p=False
):
    """Check whether the mutation is going to be masked due to first/last 10% or not."""
    # TODO: optimise this part, the only ATG position needed is the closest to 10%
    # which is not above 10%, that's it.
    if num >= right_t:
        return True
    elif left_t < num < right_t:
        return False
    if num <= left_t and mask_all_first_10p is True:
        # automatically mask mut in first 10%
        # don't account for ATG codons distribution
        return True
    # num in first 10%, need to find the next start
    next_atg_pos = _get_next_bigger_num(num, atg_codon_nums)
    return next_atg_pos <= left_t


def make_atg_data(codon_table):
    codons_num = len(codon_table)

    # get first or last 10% if CDS: mask mutations in this region:
    perc_10 = codons_num // 10
    left_t = perc_10
    right_t = codons_num - perc_10

    atg_codon_nums = _find_atg_codons(codon_table)

    # TODO: refactor this
    atg_codon_nums_data = {
        "left_t": left_t,
        "right_t": right_t,
        "atg_codon_nums": atg_codon_nums,
    }
    return atg_codon_nums_data


def scan_rf(
    codon_table,
    gene,
    chain,
    atg_codon_nums_data,
    exon_stat=None,
    v=False,
    big_indel_thrs=None,
    sec_codons=None,
    no_fpi=False,
    mask_all_first_10p=False,
):
    """Scan codon table for inactivating mutations."""
    # sec_codons -> selenocysteine-coding codons in reference
    sec_codons_set = sec_codons if sec_codons else set()
    codons_num = len(codon_table)
    in_mut_report = []  # save mutations here

    # init mutation counters: for IDs
    mut_number_stop = 1
    mut_number_atg = 1  # not inact, but need to track
    mut_number_fs = 1
    mut_number_big_indel = 1
    q_dels_in_a_row = 0  # number of deletions in a row: for big indel detection

    for num, codon in enumerate(codon_table, 1):
        # go codon-by-codon
        # if in first/last 10%: it will be masked
        # mask = True if num <= left_t or num >= right_t else False
        mask = _define_whether_mask(
            num,
            atg_codon_nums_data["left_t"],
            atg_codon_nums_data["right_t"],
            atg_codon_nums_data["atg_codon_nums"],
            mask_all_first_10p=mask_all_first_10p,
        )
        # determine whether it's the first or last exon:
        last_codon = True if num == codons_num else False
        first_codon = True if num == 1 else False

        ex_num = codon["t_exon_num"] + 1  # need to be 1-based exon num
        ref_codon = codon["ref_codon"]  # ref and query codon sequences
        que_codon = codon["que_codon"]
        codon_is_split = codon["split_"] > 0  # if > 0 - codon split between exons
        # it depends on the exon length:
        big_indel_thr = big_indel_thrs[ex_num] if big_indel_thrs else BIG_INDEL_SIZE

        if que_codon == "---" and codon_is_split is False:
            # count number of deletions in a row
            q_dels_in_a_row += 1
        else:
            # not a deletion, but it's possible that a sequence of deletions
            # is stopped -> need to check how big it was
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
                mut = Mutation(
                    gene=gene,
                    chain=chain,
                    exon=ex_num,
                    position=position,
                    mclass=mclass,
                    mut=mut_,
                    masked=mask,
                    mut_id=mut_id,
                )
                in_mut_report.append(mut)
            # anyway update counter of deleted codons
            q_dels_in_a_row = 0

        if first_codon and que_codon != "ATG":
            # start codon is missing; we don't use it in the GLP pipe
            # but still detect this sort of mutations
            mut_id = "START_1"
            que_codon_ = que_codon[:3]
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=ex_num,
                position=1,
                mclass=START_MISSING,
                mut=que_codon_,
                # this mutation should not affect the classification
                # so it's always masked, but saved to indicate that
                # 1st ATG is missing:
                masked=True,
                mut_id=mut_id,
            )
            in_mut_report.append(mut)

        # check for FS; count number of dashes in both codons
        ref_gaps = ref_codon.count("-")
        que_gaps = que_codon.count("-")
        delta = ref_gaps - que_gaps
        fs = abs(delta) % 3 != 0  # if so, it's a frameshift

        if fs:
            # we have a frame-shifting indel!
            mut_ = f"+{delta}" if delta > 0 else f"{delta}"
            mclass = FS_INS if delta > 0 else FS_DEL
            mut_id = f"FS_{mut_number_fs}"
            mut_number_fs += 1
            if not codon_is_split:
                # if FS occur in non-split codon: exon num determination is simple
                fs_ex_num = ex_num
            else:
                # FS happened in a split codon: need another procedure
                # to decide, in which one
                fs_ex_num = corr_exon_num_or_no_fs(codon, ex_num)
            # add mutation object to the list
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=fs_ex_num,
                position=num,
                mclass=mclass,
                mut=mut_,
                masked=mask,
                mut_id=mut_id,
            )
            in_mut_report.append(mut)
            if v:  # verbose
                eprint("Detected FS")
                eprint(codon)
            # we didn't stop detecting inactivating mutations in this codon
            # so we can find them (one codon might have FS together with inframe-stop)
            # mask = True to avoid counting these codons as inactivated twice
            mask = True

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
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=bi_ex_num,
                position=num,
                mclass=mclass,
                mut=mut_,
                masked=mask,
                mut_id=mut_id,
            )
            in_mut_report.append(mut)
            mask = True  # again, to prevent counting affected codons twice

        # check for inframe stop codons
        que_codon_no_gap = que_codon.replace("-", "")
        if not que_codon_no_gap:
            # there are only gaps -> nothing to catch
            continue
        triplets = parts(que_codon_no_gap, n=3)
        # one codon objects corresponds to a single reference codon
        # however, there migth be several codons in query
        # need to split query sequence in triplets
        # check that any of them is a stop-codon
        stop_triplets = [x for x in triplets if x in Constants.STOP_CODONS]
        start_triplets = [x for x in triplets if x == "ATG"]

        if len(stop_triplets) > 0 and not last_codon:
            # we have premature stop codon
            mut_ = f"{ref_codon.replace('-', '')}->{stop_triplets[0]}"
            mclass = STOP
            mut_id = f"STOP_{mut_number_stop}"
            # need to detect which exon is affected
            if not codon_is_split:
                # not split, easy to detect which exon is affected
                st_ex_num = ex_num
            else:
                # split stop codon: might be difficult
                prev_exon = ex_num - 1
                prev_exon_stat = exon_stat[prev_exon] if exon_stat else "I"
                this_exon_stat = exon_stat[ex_num] if exon_stat else "I"
                if this_exon_stat in D_M or prev_exon_stat in D_M:
                    # one of the exons is deleted: don't do anything
                    # might be a false signal
                    continue
                # codon consists of 3 letters
                # assign to exon that has 2 of them:
                st_ex_num = ex_num if codon["split_"] == 1 else ex_num - 1
            mut_number_stop += 1
            # check whether it's a selenocysteine-coding codon
            is_sec_pos = num - 1 in sec_codons_set  # 0-based in that set
            if is_sec_pos and stop_triplets[0] == "TGA":
                # mask if U-coding codon
                stop_mask = True
            else:  # apply usual rules
                stop_mask = mask
            # create mutation object, add ths to list
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=st_ex_num,
                position=num,
                mclass=mclass,
                mut=mut_,
                masked=stop_mask,
                mut_id=mut_id,
            )
            in_mut_report.append(mut)

            if v:
                eprint("Detected STOP")
                eprint(codon)
        if len(start_triplets) > 0:
            # not an inactivating mutation but needs to be saved
            mut_ = f"{ref_codon.replace('-', '')}->ATG"
            mclass = ATG
            mut_id = f"ATG_{mut_number_atg}"
            st_ex_num = ex_num
            mut_number_atg += 1
            mut = Mutation(
                gene=gene,
                chain=chain,
                exon=st_ex_num,
                position=num,
                mclass=mclass,
                mut=mut_,
                masked=True,
                mut_id=mut_id,
            )
            in_mut_report.append(mut)
            if v:
                eprint("Detected ATG (non-inactivating)")
                eprint(codon)
    # I can infer number of codons in each exon directly from codon table
    return in_mut_report


def detect_compensations(inact_mut, codon_table):
    """Detect FS compensation events.

    We call frameshifts compensated if:
    1) After a compensated FS the original reading frame is preserved.
    Such as: +2 and +1 in result give us +3.
    2) There is no stop codon in the alternative frame.
    """
    m_counter = 1  # for mut IDs
    answer = []  # collect compensated events here
    fs = [mut for mut in inact_mut if mut.mclass in {FS_DEL, FS_INS}]
    fs_num = len(fs)
    if fs_num <= 1:
        # need at least 2 frameshifs, otherwise there is nothing to compensate
        return answer, []
    # first iteration: detect potential compensatory events, based only
    # on their sizes
    potent_compensations = []
    for i_num in range(fs_num - 1):
        # get compensated runs of FS
        # pick them one by one and check whether the following FS
        # can potentially compensate them
        # skip the last mutation, it cannot start compensation
        init_mut = fs[i_num]
        init_mut_val = int(init_mut.mut)
        # use mut ID to determine compensated FS
        init_mut_id = init_mut.mut_id
        # create two lists for this run: mut sizes and IDs
        fs_values = [
            init_mut_val,
        ]
        fs_ids = [
            init_mut_id,
        ]
        # iter over next mutations
        for j_num in range(i_num + 1, fs_num):
            j_mut = fs[j_num]
            j_mut_val = int(j_mut.mut)
            j_mut_id = j_mut.mut_id
            # append size and ID to this run list
            fs_values.append(j_mut_val)
            fs_ids.append(j_mut_id)
            # check whether it's potential compensation
            if sum(fs_values) % 3 == 0:
                potent_compensations.append(fs_ids)
                break
    if len(potent_compensations) == 0:
        return [], []  # no potential compensations, skip this

    # verify potential compensations, check for stop codons in alt frame
    what_is_compensated = set()  # to avoid twice compensated FS
    alt_frame_codons = []  # lenghts of codons in alt frame

    for comp in potent_compensations:
        # comp -> a list of potentially compensated FS IDs
        if len(what_is_compensated.intersection(comp)) > 0:
            # it means that some of these FS are already compensated
            continue
        # get mutations itself using their IDs
        comp_muts = sorted(
            [m for m in fs if m.mut_id in comp], key=lambda x: x.position
        )
        start_pos = comp_muts[0].position
        end_pos = comp_muts[-1].position
        # positions are 1-based, need to correct
        # get alt frame sequence
        codons_seq = codon_table[start_pos - 1: end_pos]
        que_seq = "".join([c["que_codon"] for c in codons_seq]).replace("-", "")
        # split this sequence in codons:
        upd_codons = parts(que_seq, n=3)
        if len(Constants.STOP_CODONS.intersection(upd_codons)) > 0:
            # there are stops in the compensated sequence
            continue
        # add compensation track
        # no stop codons
        # to create a mut object we need gene name, chain id, exon num etc
        alt_frame_codons.append((start_pos - 1, end_pos))
        ethalon_mut = comp_muts[0]
        gene = ethalon_mut.gene
        chain = ethalon_mut.chain
        exon = ethalon_mut.exon  # not applicable really, but let it be
        position = ethalon_mut.position
        mclass = MutClasses.COMPENSATION
        mut_id = f"C_{m_counter}"
        # mutation -> comma-separated list of compensated mutation IDs
        fs_ids = [x.split("_")[1] for x in comp]
        # mut = "FS_{}".format(",".join(fs_ids))
        # changed format FEB 2022: FS_{start_num}-{end_num}
        # just a comma-separated list can be too long
        _ids_range = f"{fs_ids[0]}-{fs_ids[-1]}"
        mut = f"FS_{_ids_range}"
        comp_mut = Mutation(
            gene=gene,
            chain=chain,
            exon=exon,
            position=position,
            mclass=mclass,
            mut=mut,
            masked=False,
            mut_id=mut_id,
        )
        answer.append(comp_mut)
        # calculieren positions of first and last alt frame codons
        for c in comp:  # add compensated muts to comp muts set
            # to avoid adding comp mutations twice
            what_is_compensated.add(c)
        m_counter += 1
    return answer, alt_frame_codons


def mask_compensated_fs(mut_list):
    """Mask compensated mutations."""
    comp_muts = [m for m in mut_list if m.mclass == MutClasses.COMPENSATION]
    if len(comp_muts) == 0:
        # no compensations, no worries
        return mut_list
    comp_fs_ids = []  # there are compensatory events
    for cmp in comp_muts:
        # get list of comp FS IDs
        # fs_nums = cmp.mut.split("_")[1].split(",")
        comp_ids_range_str = cmp.mut.split("_")[1].split("-")
        _comp_start = int(comp_ids_range_str[0])
        _comp_end = int(comp_ids_range_str[1])
        fs_nums = list(range(_comp_start, _comp_end + 1))
        fs_ids = [f"FS_{x}" for x in fs_nums]
        comp_fs_ids.extend(fs_ids)
    comp_fs_ids = set(comp_fs_ids)
    # not touch -> other mutations, not FS and not compensated
    not_touch = [m for m in mut_list if m.mut_id not in comp_fs_ids]
    # to_mask -> these we'd like to mask
    to_mask = [m for m in mut_list if m.mut_id in comp_fs_ids]
    masked = [create_masked_mut(m) for m in to_mask]
    filtered = not_touch + masked
    return filtered


def compute_percent_id(seq_1, seq_2):
    """Return % identity for two sequences."""
    assert len(seq_1) == len(seq_2)  # otherwise it is a bug
    matches = sum(
        [
            1
            for i in range(len(seq_1))
            if seq_1[i] == seq_2[i]
            and seq_1[i] != "N"
            and seq_2[i] != "N"
            and seq_1[i] != "-"
            and seq_2[i] != "-"
        ]
    )
    length = len(seq_1.replace("N", "").replace("-", ""))  # we ignore N's
    pid = matches * 100 / length if length != 0 else 0
    return pid


def _get_last_codon_for_each_exon(codon_table):
    ret = {}
    for num, elem in enumerate(codon_table, 1):
        ret[elem["t_exon_num"] + 1] = num
    return ret


def classify_exons(
    gene,
    que,
    codon_table,
    exon_class,
    exon_gap,
    exon_pid,
    exon_blosum,
    missing_exons,
    ex_inc,
    atg_codons_data,
    mask_all_first_10p=False,
    v=False,
):
    """Classify exons as intact, deleted and missing."""
    del_num, miss_num = 1, 1  # counters for mutation IDs
    # get a list of exon numbers:
    exon_nums = list(range(codon_table[-1]["t_exon_num"] + 1))
    exon_to_last_codon_of_exon = _get_last_codon_for_each_exon(codon_table)
    exons_report = []  # save data heve
    exon_stat = [
        "X",
    ]  # exon status, start with 1, X - placeholder

    del_miss_nums = []

    for exon_num in exon_nums:
        # 0-based to 1-based
        ex_num_ = exon_num + 1
        if exon_num in missing_exons:
            print(f"Exon num {exon_num} in missing exons list") if v else None
            # a priori missing
            # missing exons list is 0-based
            # exon_num - 0 based, exon_num_ (with underscore) is 1-based
            mut = Mutation(
                gene=gene,
                chain=que,
                exon=ex_num_,
                position=0,
                mut="-",
                mclass=MutClasses.MISS_EXON,
                masked=False,
                mut_id=f"MIS_{miss_num}",
            )
            miss_num += 1
            exons_report.append(mut)
            exon_stat.append("M")
            del_miss_nums.append(exon_num)
            continue
        # parse data from CESAR wrapper output
        ex_class = exon_class.get(exon_num, None)  # exon class
        ex_gap = exon_gap.get(exon_num, None)  # intersect an asm gap
        ex_pid = exon_pid.get(exon_num, 0)  # nucleotide %ID
        ex_blosum = exon_blosum.get(exon_num, 0)  # blosum score
        exon_excl = ex_inc.get(exon_num, None)  # detected outside expected region
        # classify whether it's deleted or not:
        # print(f"Exon {exon_num} classification with the following params: ") if v else None
        # print(ex_class, exon_excl, ex_pid, ex_blosum) if v else None
        ex_non_del, q = classify_exon(ex_class, exon_excl, ex_pid, ex_blosum)
        # print(f"Results are: {del_} {q}") if v else None

        if ex_class == "M" or ex_gap:
            # if intersects assembly gap or M: write a mutation
            mut = Mutation(
                gene=gene,
                chain=que,
                exon=ex_num_,
                position=0,
                mut="-",
                mclass=MutClasses.MISS_EXON,
                masked=False,
                mut_id=f"MIS_{miss_num}",
            )
            miss_num += 1
            # add to mut list, add new exon status
            exons_report.append(mut)
            exon_stat.append("M")
            del_miss_nums.append(exon_num)
            continue
        elif ex_non_del is False:
            # exon is deleted: need to write about this
            last_codon_num = exon_to_last_codon_of_exon.get(ex_num_, 0)

            atg_mask = _define_whether_mask(
                last_codon_num,
                atg_codons_data["left_t"],
                atg_codons_data["right_t"],
                atg_codons_data["atg_codon_nums"],
                mask_all_first_10p=mask_all_first_10p,
            )

            mut = Mutation(
                gene=gene,
                chain=que,
                exon=ex_num_,
                position=0,
                mut="-",
                mclass=MutClasses.DEL_EXON,
                masked=atg_mask,
                mut_id=f"DEL_{del_num}",
            )
            del_num += 1
            # add to mut list, append new exon status
            exons_report.append(mut)
            exon_stat.append("D")
            del_miss_nums.append(exon_num)
            continue
        else:
            # something else -> exon is not deleted
            exon_stat.append("I")  # add I status to exon
            pass
    # return list of mutation objects + list of exon statuses
    return exons_report, exon_stat, del_miss_nums


def muts_to_text(
    mutations, perc_intact_1, perc_intact_2, i_prop, oub, m_80_i, m_80_p, gene
):
    """Convert mutations array to text."""
    ordered = sorted(mutations, key=lambda x: (x.chain, x.exon, x.position))
    strings = []  # save string representations here
    # first, save all inactivating mutations:
    for elem in ordered:
        # convert named tuple to list:
        elem_values = list(asdict(elem).values())
        # last field masked values False and True might be confusing
        # make explicit "masked" and "unmasked"
        if elem_values[6] is False:
            elem_values[6] = "unmasked"
        else:
            elem_values[6] = "masked"
        elem_values = [str(x) for x in elem_values]
        # hashtag to distinguish CESAR output from GeneLossScanner
        elem_string = "# " + "\t".join(elem_values)
        strings.append(elem_string)
    # then same %intact-related features, just fill the template
    # again, hashtag to distinguish from CESAR output
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
        m_80_line = f"# {gene}\t{q}\tMIDDLE_IS_INTACT {val_str}"
        strings.append(m_80_line)
    for q, val in m_80_p.items():
        val_str = "TRUE" if val is True else "FALSE"
        m_80_line = f"# {gene}\t{q}\tMIDDLE_IS_PRESENT {val_str}"
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


def compute_intact_perc(
    codon_table,
    mutations,
    q_name,
    alt_frame_ranges,
    alt_f_del=False,
    v=False,
    mask_all_first_10p=False,
):
    """Compute intact %ID-related features."""
    # compute per query
    query_muts = [m for m in mutations if m.chain == q_name]
    gene_len = len(codon_table)  # num of codons in reference
    # initiate codon_status, mark deleted codons with D, the rest with I
    codon_status = ["I" if c["que_codon"] != "---" else "D" for c in codon_table]
    # get numbers of deleted/missing exons
    del_exons = {
        m.exon - 1 for m in query_muts if m.mclass == MutClasses.DEL_EXON and m.masked is False
    }
    safe_del_exons = {
        m.exon - 1 for m in query_muts if m.mclass == MutClasses.DEL_EXON and m.masked is True
    }
    miss_exons = {m.exon - 1 for m in query_muts if m.mclass == MutClasses.MISS_EXON}
    # using this data, get numbers of codons in missing/deleted exons
    del_codon_nums = [
        n for n, c in enumerate(codon_table) if c["t_exon_num"] in del_exons
    ]
    # read codons in alt frame as deleted

    safe_del_codon_nums = [
        n for n, c in enumerate(codon_table) if c["t_exon_num"] in safe_del_exons
    ]
    miss_codon_nums = [
        n for n, c in enumerate(codon_table) if c["t_exon_num"] in miss_exons
    ]

    if alt_f_del is True:
        # if so, we consider codons in alternative frame deleted
        # alt frame -> between compensated frameshifts
        for elem in alt_frame_ranges:
            s_, e_ = elem
            for codon_num in range(s_, e_):
                safe_del_codon_nums.append(codon_num)

    # update codons status -> what is Missing, Deleted or Lost (not-safely deleted exons)
    for del_codon in del_codon_nums:
        codon_status[del_codon] = "L"
    for del_codon in safe_del_codon_nums:
        codon_status[del_codon] = "D"
    for miss_codon in miss_codon_nums:
        codon_status[miss_codon] = "M"

    # get IDs of compensated FS
    compensations = [m for m in query_muts if m.mclass == MutClasses.COMPENSATION]
    comp_fs = []  # list to keep compensated frameshifts
    for comp in compensations:
        comp_field = comp.mut
        comp_ids_range_str = comp_field.split("_")[1].split("-")
        # # fmt: FS_{start}-{end}
        _comp_start = int(comp_ids_range_str[0])
        _comp_end = int(comp_ids_range_str[1])
        comp_ids = list(range(_comp_start, _comp_end + 1))
        comp_fs_strings = [f"FS_{i}" for i in comp_ids]
        comp_fs.append(comp_fs_strings)

    # put inactivating mutations coordinates in codon status table
    for m in query_muts:
        # print(m)
        # go mutation-by-mutation
        if m.mut_id in comp_fs:
            # if compensated FS -> does not affect %intact
            continue
        elif m.mclass in MutClasses.DEL_MISS:
            # exon missing or deletion -> already considered
            continue
        elif m.mclass == MutClasses.COMPENSATION:
            # compensation is also an event in this list -> skip
            continue
        elif m.masked is True:
            continue
        elif m.mclass == MutClasses.SSM_A or m.mclass == MutClasses.SSM_D:
            # deal with splice site mutations
            if m.masked is True:
                # U12 or N-containing splice site -> do not consider
                continue
            to_what = m.mut.split("->")[1]  # what changed
            if to_what == "??" or to_what.upper() == "NN":
                # we don't know what happened there
                # I left this condition for compatibility with old GLP
                continue
            # assign to the last / first codon of the exon (depends on what splice site is affected)
            ssm_exon = m.exon - 1  # switch to 0-based now
            # get list of codon numbers of this exon
            codon_pos_at_exon = [
                n for n, c in enumerate(codon_table) if c["t_exon_num"] == ssm_exon
            ]
            if len(codon_pos_at_exon) == 0 and ssm_exon == 0:
                # very short (<3bp) 1st exon affected: simpler to assign 0
                affected_num = 0
            elif len(codon_pos_at_exon) == 0:
                # VERY strange situation
                # the exon is entirely absent -> seems to be really corrupted output
                continue  # better to skip this
            elif m.position == 0:  # left side, first codon of this exon
                affected_num = codon_pos_at_exon[0]
            else:  # right side: last codon in this exon
                affected_num = codon_pos_at_exon[-1]
            # if not Missing codon: say it's Lost (update codon status)
            # seems to be a redundant check for M (I filter those mutations)
            if codon_status[affected_num] != "M":
                codon_status[affected_num] = "L"
            continue
        elif m.mclass == START_MISSING:
            # start missing -> do not participate
            continue
        # the rest of mutations: mark corresponding codon
        affected_num = m.position - 1
        # update codon status if it's not missing
        if codon_status[affected_num] != "M":
            codon_status[affected_num] = "L"

    if all(x == "I" for x in codon_status):
        # nearly impossible, all codons are intact
        return 1.0, 1.0, 1.0, True, True

    if v:
        eprint(f"Num of codons: {len(codon_status)}")
        eprint(f"Num of I: {codon_status.count('I')}")
        eprint(f"Num of D: {codon_status.count('D')}")
        eprint(f"Num of M: {codon_status.count('M')}")
        eprint(f"Num of L: {codon_status.count('L')}")

    # simpler to convert this codon status list to string
    codon_status_string = "".join(codon_status)
    # compute features
    # remove deleted codons from computations
    repl_d = codon_status_string.replace("D", "")
    # to compute %intact in two modes:
    ignore_m = repl_d.replace("M", "")  # 1. - ignoring missing sequence
    m_intact = repl_d.replace("M", "I")  # 2. - considering missing sequence as intact

    # get spans of continuous non-lost reading frame
    # For example:
    # L-I-I-I-I-I-I-I-L-L-L-I-I-I-I
    # spans are: I-I-I-I-I-I-I & I-I-I-I
    ignore_m_spans = ignore_m.split("L")  # ignoring M sequence
    m_intact_spans = m_intact.split("L")  # considering M sequence is intact

    # get the longest span of non-lost sequence
    ignore_m_max_span = max(len(x) for x in ignore_m_spans)  # ignoring M sequence
    m_intact_max_span = max(
        len(x) for x in m_intact_spans
    )  # considering M sequence is intact

    # compute %intact: two modes again
    p_intact_ignore_m = ignore_m_max_span / gene_len  # ignoring M sequence
    p_intact_intact_m = m_intact_max_span / gene_len  # considering M sequence is intact

    # compute features related to middle 80% CDS
    ten_perc = gene_len // 10
    # cut middle 80% codons states:
    middle = codon_status_string[ten_perc:-ten_perc]
    first_90 = codon_status_string[:-ten_perc]

    if mask_all_first_10p is True:
        middle_is_intact = False if "L" in middle else True
    else:
        middle_is_intact = False if "L" in first_90 else True
    middle_is_present = False if "M" in middle else True

    # compute non-missing sequence
    not_m = len(codon_status) - codon_status.count("M")
    if not_m > 0:  # beware of zero division error
        num_of_I_codons = codon_status.count("I") / not_m
    else:  # meaning all codons are missing
        num_of_I_codons = 1.0
    # return everything computed
    return (
        p_intact_ignore_m,
        p_intact_intact_m,
        num_of_I_codons,
        middle_is_intact,
        middle_is_present,
    )


def filter_mutations(mut_list):
    """Remove mutations of deleted or missed exons."""
    # get mutations of deleted and missing exons
    # select exon deletions and missing events itself
    # we don't want to mask them, only the rest
    del_exon_muts = [m for m in mut_list if m.mclass == MutClasses.DEL_EXON]
    mis_exon_muts = [m for m in mut_list if m.mclass == MutClasses.MISS_EXON]
    del_exons = [m.exon for m in del_exon_muts]
    missing_exons = [m.exon for m in mis_exon_muts]
    del_and_missed = set(missing_exons + del_exons)
    # mut to keep: mutations that are NOT in the del/missing exons
    mut_to_keep = [m for m in mut_list if m.exon not in del_and_missed]
    # so we keep: mutations which state that some exons are M/D
    # + mutations that affect other exons, not M/D
    # mutations in the Del/Miss exons are omitted
    filtered = mut_to_keep + del_exon_muts + mis_exon_muts
    return filtered


def get_d_runs(ex_stat):
    """Get D runs."""
    d_inds = [n for n, v in enumerate(ex_stat) if v == "D" or v == "mD"]
    if len(d_inds) <= 1:
        # nothing we can do
        return []
    result = []
    curr_list = [
        d_inds[0],
    ]
    for elem in d_inds[1:]:
        prev_elem = curr_list[-1]
        if prev_elem + 1 == elem:
            # append to the current list
            curr_list.append(elem)
        else:
            # start new list
            curr_copy = curr_list.copy()
            result.append(curr_copy)
            curr_list = [
                elem,
            ]
    result.append(curr_list)
    final_res = [x for x in result if len(x) > 1]
    return final_res


def find_safe_ex_dels(mut_list, ex_stat_, ex_lens, no_fpi=False):
    """Select safe exon deletions.

    Safe exon deletions:
    1) If exon length % 3 == 0: frame is preserved after this deletion.
    2) N exons deleted in a row, sum of their lengths % 3 == 0: the same.
    We call case (2) compensated exon deletion: as compensated frameshifts.
    3) Deletion of first or last exon (if it's too short)
    We redefine them as Missing exons.
    """
    if ex_lens is None:
        # no exon lengths -> no data for decision
        return mut_list, ex_stat_
    upd_mut_list = []  # save updated mutations here
    mdel_num = 1
    compensated_ex = set()
    # get compensated exon deletions
    # at first get runs of exon deletions in a row:
    D_runs = get_d_runs(ex_stat_)
    for D_run in D_runs:
        # in D run we can detect a compensation
        if len(D_run) < 2:
            # single ex del cannot compensate itself
            continue
        # bulky solution, get lengths of deleted exon first
        D_lens = [ex_lens[d] for d in D_run]
        D_run_len = len(D_run)  # num of exons deleted

        for i_num in range(D_run_len - 1):
            # N^2 algorithm, start with each exon
            # try to find whether next exons can compensate this
            d_elem = D_run[i_num]
            if d_elem in compensated_ex:
                # already compensated: cannot be compensated twice
                continue
            d_len = D_lens[i_num]
            # fill these lists with next exons:
            d_elems = [
                d_elem,
            ]
            d_lens = [
                d_len,
            ]
            # maybe there is a bit more elegant solution?
            for j_num in range(i_num + 1, D_run_len):
                # add next exons one-by-one
                j_elem = D_run[j_num]
                j_len = D_lens[j_num]
                # add their lengths and numbers
                d_elems.append(j_elem)
                d_lens.append(j_len)
                if sum(d_lens) % 3 != 0:
                    # not compensated yet: we can continue
                    continue
                # they are compensated: add them to the compensated set
                for x in d_elems:
                    compensated_ex.add(x)
                # compensated: break J look, back to I loop
                break

    # filter exon D/M mutations
    # deal with small first/last exon deletions
    for m in mut_list:
        ex_len = ex_lens[m.exon]
        is_first = m.exon == 1
        is_last = m.exon == len(ex_lens)
        last_or_first = is_first or is_last
        frame_pres = ex_len % 3 == 0 or m.exon in compensated_ex
        # if we ignore FP indels at all: then any exon size fits, always true
        # no fpi is True -> frame-pres indels not ignored -> ignore deletions of short exons only
        # no fpi is False -> ignore all frame-preserving exon deletions
        fp_ex_len_cond = ex_len < SAFE_EXON_DEL_SIZE if no_fpi is True else True

        if last_or_first and ex_len < FIRST_LAST_DEL_SIZE:
            # then we say it's missing, short last or first exon
            # create a new mutation object
            new_m = Mutation(
                gene=m.gene,
                chain=m.chain,
                exon=m.exon,
                position=m.position,
                mclass=MutClasses.MISS_EXON,
                mut_id=f"MDEL_{mdel_num}",
                mut=m.mut,
                masked=m.masked,
            )
            mdel_num += 1
            upd_mut_list.append(new_m)
            ex_stat_[m.exon] = "M"
            continue
        elif frame_pres and fp_ex_len_cond:
            # exon del is frame-preserving and short
            masked_m = create_masked_mut(m)  # mask this mutation then
            upd_mut_list.append(masked_m)
            ex_stat_[m.exon] = "mD"  # mD -> minor deletion
        else:
            # just add the rest of mutations
            upd_mut_list.append(m)
    # return filtered exon Del/Mis list + updated exon status
    return upd_mut_list, ex_stat_


def infer_big_indel_thresholds(ex_lens):
    """For each exon define big indel threshold."""
    ex_T = {}
    if ex_lens is None:
        # no exon lengths data -> no result
        return {}
    for ex_num, ex_len in ex_lens.items():
        # if exon is small enough: use standard 40bp
        if ex_len <= BIG_EXON_THR:
            ex_T[ex_num] = BIG_INDEL_SIZE
            continue
        # exon is quite long, require 20% of it's length
        # to call indel inactivating
        thr = int(ex_len / 5)
        ex_T[ex_num] = thr
    return ex_T


def get_exon_pairs(exon_stat):
    """Get pairs of I exons separated by Deleted."""
    # for example, it we have exon start like this:
    # exon_num:  0-1-2-3-4-5-6-7-8-9-10-11
    # exon_stat: X-I-D-D-I-I-D-I-I-M-D---I
    # then we need the following output:
    # (1, 4), (5, 7). Pair (8, 11) is skipped
    # because there is a M exon between them.
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


def detect_split_stops(
    codon_table, gene, q_name, exon_stat, atg_codons_data, mask_all_first_10p=False
):
    """Considering all exon deletions find all split stop codons."""
    # we need to get pairs of intact (not deleted or missing) exons
    # between what there is a row of Deleted exons (but not missing)
    # if there is a missing exon in between -> do not make any conclusions
    # split stop codons may occur between these pairs
    i_exon_pairs = get_exon_pairs(exon_stat)
    exon_to_last_codon_of_exon = _get_last_codon_for_each_exon(codon_table)

    if len(i_exon_pairs) == 0:
        # no such pairs -> no way to detect split stop
        return []  # return nothing

    mut_num = 1
    muts = []
    for pair in i_exon_pairs:
        # those are 1-based
        first_exon = pair[0]
        second_exon = pair[1]
        # in codon table numbers are 0-based, so correct
        c_f_exon = first_exon - 1
        c_s_exon = second_exon - 1
        # get split codon for first exon
        # there are two exons, N and M, between them - deleted guys
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
        f_ex_seq = f_ex_split["que_codon"][: f_ex_split["split_"]]
        s_ex_seq = s_ex_split["que_codon"][s_ex_split["split_"]:]
        split_codon_seq = f_ex_seq + s_ex_seq
        split_triplets = parts(split_codon_seq, 3)
        stops_in = [x for x in split_triplets if x in Constants.STOP_CODONS]
        # Maybe add ATGs too?
        if len(stops_in) == 0:
            # no stops on split
            continue
        # there are stop codons!
        mut_ = stops_in[0]
        # ex_num = second_exon if s_part_len > f_part_len else first_exon
        mut_id = f"SP_STOP_{mut_num}"
        mut_num += 1
        position = exon_to_last_codon_of_exon[first_exon]
        atg_mask = _define_whether_mask(
            position,
            atg_codons_data["left_t"],
            atg_codons_data["right_t"],
            atg_codons_data["atg_codon_nums"],
            mask_all_first_10p=mask_all_first_10p,
        )

        mut = Mutation(
            gene=gene,
            chain=q_name,
            exon=second_exon,
            position=position,
            mclass=STOP,
            mut=mut_,
            masked=atg_mask,
            mut_id=mut_id,
        )
        muts.append(mut)
    return muts


def get_out_of_borders_prop(codon_table, miss_exons):
    """Compute a proportion of out-of-chain-borders bases."""
    gene_len = len(codon_table)
    m_codons_len = 0
    for m_exon in miss_exons:
        # ++ lengths of missing exons
        m_codons_len += len([c for c in codon_table if c["t_exon_num"] == m_exon])
    if gene_len == 0:
        # to avoid zerodivision error
        return 0.0
    m_prop = m_codons_len / gene_len
    return m_prop


def inact_mut_check(
    cesar_data,
    u12_introns=None,
    v=False,
    gene="None",
    ex_prop=None,
    ref_ss=None,
    sec_codons=None,
    no_fpi=False,
    alt_f_del=False,
    mask_all_first_10p=False,  # mask all inact mutations in 1st 10% regardless on ATGs
):
    """Detect inactivating mutations in the CESAR output."""
    # read cesar output
    # note that CESAR accepts only one set of reference exons but
    # also it can process multiple query sequences
    # so we need check them one-by-one
    cesar_fractions = read_cesar_out(cesar_data)
    # parse U12 data
    u12_introns_data = parse_u12_opt(gene, u12_introns)
    if ref_ss:
        # non-canonical splice sites in the reference
        # process them as U12 splice sites -> the same logic
        u12_introns_data = u12_introns_data.union(ref_ss)

    # for each projection (CESAR unit) we extract 6
    # features + list of inactivating mutations (if they exist)
    # initiate lists/ dicts for them
    p_intact_ignore_M = {}
    p_intact_intact_M = {}
    middle_is_intact = {}
    middle_is_present = {}
    i_codons_prop = {}
    out_of_b_vals = {}
    mutations = []
    del_miss_exons = {}

    for cesar_fraction in cesar_fractions:
        # analyze cesar fractions one-by-one
        fraction_mutations = []  # save inact mutations for this fraction here
        # cesar_fraction: (query_name, ref_sequence, query_sequence)
        q_name = cesar_fraction[0]  # need to distinguish with other queries
        # if called by TOGA: q_name is numeric (basically just chainID)
        q_name_d_key = int(q_name) if q_name.lstrip("-").isdigit() else q_name
        ref = cesar_fraction[1]
        query = cesar_fraction[2]

        if v:
            eprint(
                f"Detecting inactivating mutations for query: {q_name}  ({q_name_d_key})"
            )
            eprint(
                f"Types of q_name: {type(q_name)}/ of q_name_d_key: {type(q_name_d_key)}"
            )

        # parse additional information provided by CESAR wrapper
        # chain_to_exon_to_properties = (chain_exon_class, chain_exon_gap, pIDs, pBl, chain_missed)
        if ex_prop is None:
            # if not provided: create empty dicts
            exon_class = {}
            exon_gap = {}
            exon_pid = {}
            exon_blosum = {}
            missing_exons = {}
            ex_inc = {}
            ex_lens = {}
        else:
            # if provided: extract data related to this query
            exon_class = ex_prop[0].get(q_name_d_key, {})
            exon_gap = ex_prop[1].get(q_name_d_key, {})
            exon_pid = ex_prop[2].get(q_name_d_key, {})
            exon_blosum = ex_prop[3].get(q_name_d_key, {})
            missing_exons = ex_prop[4].get(q_name_d_key, set())
            ex_inc = ex_prop[5].get(q_name_d_key, {})
            ex_lens = ex_prop[6]

        #  create codon table to extract mutations
        # codon table: list of objects, describing a codon
        # such as sequence in reference and query, exon number and so on
        codon_table = parse_cesar_out(ref, query)
        atg_codons_data = make_atg_data(codon_table)

        # next loop -> for deleted/missed exons
        if ex_prop:  # if extra data provided by CESAR wrapper we can classify exons
            # dm list -> list of 0-based exon nums which are del or missing
            exon_del_miss_, exon_stat_, dm_list = classify_exons(
                gene,
                q_name,
                codon_table,
                exon_class,
                exon_gap,
                exon_pid,
                exon_blosum,
                missing_exons,
                ex_inc,
                atg_codons_data,
                mask_all_first_10p=mask_all_first_10p,
                v=v,
            )
            # get lists of deleted/missing exons
            # also find "safe" exon deletions: if they are in-frame
            # or a series of exon deletions in a row is frame-preserving
            exon_del_miss, exon_stat = find_safe_ex_dels(
                exon_del_miss_, exon_stat_, ex_lens, no_fpi=no_fpi
            )
            fraction_mutations.extend(exon_del_miss)  # also add this
        else:
            # we don't have extra exons data
            # will just skip this part
            exon_stat = None
            dm_list = []
            pass

        # big indels may be classified as inactivating mutations
        # but the bigger the exon: the bigger an indel should be
        # define the thresholds
        del_miss_exons[q_name] = set(dm_list)
        big_indel_thrs = infer_big_indel_thresholds(ex_lens)

        # scan reading frame (codon table) for the rest on inact mutations
        inact_muts = scan_rf(
            codon_table,
            gene,
            q_name,
            atg_codons_data,
            exon_stat=exon_stat,
            v=v,
            big_indel_thrs=big_indel_thrs,
            sec_codons=sec_codons,
            no_fpi=no_fpi,
            mask_all_first_10p=mask_all_first_10p,
        )
        # save this data
        fraction_mutations.extend(inact_muts)

        # now we extract inactivation mutations
        # then add them to fraction_mutations list
        # extract splice site mutations
        sps_mutations = analyse_splice_sites(
            ref,
            query,
            gene,
            q_name,
            codon_table,
            atg_codons_data,
            mask_all_first_10p,
            u12_introns_data,
            v=v,
        )
        fraction_mutations.extend(sps_mutations)

        # get a list of split stop codons: stop codons that appear after exon deletions
        # such as:
        # GCAAACGCAGCt-------------[DELETED EXON]-------agTCCCATTTCCAACTGATC
        # exon deletion raises an inframe stop codon: t + ag
        split_stop_codons = detect_split_stops(
            codon_table,
            gene,
            q_name,
            exon_stat,
            atg_codons_data,
            mask_all_first_10p=mask_all_first_10p,
        )
        fraction_mutations.extend(split_stop_codons)

        # detect compensated frameshifts
        compensations, alt_frame_ranges = detect_compensations(inact_muts, codon_table)
        fraction_mutations.extend(compensations)
        # also mask compensated frameshifts:
        fraction_mutations = mask_compensated_fs(fraction_mutations)
        # filter inactivating events
        # for example, if an exon is deleted: remove all frameshifts,
        # stop codons etc associated with this exon
        fraction_mutations = filter_mutations(fraction_mutations)

        # extract and save %intact-related features
        pintact_features = compute_intact_perc(
            codon_table,
            fraction_mutations,
            q_name,
            alt_frame_ranges,
            alt_f_del=alt_f_del,
            v=v,
            mask_all_first_10p=mask_all_first_10p,
        )
        p_intact_ignore_M[q_name] = pintact_features[0]
        p_intact_intact_M[q_name] = pintact_features[1]
        i_codons_prop[q_name] = pintact_features[2]
        middle_is_intact[q_name] = pintact_features[3]
        middle_is_present[q_name] = pintact_features[4]

        # compute %of gene that lies outside chain borders
        out_of_borders_prop = get_out_of_borders_prop(codon_table, missing_exons)
        out_of_b_vals[q_name] = out_of_borders_prop
        # add it to a overall list of inactivating mutations
        mutations.extend(fraction_mutations)
    # create a string that could be saved to a file (or stdout)
    report = muts_to_text(
        mutations,
        p_intact_ignore_M,
        p_intact_intact_M,
        i_codons_prop,
        out_of_b_vals,
        middle_is_intact,
        middle_is_present,
        gene,
    )
    return report, del_miss_exons


def main():
    """Entry point of a standalone script."""
    args = parse_args()
    f = open(args.cesar_output, "r")
    # inact_mut_check() accepts a string containing raw CESAR output
    cesar_line = f.read()
    f.close()
    report, _ = inact_mut_check(
        cesar_line, u12_introns=args.u12, gene=args.gene, v=args.verbose
    )
    print(report)
    sys.exit(0)


if __name__ == "__main__":
    main()
