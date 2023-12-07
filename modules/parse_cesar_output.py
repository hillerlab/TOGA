#!/usr/bin/env python3
"""Parse raw CESAR output for one query, get codon table."""
import argparse
import sys
from copy import deepcopy
from version import __version__

try:  # for robustness
    from modules.common import eprint
    from modules.common import die
except ImportError:
    from common import eprint
    from common import die


# constants
STOPS = {"TAG", "TGA", "TAA"}
# nucleotide %ID and blosum thresholds
# for exons classification
HQ_PID = 75
HQ_BLOSUM = 65

AB_INCL_PID = 25
AB_INCL_BLOSUM = 25

C_INCL_PID = 65
C_INCL_BLOSUM = 50

A_T_PID = 65
A_T_BLOSUM = 40

LO_T_PID = 45
LO_T_BLOSUM = 25


def parse_cesar_out(target, query, v=False):
    """Convert raw CESAR output into a codon table."""
    # CESAR output structure notes:
    # >reference
    #       ATGGCAa            aaGTC>>>CTGGGGAtt          cCCC...
    # >query
    # atcagcATGGGAAagtacgtagcgtAAGTCCCCCTACCGATAaggatcgtgtCCCC...
    # This is a pairwise alignment between the reference and the query
    # CESAR aligns reference exon on the query sequence
    # In reference:
    #    space means -> no exon aligns, either an intron or intergenic region
    #    uppercase letter -> complete codon
    #    lowercase letter -> codon split between two exons
    #    > sign -> marks intron deletion
    #    - -> insertion in query
    # In query:
    #    uppercase letter -> coding sequence (where reference exons align)
    #    lowercase letter -> non-coding sequence (introns, UTR, intergenic regions)
    #    - -> deletion in query
    # CESAR accepts intact coding reference exons as input

    # CESAR output parsing assumes that target sequence starts with a space
    # but in 0.1% cases it is not true, so for those cases I'll add 1 space:
    target = " " + target
    query = "-" + query
    # get codons number
    letters_num = len([c for c in target if c.isalpha()])
    codons_num = letters_num // 3
    eprint(f"Expecting {codons_num} target codons") if v else None
    codon_data = []
    # this "struct" contains codon information
    next_elem_box = {
        "ref_codon": "",  # keep reference codon; str
        "que_codon": "",  # keep query codon; str
        "q_exon_num": 0,  # exon number of codon, in query; int
        "t_exon_num": 0,  # exon num in target; int
        "split_": 0,  # 0 - not split
        "que_coords": [],  # relative coordinates in query
    }
    for _ in range(codons_num):
        # fill codon table with struct-like objects
        codon_data.append(deepcopy(next_elem_box))

    # initiate some values before parsing
    was_space = False  # previous char in reference is space
    intr_del_switch = False  # previous char in reference was >
    exon_num = -1  # (in query) to start with 0 for consistence
    t_exon_num = -1  # exon num in reference, initial value
    codon_num = 0  # codon num, start with 0
    codon_counter = 0  # counter for codons
    is_split_now = False  # flag means that we are reading a split codon
    q_coord = -1  # coordinate in query (relative)

    for t, q in zip(target, query):
        # read a pair of characters, one from ref, another from query
        # intron deletion-related conditions
        if t == ">":
            # going through an intron deletion
            if not intr_del_switch:
                # first char in this deletion
                # exon number in reference changes, in query -> does not
                t_exon_num += 1
            # flag ON -> not to increase ref exon number in this >>>> run
            intr_del_switch = True
        else:
            # not the intron deletion: put intron_del flag off
            intr_del_switch = False

        if q != "-":
            # not a gap in the query: next letter -> inc coordinate
            q_coord += 1

        if t == " " and not was_space:
            # if space -> intron
            # new intron starts
            was_space = True
            exon_num += 1  # exon just ended in both ref and query
            t_exon_num += 1  # so increment the numbers
            # fill the codons
            curr_codon = codon_data[codon_num]
            curr_codon_letters = [c for c in curr_codon["ref_codon"] if c.isalpha()]
            curr_codon_gaps = [c for c in curr_codon["ref_codon"] if c == "-"]

            if len(curr_codon_gaps) > 0 and len(curr_codon_letters) == 0:
                # empty codon terminates here
                codon_counter = 0
                codons_num += 1
                codon_num += 1
                codon_data.append(deepcopy(next_elem_box))

            if is_split_now:
                # this is a split codon: compute split_ value
                codon_data[codon_num]["split_"] = len(
                    codon_data[codon_num]["ref_codon"]
                )
            continue

        elif t == " ":
            # space continues (was_space == True)
            continue
        else:  # not a space, intron flag off
            was_space = False

        # t is not a space if we are here
        # add codon letters per codon
        t_fixed = t if t != ">" else "-"  # in case it's > we would like to add a -
        codon_data[codon_num]["ref_codon"] += t_fixed
        codon_data[codon_num]["que_codon"] += q
        # update exon numbers
        codon_data[codon_num]["q_exon_num"] = exon_num
        codon_data[codon_num]["t_exon_num"] = t_exon_num

        if q != "-":  # it not gap -> a letter -> a valid query coordinate
            codon_data[codon_num]["que_coords"].append(q_coord)

        # split/ non split flags
        if t.isupper():
            # increase number of items in codon
            is_split_now = False
        elif t.islower():
            # we are going through a split codon
            is_split_now = True
        elif t == "-":
            # change nothing
            continue
        elif t == ">":
            # the same as a gap
            continue
        else:  # should never happen, unexpected character in the reference seq
            eprint(f"Broken codon:{str(codon_data[codon_num])}") if v else None
            die(f"CESAR output is corrupted - char {t} in query seq", 1)
        # not a gap -> switch codon counters
        codon_counter += 1
        if codon_counter == 3:
            # we got 3 letters in the reference codon: goto next one
            codon_counter = 0
            codon_num += 1
        if codon_num == codons_num:
            # must be no codons anymore
            break

    for elem in codon_data:
        # show codon table if required
        eprint(str(elem)) if v else None

    # check and finalize
    for i in range(codons_num):
        i_codon = codon_data[i]
        is_first = i == 0
        is_last = i == codons_num - 1
        t_codon_seq = i_codon["ref_codon"].replace("-", "")
        # sanity checks, write error messages if something is wrong with input
        if t_codon_seq.upper() != "ATG" and is_first:  # starts with ATG?
            eprint("Error! CESAR output is corrupted, target must start with ATG!")
        elif t_codon_seq.upper() not in STOPS and is_last:  # ends with STOP?
            eprint(
                "Error! CESAR output is corrupted, target must end with a stop codon!"
            )
        elif t_codon_seq.upper() in STOPS and not is_last:  # Stop in frame?
            eprint(
                "Error! CESAR output is corrupted, found in-frame STOP codon in reference!"
            )
        # all ref letters in codon must be either lower or uppercase, not a mixture
        all_hi = all(x.isupper() for x in t_codon_seq)
        all_lo = all(x.islower() for x in t_codon_seq)
        if not all_hi and not all_lo:
            eprint(f"Broken codon:{str(i_codon)}")
            eprint("Error! CESAR output is corrupted, wrong split codon mapping!")
        # make all uppercase
        codon_data[i]["ref_codon"] = codon_data[i]["ref_codon"].upper()
    return codon_data


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("cesar_out", help="CESAR output file, this script input")
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def classify_exon(
            ex_class: str, incl: bool, pid: float, blosum: float
        ) -> tuple[bool, str]:
    """Decide what do we do with this exon."""
    # ex_class: how exon aligns
    # A - chain aligns exon perfectly
    # B - an exon flank (left or right) is not aligned
    # C - chain blocks don't intersect the exon
    # M - exon located outside the chain borders
    # we look at:
    # 1) exon class
    # 2) nucleotide %ID and BLOSUM scores
    # 3) was exon detected in the expected region?
    # and then classify it or say it's deleted/missing

    # A / B class branch
    if ex_class == "A" or ex_class == "B" or ex_class == "A+":
        if incl:
            # CHECK FOR A+ CLASS
            if pid > A_T_PID and blosum > A_T_BLOSUM:
                return True, "HQ"
                # ok, if class A, True, HQ -> check for A+
            else:
                return True, "AQ"
        else:
            # as well as C included
            if pid > A_T_PID and blosum > A_T_BLOSUM:
                return True, "AQ"
            elif pid > LO_T_PID and blosum > LO_T_BLOSUM:
                return True, "LQ"
            else:
                return False, "NA"
    # class C branch
    elif ex_class == "C":
        if incl:
            # the same with A/B class excluded
            if pid > A_T_PID and blosum > A_T_BLOSUM:
                return True, "AQ"
            elif pid > LO_T_PID and blosum > LO_T_BLOSUM:
                return True, "LQ"
            else:
                return False, "NA"
        else:
            # class C excluded -> we do not annotate this a priori
            return False, "NA"
    # class M branch
    else:
        # class M -> no classification actually
        return False, "NA"


def main():
    """Entry point."""
    args = parse_args()
    f = open(args.cesar_out, "r")
    lines = f.readlines()
    f.close()
    ref_seq = lines[1]
    que_seq = lines[3]
    codons_data = parse_cesar_out(ref_seq, que_seq)
    for elem in codons_data:
        print(elem)


if __name__ == "__main__":
    main()
