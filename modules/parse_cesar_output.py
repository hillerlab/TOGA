#!/usr/bin/env python3
"""Parse CESAR output for one query, get codon table."""
import argparse
import sys
from copy import deepcopy

STOPS = {"TAG", "TGA", "TAA"}
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


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_cesar_out(target, query, v=False):
    """Convert raw CESAR output into a codon table."""
    # CESAR output parsing assumes that target sequence starts with a space
    # but in 0.1% cases it is not true, so for those cases I'll as 1 space:
    target = " " + target
    query = "-" + query
    # get codons number
    letters_num = len([c for c in target if c.isalpha()])
    codons_num = letters_num // 3
    eprint(f"Expecting {codons_num} target codons") if v else None
    codon_data = []

    next_elem_box = {"ref_codon": "",  # keep reference codon; str
                     "que_codon": "",  # keep query codon; str
                     "q_exon_num": 0,  # exon number of codon, in query; int
                     "t_exon_num": 0,  # exon num in target; int
                     "split_": 0,  # 0 - not split 
                     "que_coords": []  # relative coordinates in query
                     }
    for _ in range(codons_num):
        # fill codon table with struct-like objects
        codon_data.append(deepcopy(next_elem_box))

    was_space = False
    intr_del_switch = False
    exon_num = -1  # to start with 0 for consistence
    t_exon_num = -1
    codon_num = 0
    codon_counter = 0
    is_split_now = False
    q_coord = -1
    
    for t, q in zip(target, query):
        if t == ">":
            if not intr_del_switch:
                t_exon_num += 1
            intr_del_switch = True
            continue
        else:
            intr_del_switch = False
        if q != "-":
            q_coord += 1
        if t == " " and not was_space:
            # if space -> intron
            # new exon starts
            was_space = True
            exon_num += 1
            t_exon_num += 1
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
                codon_data[codon_num]["split_"] = len(codon_data[codon_num]["ref_codon"])
                # assert codon_data[codon_num]["split_"] in {0, 1, 2}
            continue
        elif t == " ":
            continue
        else:  # not a space
            was_space = False
        # t is not a space
        codon_data[codon_num]["ref_codon"] += t
        codon_data[codon_num]["que_codon"] += q
        codon_data[codon_num]["q_exon_num"] = exon_num
        codon_data[codon_num]["t_exon_num"] = t_exon_num
        # print(codon_num, codon_data[codon_num])
        if q != "-":
            codon_data[codon_num]["que_coords"].append(q_coord)
        # split/ non split flags
        if t.isupper():
            # increase number of items in codon
            is_split_now = False
        elif t.islower():
            # we are going throw a split codon
            is_split_now = True
        elif t == "-":
            continue
        else:  # should never happen
            eprint(f"Broken codon:{str(codon_data[codon_num])}") if v else None
            die(f"CESAR output is corrupted - char {t} in query seq", 1)
        # not a gap -> switch codon counters
        codon_counter += 1
        if codon_counter == 3:
            codon_counter = 0
            codon_num += 1
        if codon_num == codons_num:
            # must be no codons anymore
            break
        
    for elem in codon_data:
        eprint(str(elem)) if v else None

    # check and finalize
    for i in range(codons_num):
        i_codon = codon_data[i]
        is_first = i == 0
        is_last = i == codons_num - 1
        t_codon_seq = i_codon["ref_codon"].replace("-", "")
        if t_codon_seq.upper() != "ATG" and is_first:
            eprint("Error! CESAR output is corrupted, target must start with ATG!")
        elif t_codon_seq.upper() not in STOPS and is_last:
            eprint("Error! CESAR output is corrupted, target must end with a stop codon!")
        elif t_codon_seq.upper() in STOPS and not is_last:
            eprint("Error! CESAR output is corrupted, found in-frame STOP codon in reference!")
        all_hi = all(x.isupper() for x in t_codon_seq)
        all_lo = all(x.islower() for x in t_codon_seq)
        if not all_hi and not all_lo:
            eprint(f"Broken codon:{str(i_codon)}")
            eprint("Error! CESAR output is corrupted, wrong split codon mapping!")
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


def classify_exon(ex_class, incl, pid, blosum):
    """Decide what do we do with this exon."""
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
    # TODO: make a better output format
    for elem in codons_data:
        print(elem)

if __name__ == "__main__":
    main()
