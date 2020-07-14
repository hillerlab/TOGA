#!/usr/bin/env python3
"""Fasta oneLine <-> N, trim, rm sequences, build tree."""
import argparse
import sys
import os
import re
import subprocess
from collections import Counter, defaultdict

__author__ = 'Bogdan Kirilenko, 2018'

# genetic code if translation is needed
nta = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
       "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
       "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
       "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
       "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
       "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
       "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
       "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
       "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
       "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
       "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
       "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
       "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
       "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
       "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
       "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
       "---": "-", "NNN": "X"}
compl = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "-": "-"}


def die(msg):
    """Write a message and die."""
    sys.stderr.write(msg + "\n")
    sys.stderr.write("Program finished with exit code 1.\n")
    sys.exit(1)


def parts(lst, n=25):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def test_reachable(path):
    """Test if a pathway exists."""
    try:  # open and close the target file
        test_open = open(path, "w")
        test_open.close()
    except FileNotFoundError:  # catch exception if the path is unreachable
        die("Path {0} is unreachable!".format(path))
    except IsADirectoryError:
        die("Path {0} is a directory! Need a text file.".format(path))


def read_fasta(fasta_file, show_headers=False, rm=""):
    """Read fasta, return dict and type."""
    # open the file
    input_stream = open(fasta_file, "r") if fasta_file != "stdin" else sys.stdin
    # with open(input_stream, "r") as f:
    fasta_data = input_stream.read().split(">")
    input_stream.close()
    # TODO check if fasta file is correct
    assert fasta_data[0] == ""  # if a file starts with > it should be empty
    del fasta_data[0]  # remove it "" we don't need that
    sequences = {}  # accumulate data here
    order = []  # to have ordered list
    # read line by line

    to_rm = rm.split(",")  # make removal list

    for elem in fasta_data:
        raw_lines = elem.split("\n")
        header = raw_lines[0]  # it must be first ['capHir1', 'ATGCCGCGCCAATTCCCCAAGCTGA... ]
        if header in to_rm:
            continue  # do not add if we don't need it
        lines = [x for x in raw_lines[1:] if x != ""]  # separate nucleotide-containing lines
        if len(lines) == 0:  # it is a mistake - empty sequene --> get rid of
            continue
        fasta_content = "".join(lines)
        sequences[header] = fasta_content
        order.append(header)
    if len(sequences) == 0:
        die("There are not fasta-formatted sequences in {0}!".format(fasta_file))
    if len(sequences.keys()) != len(order):  # it is possible in case of non-unique headers
        err = "Error! Sequences names must be unique! There are" \
              " {0} sequences and {1} unique names!\n".format(len(sequences.keys()), len(order))
        intersect = [k for k, v in Counter(order).items() if v > 1]
        err += "Intersections are:\n{0}".format(",".join(intersect))
        die(err)
    if show_headers:  # just print all the >'s and interrupt
        sys.stdout.write(",".join(order) + "\n")
        sys.exit(0)

    return sequences, order


def read_phylip(input_file, show_headers=False, rm=""):
    """Read phylip format."""
    f = open(input_file, "r")
    block_counter = 0
    seq_name_blocks = defaultdict(list)
    seq_num_name = {}
    skip_header = False
    in_counter = 0
    order = []
    for line in f:
        line_data = line.split()
        if block_counter == 0 and skip_header is False:
            skip_header = True
            continue

        if line_data == []:
            block_counter += 1
            in_counter = 0
            continue

        if block_counter == 0:
            seq_name = line_data[0]
            seq_start = line_data[1]
            order.append(seq_name)
            seq_num_name[in_counter] = seq_name
            seq_name_blocks[seq_name].append(seq_start)
            in_counter += 1
            continue

        if block_counter > 0:
            seq_block = line_data[0]
            seq_name = seq_num_name.get(in_counter)
            seq_name_blocks[seq_name].append(seq_block)
            in_counter += 1
            continue
    data = {}
    for name, seq_blocks in seq_name_blocks.items():
        seq = "".join(seq_blocks)
        data[name] = seq
    return data, order


def copy_paste(data, order, copy, paste):
    """Copy sequence with name copy to sequence with name to paste."""
    if copy not in order:  # nothing to copy
        die("Error! There is no sequence {0} to copy. Use one of:\n{1}".format(copy, " ".join(order)))
    if paste in order:  # not forbidden but may be dangerous
        sys.stderr.write("Warning! Sequence {0} to --paste exists. Rewriting...\n".format(paste))
    if copy == paste:
        sys.stderr.write("Warning! --copy and --paste are the same sequence.")
    order.append(paste) if paste not in order else None
    data[paste] = data[copy]
    return data, order


def get_branches(tree, only_term=True):
    """Return a list of branches."""
    # remove all special symbols
    filtered_str = re.sub(r'[^\w]', ' ', tree)
    # get rid of numbers
    filtered = set([x for x in filtered_str.split() if not x.isdigit()])
    filtered = [x for x in filtered if "-" not in x and "_" not in x] if only_term else filtered
    return filtered


def build_tree(elems, output, und, not_anc=False):
    """Build tree using a list of elements given, save to output."""
    # locate tree file
    tree_cwd = os.path.dirname(__file__)  # tree in the same dir with the file
    tree_file = "all_phylo.tree"  # MAIN tree file
    tree_path = os.path.join(tree_cwd, "data", tree_file)
    assert os.path.isfile(tree_path)  # check that tree file exists

    # prune the tree
    P_query = ",".join(elems)
    P_cmd = "tree_doctor -n -P {0} {1}".format(P_query, tree_path)
    P_cmd = P_cmd + " -a" if not not_anc else P_cmd
    try:
        tree_info = subprocess.check_output(P_cmd, shell=True).decode("utf-8")
    except subprocess.CalledProcessError:
        tree_info = None
        die("Error! Command {0} failed.".format(P_cmd))
    # apply formatting
    tree_nodes = get_branches(tree_info)
    if und:  # for now it is only one way to format
        tree_info = tree_info.replace("-", "_")
    # write to file or stdout
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write(tree_info)
    f.close if output != "stdout" else None
    return tree_nodes


def fill(data):
    """Replace spaces with NNN's."""
    filled = {}
    for name, seq in data.items():
        filled_seq = seq.replace(" ", "N")
        filled[name] = filled_seq
    return filled


def rm_ref_cols(data, ref_name):
    """Remove cols where ref has gaps."""
    ref_seq = data.get(ref_name)
    updated = {}
    die("Error! Seq {} not found!".format(ref_name)) if not ref_seq else None
    ref_ok = [n for n, c in enumerate(ref_seq) if c != "-"]
    for name, seq in data.items():
        new_seq = "".join([seq[i] for i in ref_ok])
        updated[name] = new_seq
    return updated


def change_case(data, up, lo):
    """Make up/lower case."""
    changed = {}
    for name, seq in data.items():
        new_line = seq.upper() if up else seq.lower()
        changed[name] = new_line
    return changed


def translate(data, force=False):
    """Translate NT to AA sequences."""
    translated = {}  # accumulate the result

    for name, seq in data.items():
        # len must be % 3 == 0!
        if not len(seq) % 3 == 0 and not force:
            die("Error! Codon alignment is required for translation! {} is out of frame".format(name))
        codons = parts(seq, n=3)  # split to codons
        # get corresponding AA for each
        aa_seq = ""  # accumulate here
        for codon in codons:
            AA = nta.get(codon.upper())
            if AA:  # continue if it is OK here
                aa_seq += AA
                continue
            # kill if gaps in codon:
            if "-" in codon and not force:  # something like AT- | must be ATG or ---
                sys.stderr.write("Sequence {} contains frameshifts!\n".format(name))
                die("Error! Codon alignment is required!")
            if "-" in codon and force:
                aa_seq += "X"
            elif "N" in codon:  # if ATN for example - don't know whatta AA
                aa_seq += "X"
            elif "!" in codon:
                aa_seq += "X"

        # save new seq
        translated[name] = aa_seq
    return translated


def rearrange_fasta(data, scale):
    """Apply scale to set of fasta lines."""
    # check if we should change anything
    if scale == 0:  # 0 means "return one line fasta"
        return data  # if so, stop here, we already done it
    rearranged = {}  # save rearranged data here
    for head, sequence in data.items():
        # apply new scale and add to new dictionary
        new_seq = "\n".join(parts(sequence, n=scale))
        rearranged[head] = new_seq
    return rearranged


def trim_data(data, t_from, t_to):
    """Trim the sequences."""
    # check is the limits are violated
    trimmed = {}  # save the result here
    seq_len = len(list(data.values())[0])
    t_to = seq_len if t_to == 0 else t_to  # 0 is default
    if t_from >= seq_len or t_to > seq_len:  # otherwise it is index error
        die("Error! Trim borders are outside the sequence length! {0} letters".format(seq_len))
    for sp, seq in data.items():
        new_seq = seq[t_from: t_to]
        trimmed[sp] = new_seq
        # print(new_seq)
    return trimmed


def misalign(data):
    """Just misalign."""
    new_data = {k: v.replace("-", "") for k, v in data.items()}
    return new_data


def save_fasta(data, order, output):
    """Save one_line fasta in a file."""
    if output == "0":  # skip saving
        return
    f = open(output, "w") if output != "stdout" else sys.stdout  # open the output file
    for head in order:
        seq = data.get(head)
        f.write(">{0}\n".format(head))
        f.write("{0}\n".format(seq))
    f.close() if output != "stdout" else None


def invert_complement(data):
    """Invert complement sequences."""
    inverted = {}
    for name, seq in data.items():
        seq_inv = seq[::-1]
        seq_compl = "".join([compl.get(c) if compl.get(c) else c for c in seq_inv])
        inverted[name] = seq_compl
    return inverted


def unwrap(data, order, lw=20, sw=100):
    """Convinient look at sequences."""
    seq_to_pieces = {}
    pieces_num = 0
    for label, seq in data.items():
        seq_parts = parts(seq, sw)
        seq_to_pieces[label] = seq_parts
        pieces_num = len(seq_parts)
    seq_lens = [len(v) for v in data.values()]
    assert all(x == seq_lens[0] for x in seq_lens)
    for i in range(pieces_num):
        for label in order:
            piece = seq_to_pieces[label][i]
            print("{0: <{1}}: {2}".format(label, lw, piece))
        print()
    

def main(args):
    """Entry point."""
    # test if output files are reachable | if needed
    test_reachable(args.output) if args.input != args.output and args.output != "0" else None
    test_reachable(args.tree) if args.tree else None

    # read initial fasta and check the format
    if not args.phylip:
        data, order = read_fasta(args.input, args.vars, args.rm)  # interrupt if -v
    else:  # it is phylip
        data, order = read_phylip(args.input, args.vars, args.rm)

    if args.up or args.lo:  # up/lo case required
        data = change_case(data, args.up, args.lo)
    if args.copy and args.paste:
        data, order = copy_paste(data, order, args.copy, args.paste)
    elif args.copy or args.paste:
        sys.stderr.write("Warning! Use --copy and --paste together, not separatly.\n")
    if args.unwrap:
        unwrap(data, order)
        sys.exit(0)

    # sort if required
    order = list(sorted(order)) if args.sort else order
    data = data if not args.vs_ref else rm_ref_cols(data, args.vs_ref)
    # fill if requered
    data = fill(data) if args.fill else data

    # translate if required
    data = translate(data, args.force) if args.trans else data

    # build tree if needed
    if args.tree:
        tree_nodes = build_tree(order, args.tree, args.tree_und, args.tree_no_anc)
        not_found = set(order).difference(set(tree_nodes))
        sys.stderr.write("Warning! Not found tree nodes for:\n{0}\n".
                         format(",".join(not_found))) if len(not_found) > 0 else None

    # apply trimming if requered
    data = trim_data(data, args.trim_from, args.trim_to) if args.trim_to > 0 or args.trim_from > 0 else data

    # misalign if required
    if args.misalign:
        data = misalign(data)

    if args.inv:
        data = invert_complement(data)

    # apply scale required
    rearranged = rearrange_fasta(data, args.fasta_scale)

    # and save the output
    save_fasta(rearranged, order, args.output)
    sys.exit(0)


if __name__ == "__main__":
    app = argparse.ArgumentParser()
    app.add_argument("input", type=str)
    app.add_argument("output", type=str, help="use 0 to replace with /dev/null")
    app.add_argument("-n", "--fasta_scale", type=int, default=0, help="number of bases in a fasta line")
    app.add_argument("-v", "--vars", action="store_true", dest="vars", help="Show a list of sequence names")
    app.add_argument("--trim_from", "--tf", type=int, default=0, help="trim all the sequences from")
    app.add_argument("--trim_to", "--tt", type=int, default=0, help="trim all the sequences from")
    app.add_argument("--tree", "-t", type=str, default=None, help="save tree to")
    app.add_argument("--rm", "-r", type=str, default="", help="comma-separated list of sequences to remove")
    app.add_argument("--up", action="store_true", dest="up", help="Make all sequences uppercase")
    app.add_argument("--lo", action="store_true", dest="lo", help="Make all sequences lowercase")
    app.add_argument("--phylip", action="store_true", dest="phylip", help="Input is a phylip file. Convert it into a normal format (fasta)")
    app.add_argument("--tree_und", "-u", action="store_true", dest="tree_und", help="Replace - with _ for tree")
    app.add_argument("--tree_no_anc", "--tnc", action="store_true", dest="tree_no_anc")
    app.add_argument("--sort", "-s", action="store_true", dest="sort", help="Sort sequences in alphabetic order")
    app.add_argument("--trans", action="store_true", dest="trans", help="Translate to AA sequence")
    app.add_argument("--copy", type=str, default=None, help="Copy sequence from... Use with --paste please.")
    app.add_argument("--paste", type=str, default=None, help="Paste copied sequence as...")
    app.add_argument("--fill", action="store_true", dest="fill", help="Replace speces with N's.")
    app.add_argument("--force", "-f", action="store_true", dest="force", help="Ignore errors.")
    app.add_argument("--vs_ref", type=str, default=None, help="Remove columns where ref has gaps.")
    app.add_argument("--misalign", "--mn", action="store_true", dest="misalign", help="Return not aligned fasta.")
    app.add_argument("--inv", action="store_true", dest="inv", help="Invert complement.")
    app.add_argument("--unwrap", "--un", action="store_true", dest="unwrap", help="Unwrap sequences")

    if len(sys.argv) < 3:  # close if very few argumants
        app.print_help()
        sys.exit(0)

    args = app.parse_args()
    # check if arguments are incorrect
    assert args.fasta_scale >= 0  # num of lines must be 0 or more
    assert args.trim_from >= 0 and args.trim_from >= 0  # we don't use negative coordinates
    assert args.up is not True or args.lo is not True  # only one of these might be true)
    main(args)
