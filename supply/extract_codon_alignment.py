#!/usr/bin/env python3
"""Extract codon alignment from TOGA results."""
import argparse
import sys
import os
from collections import defaultdict
from collections import Counter
import random
import string
import tempfile
import subprocess
import getpass
from subprocess import PIPE
from twobitreader import TwoBitFile

__author__ = "Bogdan Kirilenko, 2021"
__version__ = "0.2.6"

# revision history
# v2.5
# added sequences limit, if num of sequences of seqs to realign
# is > 1500 (default) then do not call MACSE
# v2.4
# fixed a bug: script didn't reverse exon sizes list on minus strand
# v2.3
# fixed a bug: script added exons that toga classified as deleted/missing
# actually, it was nearly random CESAR alignment for this exon
# v2.2 fixed a bug: wrong handling of completely deleted exons


SEQ_NUMBER_LIMIT = 1500
CESAR_RESULTS_FILE = "cesar_results.txt"
EXON_SEQ_CLASSES = {"query_exon", "reference_exon"}
DEL_MIS_EX = {"Missing exon", "Deleted exon"}
CODON_ALI_FILE = "codon.fasta"
ORTH_CLASS_FILE = "orthology_classification.tsv"
TAB = "\t"
ONE_TO_ONE = "one2one"
MANY_TO_ONE = "many2one"
SINGLE_COPY_CLASSES = {ONE_TO_ONE, MANY_TO_ONE}
ONE_TO_ZERO = "one2zero"
MACSE_TWO = "macse2"
PRANK = "prank"
REFERENCE_EXON = "reference_exon"
INACT_MUT_DATA = "inact_mut_data.txt"
REFERENCE = "REFERENCE"
GAP = "-"
WHOAMI = "whoami"
DEVSHM = "/dev/shm"
TMP = "/tmp"
STDIN = "/dev/stdin"
UTF_8 = "utf-8"
CODON_GAP = "---"
CODON_XXX = "XXX"
TEMP_TOGA = "temp"
MIN_ARG_NUM = 3

COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "G",
    "n": "n",
}


def generate_random_string(string_len):
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=string_len))


def print_stderr(msg, end="\n"):
    """Print to stderr."""
    sys.stderr.write(f"{msg}{end}")


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def invert_complement(seq):
    """Make inverted-complement sequence."""
    reverse = seq[::-1]
    reverse_complement = "".join(
        [COMPLEMENT.get(c) if COMPLEMENT.get(c) else "N" for c in reverse]
    )
    return reverse_complement


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument(
        "input_dirs", help="File containing a list of TOGA output directories"
    )
    app.add_argument("reference_bed", help="Bed12-file containing reference annotations")
    app.add_argument("transcript_id",
                     help="ID of the aligned transcript (must be present in"
                          " the reference bed file)")

    app.add_argument("--output", "-o", default=None, help="Output file, default stdout")
    app.add_argument(
        "--use_raw_sequences", "--raw", action="store_true", dest="use_raw_sequences"
    )
    app.add_argument(
        "--save_not_aligned",
        "--sna",
        default=None,
        help="Save input sequences in fasta format, works for entire gene alignments only",
    )
    app.add_argument(
        "--allow_one2zero",
        "-z",
        dest="allow_one2zero",
        action="store_true",
        help=(
            "Process orthologous projections in case a species has no "
            "intact/PI/UL projections (class is one2zero)"
        ),
    )
    app.add_argument(
        "--skip_dups",
        "-s",
        dest="skip_dups",
        action="store_true",
        help="Skip not one-2-one orthologs",
    )
    app.add_argument(
        "--align_entirely",
        "-a",
        dest="align_entirely",
        action="store_true",
        help="Do not align exons separately; align the gene entirely",
    )
    app.add_argument(
        "--seq_number_limit",
        default=SEQ_NUMBER_LIMIT,
        help="Exit if number of sequences exceeds the threshold",
    )
    app.add_argument(
        "--temp_dir",
        default=None,
        help="Temp dir, default /dev/shm/username or /tmp/username",
    )
    app.add_argument(
        "--macse_caller",
        default=None,
        help="Macse 2 caller command. Example: java -jar /path/to/macse2.jar "
        "(not just a path to macse2.jar!)",
    )
    app.add_argument(
        "--use_prank",
        action="store_true",
        dest="use_prank",
        help="Use prank instead of MACSE",
    )
    app.add_argument(
        "--prank_executable", "--prank", default="prank", help="Prank executable"
    )
    app.add_argument("--prank_tree", default=None, help="Tree to be used in PRANK")
    app.add_argument(
        "--debug",
        "-d",
        action="store_true",
        dest="debug",
        help="Write debugging information",
    )
    app.add_argument("--reference_2bit", default=None, help="Reference 2bit file")
    app.add_argument(
        "--intermediate_data",
        default=None,
        help="For debugging: directory name to save intermediate data",
    )
    app.add_argument(
        "--min_percent_of_sp_with_one_orth",
        "--mpo",
        default=0.0,
        type=float,
        help="Minimal fraction of species with at least one ortholog "
        " that have exactly one ortholog (one2one or one2many), default 0.0, max 1.0",
    )
    app.add_argument(
        "--max_copies", default=1, type=int, help="Maximal number of copies allowed"
    )
    # Ariadna -> need to check whether it's a good idea
    app.add_argument(
        "--exclude_UL",
        dest="exclude_UL",
        action="store_true",
        help="Do not consider UL projections as paralogous (NOT IMPLEMENTED YET)",
    )

    # if no args: print help message
    if len(sys.argv) < MIN_ARG_NUM:
        app.print_help()
        sys.exit(0)

    args = app.parse_args()

    # arguments sanity checsks
    # if args.use_prank is True and args.prank_tree is None:
    #     print_stderr("Error! The prank mode requires a tree")
    #     print_stderr("Please specify the prank_tree with --prank_tree argument")
    #     sys.exit(0)

    if (
        args.min_percent_of_sp_with_one_orth < 0.0
        or args.min_percent_of_sp_with_one_orth > 1.0
    ):
        print_stderr("Error!")
        print_stderr("--min_percent_of_sp_with_one_orth must be in range 0.0 to 1.0")
        sys.exit(0)

    return args


def read_input_dirs(in_file):
    """Read input directories."""
    with open(in_file, "r") as f:
        paths = [
            os.path.abspath(x.rstrip())
            for x in f
            if x != "\n" and not x.startswith("#")
        ]
    # check whether all files are here
    req_filenames = (
        CODON_ALI_FILE,
        ORTH_CLASS_FILE,
    )
    for path in paths:
        # TODO: fix this, works bad
        if not os.path.isdir(path):
            print_stderr(f"Error! Directory {path} not found.\nAbort.")
            sys.exit(1)
        req_files = (os.path.join(path, r) for r in req_filenames)
        if all(os.path.isfile(r) for r in req_files):
            # all files found: ok
            continue

        print_stderr(f"Error! Directory {path} is incomplete")
        print_stderr("All of the following files must be presented:")
        for f in req_filenames:
            print_stderr(f)
        print_stderr("Abort")
        sys.exit(1)
    return paths


def get_orth_projections(sp_to_dir, transcript_id):
    """Get orthologous projections IDs."""
    sp_to_ids = defaultdict(list)
    for sp, directory in sp_to_dir.items():
        orth_file = os.path.join(directory, ORTH_CLASS_FILE)
        f = open(orth_file, "r")
        f.__next__()
        for line in f:
            line_data = line.rstrip().split(TAB)
            line_trans_id = line_data[1]
            if line_trans_id != transcript_id:
                continue
            # line containing our transcript was caught
            projection_id = line_data[3]
            oclass = line_data[4]
            if oclass == ONE_TO_ZERO:
                # one2zero: ortholog not found
                continue
            is_one_to_one = oclass in SINGLE_COPY_CLASSES
            item = (projection_id, is_one_to_one)
            sp_to_ids[sp].append(item)
        f.close()
    return sp_to_ids


def get_sp_without_orthologs(
    sp_id_to_dir, sp_id_to_orth_projections, transcript_id, z_arg
):
    """Get species without orthologs."""
    sp_without_orthologs = {
        k for k in sp_id_to_dir.keys() if k not in sp_id_to_orth_projections
    }
    if len(sp_without_orthologs) > 0:
        print_stderr(
            f"# Warning! TOGA didn't find {transcript_id} orthologs for the following species:"
        )
        for elem in sp_without_orthologs:
            print_stderr(elem)
    if len(sp_id_to_orth_projections.keys()) == 0 and z_arg is False:
        print_stderr(f"# Error! No species with {transcript_id} orthologs detected")
        print_stderr(
            f"Please add -z argument to also consider Lost/Missing/etc transcripts."
        )
        sys.exit(0)  # maybe 1?
    return sp_without_orthologs


def get_sp_with_single_copy(sp_to_ids):
    """Get sp with single copy."""
    ret = []
    for k, v in sp_to_ids.items():
        num_orths = len(v)
        if num_orths == 1:
            ret.append(num_orths)
    return ret


def filter_in_data(sp_id_to_orth_projections, sp_id_to_dir, skip_dups, t_, max_copies):
    """Prepare dict: sp_id: (dir, seq_ids)."""
    ret = {}
    for k, v in sp_id_to_orth_projections.items():
        seq_ids = [x[0] for x in v]
        is_one2one = v[0][1]
        seq_num = len(seq_ids)
        print_stderr(f"# Found {seq_num} {t_} orthologs for {k}")

        if skip_dups and is_one2one is False:
            print_stderr(f"# Omitting query {k}: not one2one")
            continue

        if seq_num > max_copies:  # for example > 4 copies
            print_stderr(f"# Omitting query {k}: has {seq_num} > {max_copies} copies")
            continue

        directiory = sp_id_to_dir[k]
        val = (seq_ids, directiory)
        ret[k] = val
    return ret


def get_exon_data(ref_bed, transcript_id, two_bit):
    """Extract exon sizes.

    We will need it to split codon alignments into separate exons.
    """
    coding_sequence_exons = []
    transcript_line = None
    f = open(ref_bed, "r")
    for line in f:
        line_data = line.rstrip().split(TAB)
        if line_data[3] == transcript_id:
            transcript_line = line_data
            break
    f.close()
    # sanity check
    if transcript_line is None:
        print_stderr(f"# Error! Cannot find {transcript_id} in the {ref_bed}, abort")
        sys.exit(1)
    # trim UTRs
    chrom = transcript_line[0]
    chrom_start = int(transcript_line[1])
    thick_start = int(transcript_line[6])  # CDS start and end
    thick_end = int(transcript_line[7])
    block_count = int(transcript_line[9])
    strand = transcript_line[5]
    block_sizes = [int(x) for x in transcript_line[10].split(",") if x != ""]
    block_starts = [int(x) for x in transcript_line[11].split(",") if x != ""]
    block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
    block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
    block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
    block_new_starts, block_new_ends = [], []
    ref_genome_reader = TwoBitFile(two_bit)
    chrom_seq = ref_genome_reader[chrom]

    for block_num in range(block_count):
        # go block-by-block
        block_start = block_abs_starts[block_num]
        block_end = block_abs_ends[block_num]

        # skip the block if it is entirely UTR
        if block_end <= thick_start:
            continue
        elif block_start >= thick_end:
            continue

        block_new_start = block_start if block_start >= thick_start else thick_start
        block_new_end = block_end if block_end <= thick_end else thick_end
        block_new_starts.append(block_new_start - thick_start)
        block_new_ends.append(block_new_end - thick_start)
        exon_seq = chrom_seq[block_new_start:block_new_end]
        exon_seq = exon_seq if strand == "+" else invert_complement(exon_seq)
        coding_sequence_exons.append(exon_seq)

    cds_exon_sizes = [
        block_new_ends[i] - block_new_starts[i] for i in range(len(block_new_starts))
    ]
    cds_exon_sizes_s = cds_exon_sizes if strand == "+" else cds_exon_sizes[::-1]
    coding_sequence_exons_s = (
        coding_sequence_exons if strand == "+" else coding_sequence_exons[::-1]
    )
    coding_sequence = "".join(coding_sequence_exons_s)
    codon_numbers = parts(coding_sequence, 3)
    return cds_exon_sizes_s, codon_numbers


def extract_seq_generator(fasta_file):
    """Read very big fasta, sequence by sequence."""
    with open(fasta_file) as f:
        accum = []
        for line in f:
            if line.startswith(">"):
                if len(accum) == 0:
                    accum = [
                        line,
                    ]
                    continue
                yield accum
                accum = [
                    line,
                ]
            else:
                if line.startswith(" "):
                    continue
                elif line == "\n":
                    continue
                elif line.startswith("#"):
                    continue
                accum.append(line)
        f.close()
    yield accum


def split_proj_name(proj_name):
    """Split projection name.

    Projections named as follows: ${transcript_ID}.{$chain_id}.
    This function splits projection back into transcript and chain ids.
    We cannot just use split("."), because there might be dots
    in the original transcript ID.
    """
    proj_name_split = proj_name.split(".")
    q_num_str = proj_name_split[-1]
    trans_name = ".".join(proj_name_split[:-1])
    return trans_name, q_num_str


def extract_seq_from_fasta_by_id(fasta_file, seq_ids):
    ret = []
    f = open(fasta_file, "r")
    found_ids = set()

    prev = "\n"
    while prev != "\n":
        prev = f.__next__()

    for lnum, line in enumerate(f, 1):
        curr = line

        if curr == "\n" or prev == "\n":
            # no need this
            prev = curr
            continue

        if curr.startswith(">") and not prev.startswith(">"):
            # inverted pair, go to next
            prev = curr
            continue

        if curr.startswith(">") and prev.startswith(">"):
            # error! corrupted lines
            prev_header = prev.lstrip(">").rstrip().split(" | ")
            prev_id = prev_header[0]
            print_stderr(f"Error! Line {lnum} corrupted! File {fasta_file}")
            print_stderr(f"No seq for following ID: {prev_id}")
            print_stderr(f"### Prev:\n{prev}")
            print_stderr(f"### Curr:\n{curr}")

            prev = curr
            continue

        if prev.startswith(">") and not curr.startswith(">"):
            # proper pair of sequences
            header = prev.lstrip(">").rstrip().split(" | ")
            seq_id = header[0]
            if seq_id not in seq_ids:
                # not the seq we need
                prev = curr
                continue
            # seq we need!
            # control
            if header[1] != "CODON":
                print_stderr(f"Error! Line {lnum}: non codon sequence reported!")
                prev = curr
                continue

            found_ids.add(seq_id)
            ref_que_field = header[2]
            seq = curr

            item = (seq_id, ref_que_field, seq)
            ret.append(item)
            prev = curr
            continue

        else:
            # what else is possible?
            print_stderr(f"Error! Line {lnum} File {fasta_file}")
            print_stderr(f"Prev line:\n{prev}")
            print_stderr(f"Curr line:\n{curr}")
            prev = curr
            continue
    f.close()

    not_found = seq_ids.difference(found_ids)
    if len(not_found) > 0:
        print_stderr(f"Error! Projections {not_found} not found in {fasta_file}")
    return ret


def extract_codon_sequences(sp_dir_seq_id_data, debug=False):
    """Extract codon sequences."""
    id_to_seq = {}
    proj_to_num_seq = defaultdict(set)  # each proj must appear twice; for sanity checks
    for sp, (seq_ids, directory) in sp_dir_seq_id_data.items():
        proj_id_set = set(seq_ids)
        codon_ali_file = os.path.join(directory, CODON_ALI_FILE)
        # seq_gen = extract_seq_generator(codon_ali_file)
        seq_gen = extract_seq_from_fasta_by_id(codon_ali_file, proj_id_set)

        for fasta_elem in seq_gen:
            proj_id = fasta_elem[0]
            is_ref = fasta_elem[1] == REFERENCE
            seq = fasta_elem[2].rstrip()

            seq_id = (sp, proj_id, is_ref)
            id_to_seq[seq_id] = seq
            proj_to_num_seq[(sp, proj_id)].add(is_ref)
    # sanity checks
    for seq_id, _set in proj_to_num_seq.items():
        sp_ = seq_id[0]
        proj_id_ = seq_id[1]
        # we gonna extract ref and query for all sequences
        if _set == {True, False}:
            continue
        if len(_set) == 0:
            # a crazy situation, should never happen
            raise IndexError(f"Error! Corrupted {seq_id}")
        print_stderr(f"Error! Species {sp_} projection {proj_id_}:")
        print_stderr("Expected to get 2 sequences (reference and query)")
        ref_found = list(_set)[0]  # extract the only element of list
        if ref_found:
            print_stderr("Got only the reference sequence, abort")
        else:
            print_stderr("Got only the query sequence, abort")
        sys.exit(1)
    return id_to_seq


def convert_seq_to_fasta(seq_data):
    """Convert sequences dict to fasta-formatted string."""
    reference_seq = [v for k, v in seq_data.items() if k[2] is True][0]
    query_seqs = {k: v for k, v in seq_data.items() if k[2] is False}
    ret = [
        ">REFERENCE\n",
    ]
    reference_seq_fmt = reference_seq.replace(GAP, "").replace(" ", "")
    ret.append(f"{reference_seq_fmt}\n")
    # fasta_headers = []
    for k, v in query_seqs.items():
        sp_name = k[0]
        proj_id = k[1]
        fasta_header = f"{sp_name}\t{proj_id}"
        # fasta_headers.append(fasta_header)
        seq_id = f">{fasta_header}\n"
        seq_fmt = v.replace(GAP, "").replace(" ", "")
        ret.append(seq_id)
        # if len(seq_fmt) == 0:
        #    raise ValueError("Error! Empty fasta sequence")
        ret.append(f"{seq_fmt}\n")
    fasta_line = "".join(ret)
    return fasta_line


def macse_alignment(in_fasta, temp_dir, macse_caller, v=False):
    """Align fasta with MACSE."""
    # MACSE is unable to write to stdout
    # also it forces user to create AA alignment
    out_file, to_del_file = "", ""
    try:
        print_stderr("# Calling MACSE2...")
        fd1, out_file = tempfile.mkstemp(suffix=".fa", prefix="macse_out", dir=temp_dir)
        fd2, to_del_file = tempfile.mkstemp(
            suffix=".fa", prefix="to_del_", dir=temp_dir
        )
        cmd = f"{macse_caller} -prog alignSequences -seq {STDIN} -out_NT {out_file} -out_AA {to_del_file}"
        if v:
            print_stderr("Aligning:")
            print_stderr(in_fasta)
        p = subprocess.Popen(cmd, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
        _, stderr_ = p.communicate(input=in_fasta.encode())
        rc = p.returncode
        if rc != 0:  # MACSE crashed
            print_stderr("# Error! Macse CRASHED")
            err_msg = stderr_.decode(UTF_8)
            print_stderr(err_msg)
            sys.exit(1)
        with open(out_file, "r") as f:
            aligned_fasta = f.read()
    finally:  # to delete these files anyway
        os.remove(out_file) if os.path.isfile(out_file) else None
        os.remove(to_del_file) if os.path.isfile(to_del_file) else None
        os.close(fd1)
        os.close(fd2)
    return aligned_fasta


def __fasta_replace_X_with_N(in_fasta):
    out_lines = [x.replace("X", "N") if not x.startswith(">") else x for x in in_fasta.split("\n")]
    return "\n".join(out_lines)


def prank_alignment(in_fasta, in_tree, temp_dir, prank_caller, v=False):
    """Errors that make prank crash:
    
    Only one sequence in tempdir/62JEFU1HL8LY.input.fasta!
    """
    prank_out_basename = generate_random_string(12)
    input_path = os.path.join(temp_dir, f"{prank_out_basename}.input.fasta")
    with open(input_path, "w") as f:
        fixed_fasta = __fasta_replace_X_with_N(in_fasta)
        f.write(fixed_fasta)
    out_prank_arg = os.path.join(temp_dir, prank_out_basename)
    cmd = f"{prank_caller} -d={input_path} -once -F -DNA -codon -o={out_prank_arg}"
    if in_tree:  # if tree is provided: do this
        cmd += f" -t={in_tree} -prunetree"

    p = subprocess.Popen(cmd, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
    _, stderr_ = p.communicate()
    rc = p.returncode
    if rc != 0:  # PRANK crashed
        print_stderr("# Error! Prank CRASHED")
        err_msg = stderr_.decode(UTF_8)
        print_stderr(err_msg)
        # os.remove(input_path) if os.path.isfile(input_path) else None
        sys.exit(1)

    related_files = [
        os.path.join(temp_dir, x)
        for x in os.listdir(temp_dir)
        if x.startswith(prank_out_basename)
    ]
    best_fas_arr = [x for x in related_files if x.endswith("best.fas")]
    print(related_files)
    print(best_fas_arr)
    if len(best_fas_arr) != 1:
        print_stderr("Prank error, > 1 best fas file")
        sys.exit(1)

    best_fas = best_fas_arr[0]
    with open(best_fas, "r") as f:
        aligned_fasta = f.read()

    for path in related_files:
        os.remove(path) if os.path.isfile(path) else None
    
    # TODO: return original sequence names
    return aligned_fasta


def get_temp_dir(temp_dir_arg):
    """Create temporary directory and return path."""
    if temp_dir_arg is None:
        # default temp directory name
        # use username to distinguish temp files
        # whoami = subprocess.check_output(WHOAMI, shell=True).decode(UTF_8)[:-1]
        whoami = getpass.getuser()
        # /dev/shm is better but may be not available
        dev_shm_available = os.path.isdir(DEVSHM)
        if dev_shm_available:
            temp_dir = os.path.join(DEVSHM, whoami)
        else:
            temp_dir = os.path.join(TMP, whoami)
    else:
        # not default: use what user asks for
        temp_dir = temp_dir_arg
    # create this dir in case it doesn't exist:
    os.mkdir(temp_dir) if not os.path.isdir(temp_dir) else None
    return temp_dir


def get_codon_to_exon_mapping(exon_sizes):
    """Index: codon number, value: exon number."""
    ret = []
    remaining_nucl = 0
    for num, exon_size in enumerate(exon_sizes, 1):
        exon_size_no_left = exon_size - remaining_nucl
        full_codons = exon_size_no_left // 3
        for _ in range(full_codons):
            exon_id = (num, 0)  # 0 -> full codon
            ret.append(exon_id)
        rem_ = exon_size_no_left % 3
        remaining_nucl = 3 - rem_ if rem_ > 0 else 0

        if remaining_nucl != 0:
            # there is a split codon
            exon_id = (num, 1)  # 1 -> split codon
            ret.append(exon_id)
    return ret


def fix_codon_num(curr_codon_num, exp_codon, codon_sequences):
    """Fix codon number: find a codon that match."""
    temp = curr_codon_num - 1
    _max_ind = len(codon_sequences) - 1
    temp = temp if temp <= _max_ind else _max_ind
    while temp >= 0:
        # print(exp_codon, temp, codon_sequences[temp])
        if codon_sequences[temp] == exp_codon:
            return temp
        temp -= 1
    raise IndexError("Cannot find a proper codon number")


def split_into_exons(sp_to_codon_alis, codon_to_exon, codon_to_seq, debug=False):
    """Split sequences into exons."""
    ret = {}
    ali_id_to_seqs = defaultdict(dict)
    codon_tot_num_ = len(codon_to_exon)
    if debug:
        print_stderr("Codon to exon_dict")
        print_stderr(codon_to_exon)
        print_stderr(f"In total {codon_tot_num_} codons")

        print_stderr(f"Codon sequences:\n{codon_to_seq}")
    for seq_id, seq in sp_to_codon_alis.items():
        sp, proj_id, is_ref = seq_id
        key_ = (sp, proj_id)
        ali_id_to_seqs[key_][is_ref] = seq

    for k, v in ali_id_to_seqs.items():
        print_stderr(f"processing {k}") if debug else None
        ref_codons = v[True].split()
        que_codons = v[False].split()
        codon_num = len(ref_codons)
        assert codon_num == len(que_codons)  # if failed, write to Bogdan PLS
        ref_codon_num = 0
        exon_to_sequences = defaultdict(list)
        exp_codon_lock = False  # a workaround
        codon_ali_len = len(ref_codons)

        for i_num, (r_c, q_c) in enumerate(zip(ref_codons, que_codons)):
            exon_num = codon_to_exon[ref_codon_num]
            debug_line = f"{k}: extracting codon {ref_codon_num} out of {codon_tot_num_} | exon {exon_num}"
            print_stderr(debug_line) if debug else None

            expected_ref_codon = codon_to_seq[ref_codon_num]
            debug_line = (
                f"{r_c} -> {expected_ref_codon} {q_c} {ref_codon_num} {exon_num}"
            )
            print_stderr(debug_line) if debug else None
            ref_is_gap = r_c == CODON_GAP
            ref_is_XXX = r_c == CODON_XXX
            ref_has_N = "N" in r_c
            ref_is_exp = r_c.lower() == expected_ref_codon.lower()
            ref_is_undefined = ref_is_gap or ref_is_XXX or ref_has_N

            debug_line = f"# ref gap: {ref_is_gap} ref unknown: {ref_is_XXX} ref expected: {ref_is_exp}"
            print_stderr(debug_line) if debug else None

            if not ref_is_gap:
                ref_codon_num += 1
            if (
                ref_is_exp is False
                and ref_is_undefined is False
                and exp_codon_lock is False
            ):
                print_stderr(f"Warning! {k}: Fixing ref codon_num..")
                wrong_codon_num = ref_codon_num
                ref_codon_num = fix_codon_num(
                    wrong_codon_num, expected_ref_codon, codon_to_seq
                )
                print_stderr(f"Corrected from {wrong_codon_num} to {ref_codon_num}")

            codon_pair = (r_c, q_c)
            exon_to_sequences[exon_num].append(codon_pair)
            is_last = i_num == codon_ali_len - 1
            if ref_codon_num == codon_tot_num_ and is_last is False:
                # this will cause an index error
                # need to handle this overflow correctrly
                exp_codon_lock = True
                print_stderr(f"Warning! {k}: Codon sequence overflow...")
                ref_codon_num = codon_tot_num_ - 1

        ret[k] = exon_to_sequences
    return ret


def reformat_codon_alis(sp_to_exon_codon_alis, all_exons):
    """Restore initial sequences format."""
    ret = defaultdict(dict)
    for (sp, proj_id), v in sp_to_exon_codon_alis.items():
        for exon in all_exons:
            exon_related_seq = v[exon]
            ref_codons = [x[0] for x in exon_related_seq]
            que_codons = [x[1] for x in exon_related_seq]
            ref_seq = "".join(ref_codons)
            que_seq = "".join(que_codons)
            ref_seq_key = (sp, proj_id, True)
            que_seq_key = (sp, proj_id, False)
            ret[exon][ref_seq_key] = ref_seq
            ret[exon][que_seq_key] = que_seq
    return ret


def merge_exon_fastas(fastas_list, fasta_headers):
    """Merge per-exon fastas."""
    id_to_chunks = defaultdict(list)
    fasta_lines = []
    for fasta in fastas_list:
        # ordered exon-by-exon
        fasta_headers_added = []
        fasta_chunks = fasta.split(">")[1:]  # remove element 0: void
        fasta_seq_lens = []  # for checks, must be the same
        for chunk in fasta_chunks:
            lines = chunk.split("\n")
            header = lines[0]
            fasta_headers_added.append(header)
            seq = "".join(lines[1:])
            fasta_seq_lens.append(len(seq))
            id_to_chunks[header].append(seq)
        missing_headers = fasta_headers.difference(fasta_headers_added)
        # checks
        if any(fasta_seq_lens[0] != x for x in fasta_seq_lens):
            print_stderr("Error! Broken alignment (to fill later)")
            sys.exit(1)
        if len(missing_headers) == 0:
            # everything is OK: no deleted exons and absent fasta entries
            continue
        # some fasta headers are missing: need to fill the gaps
        filler_size = fasta_seq_lens[0]
        filler = "-" * filler_size
        for header in missing_headers:
            id_to_chunks[header].append(filler)

    # generate fasta
    for k, v in id_to_chunks.items():
        fasta_lines.append(f">{k}\n")
        seq = "".join(v)
        fasta_lines.append(f"{seq}\n")
    return "".join(fasta_lines)


def get_reference_2bit(ref_2bit_arg, wd):
    """Get reference 2bit file."""
    if ref_2bit_arg:
        # set by user
        return ref_2bit_arg
    ref_2bit = os.path.join(wd, "t2bit.link")
    return ref_2bit


def split_into_codons(ref_seq, que_seq):
    """Split sequences into codons.

    Reference-based.
    """
    ref_codons = [
        [],
    ]
    que_codons = [
        [],
    ]
    curr_num = 0
    codon_num = 0
    for t, q in zip(ref_seq, que_seq):
        if curr_num == 3:
            curr_num = 0
            ref_codons.append([])
            que_codons.append([])
        if t != "-":
            curr_num += 1
        # add to the last codon:
        ref_codons[-1].append(t)
        que_codons[-1].append(q)
    # list of lists to list of strings
    ref_str_codons = ["".join(x) for x in ref_codons]
    que_str_codons = ["".join(x) for x in que_codons]
    ref_codon_seq = " ".join(ref_str_codons)
    que_codon_seq = " ".join(que_str_codons)
    return ref_codon_seq, que_codon_seq


def extract_codon_sequences_from_raw(sp_dir_seq_id_data):
    """Extract codon sequences."""
    id_to_seq = {}
    for sp, (seq_ids, directory) in sp_dir_seq_id_data.items():
        proj_id_set = set(seq_ids)
        codon_ali_file = os.path.join(directory, TEMP_TOGA, CESAR_RESULTS_FILE)
        seq_gen = extract_seq_generator(codon_ali_file)
        projection_to_seq_list = defaultdict(list)
        # extract sequence data from file
        for num, fasta_elem in enumerate(seq_gen):
            header = fasta_elem[0]
            header_data = header.lstrip(">").rstrip().split(" | ")
            seq_class = header_data[-1]
            if seq_class not in EXON_SEQ_CLASSES:
                continue
            is_ref = seq_class == REFERENCE_EXON
            trans_id = header_data[0]
            exon_num = int(header_data[1])
            chain_id = int(header_data[2])
            proj_id = f"{trans_id}.{chain_id}"
            if proj_id not in proj_id_set:
                continue
            seq = "".join(x.rstrip() for x in fasta_elem[1:]).replace("X", "N")
            seq_elem = (exon_num, is_ref, seq)
            projection_to_seq_list[proj_id].append(seq_elem)
        # merge sequences in a proper order
        for proj_id, seq_data in projection_to_seq_list.items():
            reference_sequences_dat = sorted(
                [x for x in seq_data if x[1] is True], key=lambda y: y[0]
            )
            query_sequences_dat = sorted(
                [x for x in seq_data if x[1] is False], key=lambda y: y[0]
            )
            ref_exons, que_exons = [], []
            for ref_exon, que_exon in zip(reference_sequences_dat, query_sequences_dat):
                ref_eseq = ref_exon[2]
                que_eseq = que_exon[2]
                if len(que_eseq) != len(ref_eseq):
                    que_eseq = "-" * len(ref_eseq)
                ref_exons.append(ref_eseq)
                que_exons.append(que_eseq)
            ref_sequence = "".join(ref_exons)
            que_sequence = "".join(que_exons)
            ref_codons, que_codons = split_into_codons(ref_sequence, que_sequence)
            # print(sp, proj_id)
            # print(ref_codons)
            # print(que_codons)
            ref_seq_id = (sp, proj_id, True)
            que_seq_id = (sp, proj_id, False)
            id_to_seq[ref_seq_id] = ref_codons
            id_to_seq[que_seq_id] = que_codons
    return id_to_seq


def get_fasta_headers(sp_keys):
    """Get a list of expected fasta headers.

    This is necessary in case some exon sequences are deleted.
    """
    ret = set([f"{x[0]}\t{x[1]}" for x in sp_keys])
    return ret


def _save_intermed(in_fasta_str, name_raw, intermed_data_dir):
    """Save intermediate fasta if required."""
    if intermed_data_dir is None:
        return
    name = name_raw.replace(" ", "").replace("(", "").replace(")", "")
    path = os.path.join(intermed_data_dir, name)
    with open(path, "w") as f:
        f.write(in_fasta_str)


def extract_del_exons(sp_to_dir_and_seq_id):
    """Extract deleted and missing exons."""
    ret = {}
    for sp, data in sp_to_dir_and_seq_id.items():
        projections = set(data[0])
        sp_dir = data[1]
        # TODO: check that inact mut data is presence earlier!
        inact_mut_file = os.path.join(sp_dir, INACT_MUT_DATA)
        proj_to_del_exons = defaultdict(list)
        f = open(inact_mut_file, "r")
        for line in f:
            line_data = line.rstrip().lstrip("# ").split("\t")
            if len(line_data) != 8:
                continue
            if line_data[4] not in DEL_MIS_EX:
                continue
            trans = line_data[0]
            chain_id = line_data[1]
            proj_id = f"{trans}.{chain_id}"
            if proj_id not in projections:
                continue
            exon_num = int(line_data[2])  # 0-based numbers outside
            proj_to_del_exons[proj_id].append(exon_num)
        f.close()
        ret[sp] = proj_to_del_exons
    return ret


def mask_deleted_and_missing_seq(
    exon_to_alignments_codon_lists_raw, sp_to_deleted_exons
):
    ret = {}
    for key, ex_id_to_seq in exon_to_alignments_codon_lists_raw.items():
        sp_dir = key[0]
        proj_id = key[1]
        if not sp_to_deleted_exons.get(sp_dir):
            # no deleted exons in this species
            # do not change anything
            ret[key] = ex_id_to_seq
            continue
        proj_id_to_del_exons = sp_to_deleted_exons[sp_dir]
        exons_del_in_this_proj = proj_id_to_del_exons.get(proj_id, None)
        if not exons_del_in_this_proj:
            # this projection is not affected: not do anything
            ret[key] = ex_id_to_seq
            continue
        # we got deleted exons here, need to get their ids and mask
        # (exon_num, split exon / not)
        # if exon num in exons del: delete it
        # if exon_num + 1 in del list and split exon: also del: this split exon is also affected
        exons_to_mask = set(
            [
                k
                for k in ex_id_to_seq
                if k[0] in exons_del_in_this_proj
                or (k[0] + 1 in exons_del_in_this_proj and k[1] == 1)
            ]
        )
        ex_id_to_seq_upd = {}
        for exon_id, codon_ali in ex_id_to_seq.items():
            to_mask_it = exon_id in exons_to_mask
            if to_mask_it is False:
                # keep this exon as is
                ex_id_to_seq_upd[exon_id] = codon_ali
                continue
            # need to mask each 2nd codon
            updated_codons = [(c[0], "---") for c in codon_ali]
            ex_id_to_seq_upd[exon_id] = updated_codons
        ret[key] = ex_id_to_seq_upd
    return ret


def extract_proj_ids_from_given_codon_ali_file(codon_ali_file, given_trans_id):
    """For a given codon.fasta return all projections for a given transcript present."""
    f = open(codon_ali_file, "r")
    ret = set()
    for line in f:
        if not line.startswith(">"):
            continue
        line_data = line.lstrip(">").rstrip().split(" | ")
        projection_id = line_data[0]
        transcript_id, _ = split_proj_name(projection_id)
        if transcript_id == given_trans_id:
            ret.add(projection_id)
    f.close()
    return ret


def extract_to_zero_orthologs(sp_id_to_dir, sp_without_orthologs, transcript_id):
    """For species without I/PI/UL orthologs extract at least something.

    For example, lost orthologous genes.
    Ariadna: if we need to adjust "allow one2zero" behaviour (for RELAX for instance),
    this is the function to be modified. Right now it stupidly takes ALL projections
    of a given transcript from ALL species with one2zero. Probably need to filter out
    projections that make no sense (like veeeery short etc).
    """
    ret = {}
    print_stderr("Trying to extract orthologous projections for:")
    for sp in sp_without_orthologs:
        print_stderr(f"{sp}\n")
        sp_dir = sp_id_to_dir[sp]
        codon_ali_path = os.path.join(sp_dir, CODON_ALI_FILE)
        orth_projections = list(
            extract_proj_ids_from_given_codon_ali_file(codon_ali_path, transcript_id)
        )
        if len(orth_projections) == 0:
            continue

        if len(orth_projections) > 1:
            orth_list = [(x, False) for x in orth_projections]
        else:
            orth_list = [(orth_projections[0], True)]
        ret[sp] = orth_list
    return ret


def main():
    """Entry point."""
    args = parse_args()
    temp_dir = get_temp_dir(args.temp_dir)
    print_stderr(f"Using temp dir {temp_dir}") if args.debug else None
    # define aligner binaries
    macse_caller = MACSE_TWO if args.macse_caller is None else args.macse_caller
    prank_caller = PRANK if args.prank_executable is None else args.prank_executable
    print_stderr(f"Macse caller: {macse_caller}") if args.debug else None
    print_stderr(
        f"Prank caller: {prank_caller} (if prank mode)"
    ) if args.debug else None

    # extract and check input directories:
    in_dirs = read_input_dirs(args.input_dirs)
    sp_id_to_dir = {os.path.basename(x): x for x in in_dirs}

    # extract IDs of orthologous sequences for the given gene:
    sp_id_to_orth_projections = get_orth_projections(sp_id_to_dir, args.transcript_id)
    if args.intermediate_data:
        os.mkdir(args.intermediate_data) if not os.path.isdir(
            args.intermediate_data
        ) else None
    # catch potential issues
    sp_without_orthologs = get_sp_without_orthologs(
        sp_id_to_dir, sp_id_to_orth_projections, args.transcript_id, args.allow_one2zero
    )

    # will need this for further computations
    sp_with_at_least_one_ortholog = set(sp_id_to_orth_projections).difference(
        sp_without_orthologs
    )

    # get number of species with single copy
    sp_with_single_copy = get_sp_with_single_copy(sp_id_to_orth_projections)

    # check "minimal % of species with single copy requirements"
    num_atleast_one_copy = len(sp_with_at_least_one_ortholog)
    num_sp_single_copy = len(sp_with_single_copy)

    prop = num_sp_single_copy / num_atleast_one_copy
    if prop < args.min_percent_of_sp_with_one_orth:
        print_stderr(
            f"Minimal % of species with single orthologous copy: {args.min_percent_of_sp_with_one_orth}"
        )
        print_stderr(f"Got: {prop} - less than required")
        print_stderr(f"Sp with at least one copy: {num_atleast_one_copy}")
        print_stderr(f"Sp with single copy: {num_sp_single_copy}")
        print_stderr("Abort")
        sys.exit(0)

    # if -z option and some species have no I/PI/UL orthologs
    # extract what is possible to extract
    if len(sp_without_orthologs) > 0 and args.allow_one2zero:
        extra_sp_to_projections = extract_to_zero_orthologs(
            sp_id_to_dir, sp_without_orthologs, args.transcript_id
        )
        sp_id_to_orth_projections.update(extra_sp_to_projections)

    if args.allow_one2zero and len(sp_id_to_orth_projections) == 0:
        print_stderr(
            f"Unfortunately, none of the listed species have an ortholog of {args.transcript_id}"
        )
        sys.exit(0)

    # prepare dict to work on: exclude not one2ones
    # please note that --skip_dups option also excludes many2ones!
    sp_to_dir_and_seq_id = filter_in_data(
        sp_id_to_orth_projections,
        sp_id_to_dir,
        args.skip_dups,
        args.transcript_id,
        args.max_copies,
    )

    sp_to_deleted_exons = extract_del_exons(sp_to_dir_and_seq_id)
    # final checks; extract pairwise alignments
    if len(sp_to_dir_and_seq_id) == 0:
        # nothing left after filters
        print_stderr(f"{args.transcript_id}: No species left after filters, abort")
        sys.exit(0)

    # USE RAW SEQUENCES: CONTAINING FS, STOPS, ETC
    if args.use_raw_sequences:
        sp_to_codon_alis = extract_codon_sequences_from_raw(sp_to_dir_and_seq_id)
    # USE FILTERED CODON SEQUENCES
    # in fact, not sp to something
    # key is: (species, projection_ID, IS_REF)
    # value: codon alignment sequence from codon.fasta
    else:
        sp_to_codon_alis = extract_codon_sequences(
            sp_to_dir_and_seq_id, debug=args.debug
        )

    # CHECK WHETHER WE ARE WITHIN THE LIMIT
    query_seq_num = len([k for k, v in sp_to_codon_alis.items() if k[2] is False])
    if query_seq_num > args.seq_number_limit:
        # no need to call the program
        f = open(args.output, "w") if args.output else sys.stdout
        f.write(f"#WARNING! Sequences num limit was set to {args.seq_number_limit}\n")
        f.write(f"Extracted {query_seq_num} sequences to realign, quit\n")
        f.close() if args.output else None
        sys.exit(0)

    # ALIGNING ENTIRE GENE
    if args.align_entirely:
        # if so: do not split it into codons, just align with MACSE
        in_fasta = convert_seq_to_fasta(sp_to_codon_alis)
        # save input fasta if required:
        if args.save_not_aligned:
            f = open(args.save_not_aligned, "w")
            f.write(in_fasta)
            f.close()
        # Call MACSE and save the result:
        print_stderr(f"# {args.transcript_id}: aligning the gene entirely")
        if args.use_prank:
            aligned_fasta = prank_alignment(
                in_fasta, args.prank_tree, temp_dir, prank_caller, v=args.debug
            )
        else:
            aligned_fasta = macse_alignment(
                in_fasta, temp_dir, macse_caller, v=args.debug
            )
    # ALIGNING EXON-BY-EXON
    else:
        fasta_headers_all = set([f"{x[0]}\t{x[1]}" for x in sp_to_codon_alis.keys()])
        # we like to align exons separately, first of all get
        # codon number -> exon number mapping
        reference_2bit = get_reference_2bit(args.reference_2bit, in_dirs[0])
        transcript_exon_sizes, codon_to_seq = get_exon_data(
            args.reference_bed, args.transcript_id, reference_2bit
        )
        print(transcript_exon_sizes, codon_to_seq)
        codon_to_exon = get_codon_to_exon_mapping(transcript_exon_sizes)
        exon_to_alignments_codon_lists_raw = split_into_exons(
            sp_to_codon_alis, codon_to_exon, codon_to_seq, args.debug
        )
        # print(exon_to_alignments_codon_lists)
        exon_to_alignments_codon_lists = mask_deleted_and_missing_seq(
            exon_to_alignments_codon_lists_raw, sp_to_deleted_exons
        )
        all_exons = sorted(set(codon_to_exon), key=lambda x: (x[0], x[1]))
        exon_to_alignments = reformat_codon_alis(
            exon_to_alignments_codon_lists, all_exons
        )
        # realign this exon-by-exon
        aligned_fastas = []
        for exon in all_exons:
            print_stderr(f"# {args.transcript_id}: aligning exon {exon}")
            exon_seq_data = exon_to_alignments[exon]
            in_fasta = convert_seq_to_fasta(exon_seq_data)
            _save_intermed(in_fasta, f"exon_{exon}", args.intermediate_data)
            if args.use_prank:
                aligned_fasta_exon = prank_alignment(
                    in_fasta, args.prank_tree, temp_dir, prank_caller, v=args.debug
                )
            else:
                aligned_fasta_exon = macse_alignment(
                    in_fasta, temp_dir, macse_caller, v=args.debug
                )
            _save_intermed(
                aligned_fasta_exon, f"aligned_exon_{exon}", args.intermediate_data
            )
            aligned_fastas.append(aligned_fasta_exon)
        aligned_fasta = merge_exon_fastas(aligned_fastas, fasta_headers_all)

    f = open(args.output, "w") if args.output else sys.stdout
    f.write(aligned_fasta)
    f.close() if args.output else None


if __name__ == "__main__":
    main()
