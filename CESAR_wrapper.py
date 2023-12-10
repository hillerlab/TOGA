#!/usr/bin/env python3
"""CESAR2.0 wrapper.

Retrieve exons from the query genome.
"""
import argparse
import os
import sys
import subprocess
import uuid
import math
from datetime import datetime as dt
from re import finditer, IGNORECASE
from collections import defaultdict
import ctypes
from operator import and_
from functools import reduce
from twobitreader import TwoBitFile
from modules.common import parts, bed_extract_id_text
from modules.common import bed_extract_id, chain_extract_id
from modules.common import make_cds_track
from modules.common import eprint
from modules.common import die
from modules.common import flatten
from modules.inact_mut_check import inact_mut_check
from modules.parse_cesar_output import parse_cesar_out
from constants import Constants
from constants import GENETIC_CODE
from constants import COMPLEMENT_BASE
from version import __version__

__author__ = "Bogdan M. Kirilenko"

LOCATION = os.path.dirname(os.path.abspath(__file__))
VERBOSE = False

TMP_NAME_SIZE = 10

EXP_REG_EXTRA_FLANK = 50
EXON_SEQ_FLANK = 10
ERR_CODE_FRAGM_ERR = 2

# alias; works for Hillerlab-only
two_bit_templ = "/projects/hillerlab/genome/gbdb-HL/{0}/{0}.2bit"
chain_alias_template = "/projects/hillerlab/genome/gbdb-HL/{0}/lastz/vs_{1}/axtChain/{0}.{1}.allfilled.chain.gz"

# connect shared lib; define input and output data types
chain_coords_conv_lib_path = os.path.join(
    LOCATION, "modules", "chain_coords_converter_slib.so"
)
ch_lib = ctypes.CDLL(chain_coords_conv_lib_path)
ch_lib.chain_coords_converter.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_char_p),
]
ch_lib.chain_coords_converter.restype = ctypes.POINTER(ctypes.c_char_p)

# connect extract subchain lib
extract_subchain_lib_path = os.path.join(
    LOCATION, "modules", "extract_subchain_slib.so"
)
DEFAULT_CESAR = os.path.join(LOCATION, "CESAR2.0", "cesar")
ex_lib = ctypes.CDLL(extract_subchain_lib_path)
ex_lib.extract_subchain.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
ex_lib.extract_subchain.restype = ctypes.POINTER(ctypes.c_char_p)

# blosum matrix address
BLOSUM_FILE = os.path.join(LOCATION, "supply", "BLOSUM62.txt")
GAP_SCORE = 0
INS_PEN = -1
DEL_PEN = -1
AA_FLANK = 15
SS_SIZE = 2

# CESAR2.0 necessary files location
FIRST_CODON_PROFILE = "extra/tables/human/firstCodon_profile.txt"
LAST_CODON_PROFILE = "extra/tables/human/lastCodon_profile.txt"
ACC_PROFILE = "extra/tables/human/acc_profile.txt"
DO_PROFILE = "extra/tables/human/do_profile.txt"
EQ_ACC_PROFILE = os.path.join(LOCATION, "supply", "eq_acc_profile.txt")
EQ_DO_PROFILE = os.path.join(LOCATION, "supply", "eq_donor_profile.txt")


def verbose(msg):
    """Write verbose in stderr if verbose == TRUE."""
    eprint(msg) if VERBOSE else None


def parse_args():
    """Parse args, check and return."""
    # parse args
    app = argparse.ArgumentParser(description="")
    app.add_argument("gene", type=str, help="working gene")  # mandatory
    app.add_argument(
        "chains",
        type=str,
        help="chain ids, write 'region' if you like to type a region directly",
    )
    app.add_argument(
        "bdb_bed_file", type=str, help="BDB BED FILE, text bed12 also works but slower."
    )
    app.add_argument(
        "bdb_chain_file",
        type=str,
        help="BDB CHAIN FILE" "Or chain file, or even chain.gz file, but slower.",
    )
    app.add_argument("tDB", type=str, help="target db, alias or 2bit file")
    app.add_argument("qDB", type=str, help="query db, alias or 2bit file")
    app.add_argument(
        "--ref",
        default="hg38",
        help="Use this is you like to use " "chain alias and your ref is not the hg38.",
    )
    app.add_argument(
        "--cesar_binary", type=str, default=None, help="CESAR2.0 binary address."
    )
    app.add_argument(
        "--cesar_input_save_to",
        type=str,
        default=None,
        help="Where to save the cesar input file. /dev/shm/{whoami}/XXXXXX.fa is default.",
    )
    app.add_argument(
        "--memlim", default="Auto", help="memory limit for CESAR, auto as default."
    )
    app.add_argument(
        "--mask_stops",
        "--ms",
        action="store_true",
        dest="mask_stops",
        help="In case if there are inframe stop codons replace them with NNN.",
    )
    app.add_argument(
        "--estimate_memory",
        "-e",
        action="store_true",
        dest="estimate_memory",
        help="do not compute anything, just return what amount of memory is required for this run.",
    )
    app.add_argument(
        "--shift",
        type=int,
        default=2,
        help="num of blocks in up/downstream regions, 2 is default",
    )
    app.add_argument(
        "--raw_output",
        "--ro",
        default=None,
        help="Save raw CESAR output, " "You can use either stdout of a file for that.",
    )
    app.add_argument(
        "--cesar_output",
        "--co",
        default=None,
        help="If CESAR output is already done, specify it here",
    )
    app.add_argument(
        "--save_loci", default=None, help="Save gene: chain: locus dict to a file."
    )
    app.add_argument(
        "--output",
        default="stdout",
        help="Final output, stdout as default "
        "unreachable if params --raw_output or --estimate_memory are set.",
    )
    app.add_argument(
        "--prot_out", "--po", default="stdout", help="Save protein sequence"
    )
    app.add_argument(
        "--codon_out", "--cdo", default="stdout", help="Save codon alignment"
    )
    app.add_argument("--verbose", "-v", action="store_true", dest="verbose")
    app.add_argument(  # TODO: remove everything related to this parameter
        "--exon_flank",
        default=2,
        type=int,
        help="OBSOLETE When intersecting exons and intersected blocks"
        "add flanks to extend exons range.",
    )
    app.add_argument(
        "--gap_size",
        default=10,
        type=int,
        help="How many gaps in a row consider as an assembly gap.",
    )
    app.add_argument(
        "--gene_flank",
        default=1000,
        type=int,
        help="Up and downstream in query genome.",
    )
    app.add_argument(
        "--extra_flank",
        "--ef",
        default=0,
        type=int,
        help="Add bases up and down stream.",
    )
    app.add_argument(
        "--query_len_limit",
        "--ql",
        default=0,
        type=int,
        help="Skip job if query is too long",
    )
    app.add_argument(
        "--check_loss",
        "--cl",
        action="store_true",
        dest="check_loss",
        help="Check whether a gene found has inactivating mutations.",
    )
    app.add_argument(
        "--print_inact_mut",
        "--pim",
        action="store_true",
        dest="print_inact_mut",
        help="Print gene loss data on the screen",
    )
    app.add_argument(
        "--ic",
        action="store_true",
        dest="ic",
        help="Invert complement query sequences (DANGER)",
    )
    app.add_argument(
        "--paral",
        action="store_true",
        dest="paral",
        help="Mark that this gene is called for paralogous chains",
    )
    app.add_argument("--u12", default=None, help="Add U12 exons data")
    app.add_argument("--uhq_flank", default=50, type=int, help="UHQ flank filter")
    app.add_argument(
        "--no_fpi",
        action="store_true",
        dest="no_fpi",
        help="Consider some frame-preserving mutations as inactivating. "
        "See documentation for details.",
    )
    app.add_argument(
        "--fragments",
        "--fr",
        action="store_true",
        dest="fragments",
        help="Stitch scaffolds -> in case of a fragmented genome",
    )
    app.add_argument(
        "--precomputed_orth_loci",
        default=None,
        help="Path to loci saved with --save_orth_locus"
    )
    app.add_argument(
        "--do_not_check_exon_chain_intersection",
        default=False,
        action="store_true",
        dest="do_not_check_exon_chain_intersection",
        help="Do not extract chain blocks (not recommended)"
    )
    app.add_argument(
        "--alt_frame_del",
        "--lfd",
        default=False,
        action="store_true",
        dest="alt_frame_del",
        help="Consider codons in alternative frame (between compensated FS) deleted"
    )
    app.add_argument(
        "--mask_all_first_10p",
        "--m_f10p",
        action="store_true",
        dest="mask_all_first_10p",
        help="Automatically mask all inactivating mutations in first 10 percent of "
             "the reading frame, ignoring ATG codons distribution."
    )
    app.add_argument(
        "--temp_dir",
        default=None,
        help="Temp dir to store intermediate files.\n"
             "In TOGA pipeline: $toga_output/temp/cesar_temp_files; "
             "As a standalone script, by default: cesar_temp_files in the "
             "directory where it was called"
    )

    if len(sys.argv) == 1:
        app.print_help()
        sys.exit(0)
    # call the main
    ret_args = app.parse_args()

    # check args
    die(
        f"Error! --shift must be >=0, {ret_args.shift} given", 1
    ) if ret_args.shift < 0 else None
    die(
        f"Error! --exon_flank must be >=0 {ret_args.exon_flank} given", 1
    ) if ret_args.exon_flank < 0 else None
    die("Error! --gap_size must be >= 1!", 1) if ret_args.gap_size < 1 else None
    global VERBOSE
    VERBOSE = True if ret_args.verbose else False
    verbose("Program runs with the following params:")
    for k, v in vars(ret_args).items():
        verbose('"{0}": {1}'.format(k, v))
    return ret_args


def read_bed(gene, bed_file):
    """Extract and parse gene-related bed track."""
    try:  # if it's berkeley db file
        bed_track_raw = bed_extract_id(bed_file, gene)
    except OSError:
        # h5df cannot open this file: this is likely a text file
        bed_track_raw = bed_extract_id_text(bed_file, gene)

    # left CDS only
    bed_track = make_cds_track(bed_track_raw)

    # regular parsing of a bed-12 formatted file
    bed_info = bed_track.split("\t")

    # if not 12 fields -> then input is wrong
    die("Error! Bed-12 required!", 1) if len(bed_info) != 12 else None
    bed_data = {
        "chrom": bed_info[0],
        "chromStart": int(bed_info[1]),
        "chromEnd": int(bed_info[2]),
        "name": bed_info[3],
        "strand": True if bed_info[5] == "+" else False,
    }
    block_sizes = [int(x) for x in bed_info[10].split(",") if x != ""]
    block_starts = [
        int(x) + bed_data["chromStart"] for x in bed_info[11].split(",") if x != ""
    ]
    block_ends = [x[0] + x[1] for x in zip(block_starts, block_sizes)]
    bed_data["blocks"] = list(zip(block_starts, block_ends))
    bed_data["block_sizes"] = (
        block_sizes[::-1] if not bed_data["strand"] else block_sizes
    )
    verbose("Bed data is:\n{0}".format(bed_data))
    return bed_data


def revert(line):
    """Revert string and change with complement bases."""
    line = line[::-1].upper()
    new_str = "".join([COMPLEMENT_BASE.get(c) for c in line if COMPLEMENT_BASE.get(c)])
    return new_str


def get_2bit_path(db_arg):
    """Check if alias and return a path to 2bit file."""
    if os.path.isfile(db_arg):  # not an alias
        return db_arg  # there is nothing to do
    elif os.path.islink(db_arg):
        return os.readlink(db_arg)
    aliased = two_bit_templ.format(db_arg)
    # check that it's a file
    die(f"Error! Cannot find {aliased} file", 1) if not os.path.isfile(
        aliased
    ) else None
    return aliased


def find_chain_file(ref_name, chain_arg):
    """Return chain file path."""
    if os.path.isfile(chain_arg):
        # not an alias -> just return it
        return chain_arg
    # not a chain path/ Hillerlab-specific place!
    chain_path = chain_alias_template.format(ref_name, chain_arg)
    die(
        f"Error! File {chain_path} not found! Please set the chain path explicitly.", 1
    ) if not os.path.isfile(chain_path) else None
    return chain_path


def get_exons(bed_data, t_db):
    """Extract exons sequences for reference."""
    exons_raw = (
        bed_data["blocks"][::-1] if not bed_data["strand"] else bed_data["blocks"]
    )
    exons_pos = {}  # contain exon_num : positions
    s_sites = []
    exon_flanks = {}
    for num, exon in enumerate(exons_raw):
        # start, end if strand == + end, start otherwise
        # need to know start and end to extract from 2bit file
        exons_pos[num] = (
            (int(exon[0]), int(exon[1]))
            if bed_data["strand"]
            else (int(exon[1]), int(exon[0]))
        )

    max_exon_num = max(exons_pos.keys())
    _all_positions = sorted(flatten(exons_pos.values()))
    gene_borders = {_all_positions[0], _all_positions[-1]}
    # extract sequences
    exons_seq = {}  # exon number: sequence dict
    target_genome = TwoBitFile(get_2bit_path(t_db))  # use 2bitreader library
    get_chr = bed_data["chrom"]
    try:
        chrom_seq = target_genome[get_chr]
    except KeyError:
        chrom_seq = []  # to suppress PyCharm analyzer
        die(f"Error! Cannot find chrom {get_chr} in 2bit file {t_db}")
    verbose("\nExons sequences ####\n")
    for num, pos in exons_pos.items():
        is_first_exon = num == 0
        is_last_exon = num == max_exon_num

        # for twoBitToFa start must be < end
        # determine search start and end
        # do not subtract/add SS_SIZE if gene border: no splice sites then
        min_pos = min(pos)
        max_pos = max(pos)
        start = min_pos - SS_SIZE if min_pos not in gene_borders else min_pos
        end = max_pos + SS_SIZE if max_pos not in gene_borders else max_pos
        # get exon 10 bp flanks:
        left_brd_ = min_pos - EXON_SEQ_FLANK if min_pos - EXON_SEQ_FLANK > 0 else 0
        left_flank_coord = (left_brd_, min_pos)
        right_flank_coord = (max_pos, max_pos + EXON_SEQ_FLANK)
        left_flank = chrom_seq[left_flank_coord[0]: left_flank_coord[1]].upper()
        right_flank = chrom_seq[right_flank_coord[0]: right_flank_coord[1]].upper()
        # correct for strand:
        left_flank = left_flank if bed_data["strand"] else revert(left_flank)
        right_flank = right_flank if bed_data["strand"] else revert(right_flank)
        if not bed_data["strand"]:
            left_flank, right_flank = right_flank, left_flank
        # placeholder in case we coudld not extract flanks
        left_flank = left_flank if len(left_flank) > 0 else "X"
        right_flank = right_flank if len(right_flank) > 0 else "X"

        exon_seq_raw = chrom_seq[start:end].upper()
        # revert if negative strand
        exon_seq_w_ss = revert(exon_seq_raw) if not bed_data["strand"] else exon_seq_raw
        # trim splice sites
        if not is_first_exon and not is_last_exon:
            # both splice sites are here
            exon_seq = exon_seq_w_ss[SS_SIZE:-SS_SIZE]
            acc_ = exon_seq_w_ss[:SS_SIZE]
            don_ = exon_seq_w_ss[-SS_SIZE:]
        elif is_first_exon and is_last_exon:
            # no splice sites
            exon_seq = exon_seq_w_ss
            acc_, don_ = "NN", "NN"
        elif is_first_exon:
            exon_seq = exon_seq_w_ss[:-SS_SIZE]
            acc_ = "NN"
            don_ = exon_seq_w_ss[-SS_SIZE:]
        elif is_last_exon:
            exon_seq = exon_seq_w_ss[SS_SIZE:]
            acc_ = exon_seq_w_ss[:SS_SIZE]
            don_ = "NN"
        else:
            raise RuntimeError("Unreachable branch reached")
        s_sites.append(acc_)
        s_sites.append(don_)
        exons_seq[num] = exon_seq
        exon_flanks[num] = {"L": left_flank, "R": right_flank}
        verbose(f"Exon {num} in the range {pos}; sequence:\n{exon_seq}")
    return exons_pos, exons_seq, s_sites, exon_flanks


def check_ref_exons(exon_seqs, mask_stops):
    """Check if the reference sequence is correct.

    Should start with ATG and end with a stop.
    Mask_stops controls handling of inframe stops.
    """
    sec_codons = set()  # in case there are TGA codons in the ref seq -> collect them
    gene_seq = "".join([exon_seqs[i] for i in range(len(exon_seqs.keys()))])
    codons = parts(gene_seq, n=3)  # split a seq of letters in chunks of len == 3
    if codons[0] != "ATG":
        eprint("Input is corrupted! Reference sequence should start with ATG!")
    elif codons[-1] not in Constants.STOP_CODONS:
        eprint("Input is corrupted! Reference sequence should end with a stop codon!")
    stop_codons = [(n, c) for n, c in enumerate(codons[:-1]) if c in Constants.STOP_CODONS]
    if len(stop_codons) == 0:  # no stop codons -> nothing else to do
        return exon_seqs, set()
    # there are stop codons in reference sequence:
    eprint("Warning! There are inframe stop codons!")
    for stop in stop_codons:
        eprint(f"Codon num {stop[0] + 1} - {stop[1]}")
        codons[stop[0]] = Constants.NNN_CODON if mask_stops else codons[stop[0]]
        if stop[1] == "TGA":
            # maybe a sec codon
            sec_codons.add(stop[0])

    eprint(">>>STOP_CODON>>>") if not mask_stops else None
    die("Abort, there are inframe stop codons.", 0) if not mask_stops else None
    # if stop codons in reference are allowed, then we need to mask them (rewrite as NNN)
    # otherwise CESAR will show an error
    safe_seq = "".join(codons)
    stop_masked = {}
    prev_index = 0
    for num, exon_seq in exon_seqs.items():
        exon_len = len(exon_seq)
        stop_masked[num] = safe_seq[prev_index: prev_index + exon_len]
        prev_index += exon_len
    return stop_masked, sec_codons


def prepare_exons_for_cesar(exon_seqs):
    """Cesar requires some special formatting for exons."""
    # CESAR requires a special king of reference exons formatting
    # Split codons should be written in lowercase, as follows:
    # ATGTTTa ctGTAAAGTGCc ttAGTTGA
    verbose("prepare_exons_for_cesar")
    left_pointer = 0  # init value
    cesar_input_exons = {}  # accumulate result here

    for k, exon_seq in exon_seqs.items():
        # define number of letters to lowercase at the right side
        right_pointer = (len(exon_seq) - left_pointer) % 3
        # apply left pointer
        if left_pointer != 3:  # it it were accumulated 3 bases it has no sense like 0
            exon_seq_lfix = (
                exon_seq[:left_pointer].lower() + exon_seq[left_pointer:].upper()
            )
        else:  # don't touch if left_pointer == 0 or == 3
            exon_seq_lfix = exon_seq
        # re-define left pointer
        left_pointer = 3 - right_pointer
        # apply right-side pointer
        if right_pointer != 0:
            exon_seq_rfix = (
                exon_seq_lfix[:-right_pointer] + exon_seq_lfix[-right_pointer:].lower()
            )
        else:  # save as it was
            exon_seq_rfix = exon_seq_lfix
        # save prepared exon into special dict
        cesar_input_exons[k] = exon_seq_rfix
    return cesar_input_exons


def get_chain(chain_file, chain_id):
    """Return chain string according the parameters passed."""
    chain = None  # to calm IDE down
    if chain_file.endswith(".bst"):
        # we have bdb file; extract with BDB extractor
        chain = chain_extract_id(chain_file, chain_id)
        return chain
    elif chain_file.endswith(".gz"):  # a gzipped chain file was given
        # gzip and redirect scteam to chain_filter_by_id binary
        extract_by_id_cmd = (
            f"gzip -dc {chain_file} | ./modules/chain_filter_by_id stdin {chain_id}"
        )
        try:  # check that output is OK
            chain = subprocess.check_output(extract_by_id_cmd, shell=True).decode(
                "utf-8"
            )
        except subprocess.CalledProcessError:
            # die if the command died
            die(
                "Error! Process {extract_by_id_cmd} died! Please check if input data is correct",
                1,
            )
        return chain

    else:  # just a chain file, extract the chain we need
        # the same as above, but without gzip stage
        extract_by_id_cmd = f"./modules/chain_filter_by_id {chain_file} {chain_id}"
        try:  # also check if output is OK, die otherwise
            chain = subprocess.check_output(extract_by_id_cmd, shell=True).decode(
                "utf-8"
            )
        except subprocess.CalledProcessError:
            die(
                "Error! Process {extract_by_id_cmd} died! Please check if input data is correct",
                1,
            )
        return chain


def range_corrector(g_range):
    """Swap start and end if start > end."""
    chrom, start_end = g_range.split(":")
    start_end_split = start_end.split("-")
    start, end = int(start_end_split[0]), int(start_end_split[1])
    if start < end:
        return g_range
    else:
        return f"{chrom}:{end}-{start}"


def chain_cut(chain_str, gene_range, gene_flank, extra_flank=0):
    """Call chain_cut binary.

    Project reference gene coordinates to query through a chain.
    Also add flanks if shift is > 0.
    """
    # need to get genomic region for the gene
    # also need to translate python data types to C
    # to call the shared library; I do it 2 times here
    # for shift = 0 and shifts = 2 (add flanks around gene)
    c_chain = ctypes.c_char_p(chain_str.encode())
    c_shift_2 = ctypes.c_int(2)
    c_shift_0 = ctypes.c_int(0)
    granges_num = 1
    c_granges_num = ctypes.c_int(granges_num)  # we need only one grange to analyze
    granges_arr = (ctypes.c_char_p * (granges_num + 1))()  # granges_num + 1
    granges_bytes = [gene_range.encode("utf-8")]
    # need to do this tricks to pass strings array to C
    granges_arr[:-1] = granges_bytes
    granges_arr[granges_num] = None
    raw_ch_conv_s2 = ch_lib.chain_coords_converter(
        c_chain, c_shift_2, c_granges_num, granges_arr
    )
    chain_coords_conv_out_s2 = []  # keep lines here
    # convert C output to python-readable type
    for i in range(granges_num + 1):
        chain_coords_conv_out_s2.append(raw_ch_conv_s2[i].decode("utf-8"))
    # chain', 'chr5', '+', '137889395', '148245211', 'chr18', '+', '34409342', '44120958
    chain_data = chain_coords_conv_out_s2[0].split(" ")
    t_strand = True if chain_data[2] == "+" else False
    q_strand = True if chain_data[7] == "+" else False
    t_size = int(chain_data[3])
    q_size = int(chain_data[8])

    # re-define arrays to avoid segfault
    c_chain = ctypes.c_char_p(chain_str.encode())
    granges_arr = (ctypes.c_char_p * (granges_num + 1))()  # granges_num + 1
    granges_bytes = [gene_range.encode("utf-8")]
    granges_arr[:-1] = granges_bytes
    granges_arr[granges_num] = None

    raw_ch_conv_s0 = ch_lib.chain_coords_converter(
        c_chain, c_shift_0, c_granges_num, granges_arr
    )
    chain_coords_conv_out_s0 = []  # keep lines here

    # convert C output to python-readable type
    for i in range(granges_num + 1):
        chain_coords_conv_out_s0.append(raw_ch_conv_s0[i].decode("utf-8"))
    # another approach to detect range
    # sometimes blocks go so far
    # ------------------genegene-------------------
    # block-------------blockblock------------block

    # to avoid very huge query sequences program controls it's size
    search_region_shift_str = range_corrector(
        chain_coords_conv_out_s2[1].split("\t")[1]
    )
    search_region_abs_str = range_corrector(chain_coords_conv_out_s0[1].split("\t")[1])

    chrom = search_region_shift_str.split(":")[0]
    search_reg_shift = [
        int(x) for x in search_region_shift_str.split(":")[1].split("-")
    ]
    search_reg_abs = [int(x) for x in search_region_abs_str.split(":")[1].split("-")]
    search_reg_flanked = [
        search_reg_abs[0] - gene_flank,
        search_reg_abs[1] + gene_flank,
    ]

    # define actual starts and ends
    act_start = (
        search_reg_shift[0]
        if search_reg_shift[0] > search_reg_flanked[0]
        else search_reg_flanked[0]
    )
    act_end = (
        search_reg_shift[1]
        if search_reg_shift[1] < search_reg_flanked[1]
        else search_reg_flanked[1]
    )

    if extra_flank > 0:  # add extra flanks if required
        act_start = act_start - extra_flank if act_start - extra_flank > 0 else 0
        act_end = (
            act_end + extra_flank if act_end + extra_flank < q_size else q_size - 1
        )

    act_search_range = f"{chrom}:{act_start}-{act_end}"
    # ext_search_range = f"{chrom}:{}-{}"
    del raw_ch_conv_s0  # not sure if this is necessary
    del raw_ch_conv_s2  # but just in case
    return (
        act_search_range,
        search_region_shift_str,
        (t_strand, t_size, q_strand, q_size),
    )


def extract_subchain(chain_str, search_locus):
    """Extract subchain containing only the search locus."""
    # convert python types to C types
    c_chain = ctypes.c_char_p(chain_str.encode())
    c_mode = ctypes.c_char_p("q".encode())
    c_search = ctypes.c_char_p(search_locus.encode())
    # call the library finally
    ex_out_raw = ex_lib.extract_subchain(c_chain, c_mode, c_search)
    ex_out = []

    for elem in ex_out_raw:
        line = elem.decode("utf-8")
        if line == "END":
            break
        ex_out.append(line)

    blocks = [[int(x) for x in elem.split()] for elem in ex_out]
    if len(blocks) == 0:  # die an show chain id
        chain_id = chain_str.split("\n")[0].split()[-1]
        die(f"Error! No overlapping blocks for chain {chain_id} found!", 1)
    del ex_out_raw  # maybe it's necessary
    return blocks


def orient_blocks(subchain_blocks_raw, chain_data):
    """Create block num: coordinates dict.

    Orient them in correct direction.
    Add interblock regions, like block 1_2 between blocks 1 and 2.
    """
    block_ranges = {}
    t_strand, t_size, q_strand, q_size = chain_data

    for i in range(len(subchain_blocks_raw)):
        i_block = subchain_blocks_raw[i]
        t_start = i_block[0] if t_strand else t_size - i_block[1]
        t_end = i_block[1] if t_strand else t_size - i_block[0]
        q_start = i_block[2] if q_strand else q_size - i_block[3]
        q_end = i_block[3] if q_strand else q_size - i_block[2]
        block_ranges[str(i)] = (
            min(t_start, t_end),
            max(t_start, t_end),
            min(q_start, q_end),
            max(q_start, q_end),
        )

    # if there is only one block (weird but possible) --> no gaps between blocks
    if len(block_ranges) == 1:
        return block_ranges

    # make interblock gaps
    blocks_num = len(block_ranges.keys())
    for i in range(1, blocks_num):
        prev = block_ranges[str(i - 1)]
        current = block_ranges[str(i)]
        # connect xEnd_prev to xStart_current
        inter_t_start = prev[1] if t_strand else prev[0]
        inter_t_end = current[0] if t_strand else current[1]
        inter_q_start = prev[3] if q_strand else prev[2]
        inter_q_end = current[2] if q_strand else current[3]
        interblock_data = (inter_t_start, inter_t_end, inter_q_start, inter_q_end)
        block_ranges[f"{i - 1}_{i}"] = interblock_data
    return block_ranges


def intersect_ranges(range_1, range_2):
    """If > 0: ranges intersect."""
    range_1 = range_1 if range_1[1] > range_1[0] else range_1[::-1]
    range_2 = range_2 if range_2[1] > range_2[0] else range_2[::-1]
    return min(range_1[1], range_2[1]) - max(range_1[0], range_2[0])


def get_aa_ex_cov(exon_range, i_block_coords):
    """Get number of bases covered in target and query."""
    # deal with t_cov first
    # chain blocks without _ in name - real blocks, len in T == len in Q
    # blocks with _ -> not aligned
    # will get exon len + flanks - len of blocks with _
    flanked_exon_len = max(exon_range) - min(exon_range)
    exon_range = (min(exon_range), max(exon_range))
    interblock_reg = [x[1] for x in i_block_coords if "_" in x[0]]
    if len(interblock_reg) == 0:
        # single block covers the exon + flanks: ideal case
        return flanked_exon_len, flanked_exon_len, flanked_exon_len
    # select blocks, ignore
    cov_blocks_not_trimmed = [x[1] for x in i_block_coords if "_" not in x[0]]
    t_cov_blocks_trimmed = []
    # trim blocks | convert this:
    # --------exonexonexonexonexon--------
    # --blockclobk===blockblock===========
    # to:
    # --------exonexonexonexonexon--------
    # --------block==blockblock===========
    for block in cov_blocks_not_trimmed:
        t_start, t_end = block[0], block[1]
        t_block_vals = (t_start, t_end)
        t_block = [min(t_block_vals), max(t_block_vals)]
        start_within = exon_range[0] <= t_block[0] <= exon_range[1]
        end_within = exon_range[0] <= t_block[1] <= exon_range[1]
        if not start_within and not end_within:
            # something impossible
            continue
            # die("Impossible block configuration in get_aa_ex_cov!")
        if start_within and end_within:
            t_cov_blocks_trimmed.append(t_block)
            continue
        elif not start_within:
            t_block[0] = exon_range[0]
            t_cov_blocks_trimmed.append(t_block)
            continue
        else:
            t_block[1] = exon_range[1]
            t_cov_blocks_trimmed.append(t_block)
            continue
    t_cov = sum([abs(x[1] - x[0]) for x in t_cov_blocks_trimmed])
    t_cov = t_cov if t_cov > 0 else 0

    # now compute q_cov
    q_blocks_within = []
    for block in interblock_reg:
        t_start, t_end = block[0], block[1]
        q_start, q_end = block[2], block[3]
        q_block_vals = (q_start, q_end)
        t_block_vals = (t_start, t_end)
        t_block = [min(t_block_vals), max(t_block_vals)]
        q_block = [min(q_block_vals), max(q_block_vals)]
        start_within = exon_range[0] <= t_block[0] <= exon_range[1]
        end_within = exon_range[0] <= t_block[1] <= exon_range[1]
        if not start_within or not end_within:
            continue
        q_blocks_within.append(q_block)
    q_corr_reg = sum([abs(x[1] - x[0]) for x in q_blocks_within])
    q_cov = t_cov + q_corr_reg
    q_cov = q_cov if q_cov > 0 else 0
    return flanked_exon_len, t_cov, q_cov


def intersect_exons_blocks_gaps(
    exon_coordinates, subchain_blocks, gap_coordinates, flank, uhq_flank
):
    """Intersect exons, chain blocks and gaps.

    Create the following dictionaries:
    exon_num: intersected chain blocks (list)
    chain_block: intersects gap or not (bool)
    """
    # make a pseudoblock --> for the entire chain
    # all T - all chain blocks in reference coordinates
    all_t = [b[0] for b in subchain_blocks.values()] + [
        b[1] for b in subchain_blocks.values()
    ]
    # all Q - all chain blocks in query coordinates
    all_q = [b[2] for b in subchain_blocks.values()] + [
        b[3] for b in subchain_blocks.values()
    ]
    global_block = (min(all_t), max(all_t), min(all_q), max(all_q))
    # get a dict of exons with flanks
    exon_flank_coordinates = {}
    exon_aa_flank_coordinates = {}  # just to check that AA criteria satisfied
    marginal_cases = {e: False for e in exon_coordinates.keys()}

    for exon_num, exon_coords in exon_coordinates.items():
        # get exon -> flanked coordinates
        exon_start, exon_end = min(exon_coords), max(exon_coords)
        exon_size = exon_end - exon_start
        flanked_exon_start = exon_start - flank
        flanked_exon_end = exon_end + flank
        aa_flanked_start = exon_start - uhq_flank
        aa_flanked_end = exon_end + uhq_flank
        exon_flank_coordinates[exon_num] = (flanked_exon_start, flanked_exon_end)
        exon_aa_flank_coordinates[exon_num] = (aa_flanked_start, aa_flanked_end)
        # also check if a block outside the chain
        glob_intersect = intersect_ranges(
            exon_coords, (global_block[0], global_block[1])
        )
        if glob_intersect != exon_size:
            # exon doesn't intersect chain -> put to marginal cases
            marginal_cases[exon_num] = True

    # find intersections for blocks and assembly gaps
    exon_num_blocks = defaultdict(list)
    fl_exon_num_blocks = defaultdict(list)
    aa_exon_num_blocks = defaultdict(list)
    block_gaps = {b: False for b in subchain_blocks.keys()}

    for block_num, block_coords in subchain_blocks.items():
        # try intersection with all exons
        for exon_num, exon_coords in exon_coordinates.items():
            flanked_coords = exon_flank_coordinates.get(exon_num)
            aa_flank_coords = exon_aa_flank_coordinates.get(exon_num)
            intersect_norm = intersect_ranges(
                exon_coords, (block_coords[0], block_coords[1])
            )
            intersect_flanked = intersect_ranges(
                flanked_coords, (block_coords[0], block_coords[1])
            )
            intersect_aa_class = intersect_ranges(
                aa_flank_coords, (block_coords[0], block_coords[1])
            )
            if intersect_norm >= 0:
                exon_num_blocks[exon_num].append(block_num)
            if intersect_flanked >= 0:
                fl_exon_num_blocks[exon_num].append(block_num)
            if intersect_aa_class >= 0:
                aa_exon_num_blocks[exon_num].append(block_num)
        # and gaps
        for gap_num, gap in enumerate(gap_coordinates):
            intersect = intersect_ranges(gap, (block_coords[2], block_coords[3]))
            if intersect <= 0:
                continue
            block_gaps[block_num] = gap_num
            verbose_msg = (
                f"Block num {block_num} in coords t:{block_coords[0]}-{block_coords[1]} "
                f"q:{block_coords[2]}-{block_coords[3]} intersects gap {gap[0]}-{gap[1]}"
            )
            verbose(verbose_msg)
    # get missed exons to exclude
    missing_exons = [e for e in exon_coordinates.keys() if not exon_num_blocks.get(e)]
    not_covered_str = ", ".join([str(e) for e in missing_exons])
    verbose(f"Exons:\n{not_covered_str}\nare not covered by chain") if len(missing_exons) > 0 else None
    verbose(f"Flank size is: {flank}")

    for exon_num in sorted(exon_num_blocks.keys()):
        verbose(f"Exon {exon_num} intersects blocks {exon_num_blocks.get(exon_num)}")
        verbose(
            f"Exon {exon_num} with flanks intersects blocks {fl_exon_num_blocks.get(exon_num)}"
        )
        verbose(
            f"Exon {exon_num} with AA flanks intersects blocks {aa_exon_num_blocks.get(exon_num)}"
        )
        verbose("\n")
    aa_sat = {}

    for exon_num, blocks in aa_exon_num_blocks.items():
        exon_range = exon_aa_flank_coordinates[exon_num]
        i_block_coords = [(b, subchain_blocks[b]) for b in blocks]
        fex_len, t_cov, q_cov = get_aa_ex_cov(exon_range, i_block_coords)
        # allow 20% deviation
        # if corresponding regions are in this range: it's OK
        dev_thr = uhq_flank * 0.2
        lo_, hi_ = fex_len - dev_thr, fex_len + dev_thr
        t_cov_ok = lo_ <= t_cov <= hi_
        q_cov_ok = lo_ <= q_cov <= hi_
        cond_sat_ = t_cov_ok is True and q_cov_ok is True
        aa_sat[exon_num] = True if cond_sat_ else False
    return (
        exon_num_blocks,
        fl_exon_num_blocks,
        block_gaps,
        missing_exons,
        exon_flank_coordinates,
        marginal_cases,
        aa_sat,
    )


def sort_blocks(blocks_num):
    """Return block: index and index: block dicts."""
    block_to_index, index_to_block = {}, {}
    pointer = -1
    for i in range(blocks_num):
        if i % 2 == 0:
            # if so, it is a real block
            pointer += 1
            block_id = str(pointer)
        else:  # otherwise an interblock range
            block_id = f"{pointer}_{pointer + 1}"
        block_to_index[block_id] = i
        index_to_block[i] = block_id
    return block_to_index, index_to_block


def classify_predict_exons(exon_blocks, block_coordinates, margin_cases):
    """Classify exons and get expected coordinates."""
    exon_class, exon_exp_q_region = {}, {}
    block_deleted = {}
    block_to_index, index_to_block = sort_blocks(len(block_coordinates))
    # for each block compute if it corresponds to gap in query
    for block_num, block_coords in block_coordinates.items():
        block_deleted[block_num] = True if block_coords[2] == block_coords[3] else False

    for exon_num, blocks in exon_blocks.items():
        # extract expected query region
        block_indexes = sorted([block_to_index.get(b) for b in blocks])
        start_block_id = index_to_block.get(block_indexes[0])
        end_block_id = index_to_block.get(block_indexes[-1])
        start_block_coords = block_coordinates.get(start_block_id)
        end_block_coords = block_coordinates.get(end_block_id)
        exp_q_region_start = min(start_block_coords[2], end_block_coords[2])
        exp_q_region_end = max(start_block_coords[3], end_block_coords[3])
        exon_exp_q_region[exon_num] = (exp_q_region_start, exp_q_region_end)
        margin_exon = margin_cases.get(exon_num)

        # assign the class
        if margin_exon:
            # in case of margin exon the following happens:
            #      exonexonexon
            #          chainchainchain===chainchain
            # so exon is not covered by the chain partly
            this_exon_class = "M"
        elif "_" not in start_block_id and "_" not in end_block_id:
            # both start and end of the exon are covered, class A
            # ----------exonexonexon----------
            # ---chainchain====chaincha=======
            # ==========chainchainchain-------
            this_exon_class = "A"
        # >ENST00000464162 | 1 | 4482 | KK498653:40658904-40658593 | 59.81 | OK | A |
        # exp:40655721-40655810 | EXCL
        elif start_block_id == end_block_id and "_" in start_block_id:
            # exon unaligned or completely removed, class C
            # ----------exonexonexon----------
            # ================================
            # --------------------------------
            this_exon_class = "C"
        else:
            # exon is covered somehow
            # ----------exonexonexon----------
            # -------chainchain===============
            # -------------chainch------------
            this_exon_class = "B"
        # three options
        exon_class[exon_num] = this_exon_class
        exon_exp_q_region[exon_num] = (exp_q_region_start, exp_q_region_end)

    return exon_class, exon_exp_q_region


def find_exons_gaps(
    fl_exon_coordinates,
    flanked_exon_blocks,
    block_coordinates,
    blocks_gaps,
    gap_coordinates,
):
    """For each flanked exon define if it overlaps gaps."""
    exon_gap = {}
    block_to_index, index_to_block = sort_blocks(len(block_coordinates))
    for exon_num, exon_coordinates in fl_exon_coordinates.items():
        intersected_blocks = flanked_exon_blocks.get(exon_num)
        if not intersected_blocks:  # exon not covered
            continue
        # find out if there are gaps intersected
        # asm_gaps = [blocks_gaps.get(b) for b in intersected_blocks]
        asm_gaps = {
            b: blocks_gaps.get(b)
            for b in intersected_blocks
            if blocks_gaps.get(b) is not False
        }
        gap_numbers = list(set(asm_gaps.values()))
        if len(gap_numbers) == 0:  # exon doesn't intersect any asm gaps
            exon_gap[exon_num] = False
            continue
        elif len(intersected_blocks) == 1 and "_" in intersected_blocks[0]:
            # only one block intersecting the gap
            # only one possibility to observe this:
            #              exonexonexon
            # block=====================block
            #       gapgap
            exon_gap[exon_num] = True
            continue
        elif len(intersected_blocks) == 1:
            # should not happen
            exon_gap[exon_num] = False
            continue
        # risky zone; more sophisticated cases
        block_indexes = sorted([block_to_index.get(b) for b in intersected_blocks])
        first_block, last_block = index_to_block.get(
            block_indexes[0]
        ), index_to_block.get(block_indexes[-1])
        if "_" not in first_block and "_" not in last_block:
            # it is possible in only the one case
            #      exonexonexonexonexonexon
            # blockblock==========block==block---
            #              gapgap
            exon_gap[exon_num] = True
            continue
        # check if gap inside
        not_firsts = [gapped_block != first_block for gapped_block in asm_gaps.keys()]
        not_lasts = [gapped_block != last_block for gapped_block in asm_gaps.keys()]
        if True in not_firsts or True in not_lasts:
            exon_gap[exon_num] = True
            continue
        # the most complicated case
        # exonexon---| flank
        # blockblock============
        #             gapgap
        gap_coords_to_check = [gap_coordinates[i] for i in gap_numbers]
        f_block_coords, l_block_coords = block_coordinates.get(
            first_block
        ), block_coordinates.get(last_block)
        f_block_delta = intersect_ranges(
            exon_coordinates, (f_block_coords[0], f_block_coords[1])
        )
        l_block_delta = intersect_ranges(
            exon_coordinates, (l_block_coords[0], l_block_coords[1])
        )
        # we need the direction of blocks
        # 1 - 2 - 3 - 4 - 5 OR
        # 5 - 4 - 3 - 2 - 1 ?
        direction = True if f_block_coords[0] < l_block_coords[1] else False
        exon_gap[exon_num] = False  # default value
        if direction:
            f_block_q_end = max(f_block_coords[2], f_block_coords[3])
            l_block_q_start = min(l_block_coords[2], l_block_coords[3])
            f_block_q_ext = (f_block_q_end - f_block_delta, f_block_q_end)
            l_block_q_ext = (l_block_q_start, l_block_q_start + l_block_delta)
        else:
            l_block_q_end = max(l_block_coords[2], l_block_coords[3])
            f_block_q_start = min(f_block_coords[2], f_block_coords[3])
            l_block_q_ext = (l_block_q_end - l_block_delta, l_block_q_end)
            f_block_q_ext = (f_block_q_start, f_block_q_start + f_block_delta)

        for gap_coord in gap_coords_to_check:
            f_intersect = intersect_ranges(gap_coord, l_block_q_ext)
            l_intersect = intersect_ranges(gap_coord, f_block_q_ext)
            if f_intersect > 0 or l_intersect > 0:
                exon_gap[exon_num] = True
    return exon_gap


def make_query_seq(chain_id, search_locus, q_db, chain_strand, bed_strand):
    """Extract query sequence."""
    query_genome = TwoBitFile(get_2bit_path(q_db))
    q_name, region = search_locus.split(":")
    query_start = int(region.split("-")[0])
    query_end = int(region.split("-")[1])
    chrom_seq = query_genome[q_name]
    # find the locus
    verbose(
        f"Looking for exons in range: {q_name}:{query_start}-{query_end} for chain {chain_id}"
    )
    query_seq_no_dir = chrom_seq[query_start:query_end]
    directed = chain_strand == bed_strand
    query_seq = query_seq_no_dir if directed else revert(query_seq_no_dir)
    verbose(f"Query length for {chain_id} in locus {search_locus}: {len(query_seq)}")
    die(f"Error! No query sequence for chain {chain_id}!", 1) if len(
        query_seq
    ) == 0 else None
    return query_seq, directed


def find_gaps(query_seq, search_locus, gap_size, directed):
    """Return ranges of gaps."""
    pattern = r"N{" + str(gap_size) + ",}"
    seq_start = int(search_locus.split(":")[1].split("-")[0])
    query_seq_corr = query_seq if directed else revert(query_seq)
    gap_ranges = []
    for match in finditer(pattern, query_seq_corr, IGNORECASE):
        span_start, span_end = match.span()
        gap_ranges.append((seq_start + span_start, seq_start + span_end))
    verbose("\n# Gap ranges are:") if len(gap_ranges) > 0 else None
    for gap_range in gap_ranges:
        verbose(f"{gap_range[0]}-{gap_range[1]}")
    verbose("\n# No assembly gaps detected") if len(gap_ranges) == 0 else None
    return gap_ranges


def make_in_filename(cesar_in, temp_dir):
    """Make filename for the CESAR input."""
    # TODO: refactor handling temp files
    if cesar_in:  # is defined by user
        # cesar in - file addr, False --> to mark it is not temp
        return cesar_in, False
    os.mkdir(temp_dir) if not os.path.isdir(temp_dir) else None
    filename = f"{uuid.uuid4()}.fa"
    cesar_in_path = os.path.join(temp_dir, filename)
    # mark that it is a temp file with True
    return cesar_in_path, True


def make_cesar_in(exons, queries, u12_elems, cesar_in_filename, cesar_temp_dir):
    """Make input for CESAR."""
    cesar_in_filename, is_temp = make_in_filename(cesar_in_filename, cesar_temp_dir)
    input_line = ""  # init var

    # write exons first
    exons_num = len(exons.keys())
    exons_last_ind = exons_num - 1

    for num, exon in exons.items():
        if exons_num == 1:
            # if single exon -> don't do anything with profiles
            header = f">exon_{num}"
            input_line += f"{header}\n{exon}\n"
            continue
        # not a single exon gene
        num_corr = num + 1
        acc_elem = (num_corr, 0)
        don_elem = (num_corr, 1)
        acc = ACC_PROFILE if acc_elem not in u12_elems else EQ_ACC_PROFILE
        don = DO_PROFILE if don_elem not in u12_elems else EQ_DO_PROFILE
        # if first or last: not a splice site actually
        acc = acc if num != 0 else FIRST_CODON_PROFILE
        don = don if num != exons_last_ind else LAST_CODON_PROFILE
        header = f">exon_{num}\t{acc}\t{don}"
        input_line += f"{header}\n{exon}\n"

    input_line += "####\n"  # required by CESAR to separate target and query with ####
    # write queries
    for chain_id, chain_seq in queries.items():
        input_line += f">{chain_id}\n{chain_seq}\n"

    if is_temp is True:
        # in TOGAs < 1.1.4 - we'd create a temp file in the /dev/shm
        # now, we don't do it - just pass to the CESAR /dev/stdin
        return input_line, None, True
    
    # user asked to create a specific file to store cesar's input
    with open(cesar_in_filename, "w") as f:
        f.write(input_line)
    return None, cesar_in_filename, False


def get_exon_indexes(exon_sizes):
    """For an input file given return an array of letter indexes in exons."""
    exon_indexes = {}
    exon_num, prev_max = 0, 0
    for exon_size in exon_sizes:
        inds = list(range(prev_max, prev_max + exon_size))
        exon_indexes.update({ind: exon_num for ind in inds})
        exon_num += 1
        prev_max += exon_size
    return exon_indexes


def memory_check(block_sizes, qlength_max, estimate_memory):
    """Return amount of memory required.

    Estimate memory consumption exactly as CESAR2.0 doest it.
    """
    num_states, rlength, extra = 0, 0, 100000
    for block_size in block_sizes:
        num_codons = block_size // 3
        num_states += 6 + 6 * num_codons + 1 + 2 + 2 + 22 + 6
        rlength += block_size
    mem = (
        (num_states * 4 * 8)
        + (num_states * qlength_max * 4)
        + (num_states * 304)
        + (2 * qlength_max + rlength) * 8
        + (qlength_max + rlength) * 2 * 1
        + extra
    )
    # bytes to GB
    mem_gb = mem / 1000000000
    if estimate_memory:
        sys.stdout.write(f"Expected memory consumption of:\n{mem_gb} GB\n")
        sys.exit(0)
    return mem_gb


def run_cesar(
    input_file,
    input_data,
    memory_raw,
    memlim,
    cesar_binary
):
    """Run CESAR for the input file or data.

    Normally, runs using input_data -> a string representing the input data for CESAR
    Alternatively, can be executed using input_file -> temp file created to call CESAR.
    The second option may be useful if CESAR_wrapper.py is executed as a standalone tool.
    """
    verbose(f"Running CESAR for {input_file}")
    x_param = math.ceil(memlim + 1) if memlim else math.ceil(memory_raw) + 1
    if memlim and memory_raw > memlim:
        eprint(
            f"Warning! Memory limit is {memlim} but {memory_raw} requested for this run."
        )

    # check whether the function was called correctly
    if input_file:
        cesar_cmd = [cesar_binary, input_file, '-x', str(x_param)] 
    elif input_data:
        cesar_cmd = [cesar_binary, '/dev/stdin', '-x', str(x_param)] 
    else:
        raise ValueError("run_cesar: both input_file and input_data are missing!")

    verbose(f"Calling CESAR command:\n{' '.join(cesar_cmd)}")

    p = subprocess.Popen(
        cesar_cmd, 
        shell=False,
        stdin=subprocess.PIPE,  # use PIPE for stdin
        stderr=subprocess.PIPE, 
        stdout=subprocess.PIPE
    )

    if input_file:
        # if a file was used to call CESAR -> just get the stdout and stderr
        b_stdout, b_stderr = p.communicate()
    elif input_data:
        # otherwise, feed the input_data using stdin
        b_stdout, b_stderr = p.communicate(input=input_data.encode())
    else:
        # b_stdout, b_stderr = None, None
        raise RuntimeError("Unreachable condition")

    rc = p.returncode
    if rc != 0:
        # CESAR job failed: die
        stderr = b_stderr.decode("utf-8")
        # os.remove(input_file) if is temp else None
        eprint("CESAR wrapper command crashed: {0}".format(" ".join(sys.argv)))
        die(f"CESAR failed run for {cesar_cmd}: {stderr}", 1)
        return None  # to suppress linter
    cesar_out = b_stdout.decode("utf-8")
    return cesar_out


def save(output, dest, t0_, loss_report=None):
    """Save raw CESAR output and quit."""
    f = open(dest, "w") if dest != "stdout" else sys.stdout
    f.write(output)
    if loss_report:
        f.write(loss_report)
    f.close() if dest != "stdout" else None
    verbose(f"Estimated: {dt.now() - t0_}")
    sys.exit(0)


def compute_percent_id(seq_1, seq_2):
    """Return % identity for two sequences."""
    verbose("Computing %ID for sequences:")
    verbose(seq_1)
    verbose(seq_2)
    assert len(seq_1) == len(seq_2)  # otherwise it is a bug
    matches = sum(
        [
            1
            for i in range(len(seq_1))
            if seq_1[i] == seq_2[i] and seq_1[i] != "N" and seq_2[i] != "N"
        ]
    )
    length = len(seq_1.replace("N", ""))  # we ignore N's
    pid = matches * 100 / length if length != 0 else 0
    verbose(f"PID is {pid}\n")
    return pid


def make_blosum_matrix():
    """Make python object with BLOSUM matrix."""
    # fill this matrix:
    # key will be amino_acid_1 -> amino_acid_2
    matrix = defaultdict(dict)
    num = 0
    alphabet = []
    f = open(BLOSUM_FILE, "r")
    for line in f:
        if line.startswith(" "):
            # the first line
            alphabet = line.split()
            continue
        # alphabet is given
        assert len(alphabet) > 0  # should be obtained before we reach this path
        line_data = line.split()
        del line_data[0]  # doesn't make any sense
        scores = [int(x) for x in line_data]
        row_char = alphabet[num]
        for subnum, score in enumerate(scores):
            col_char = alphabet[subnum]
            matrix[row_char][col_char] = score
        num += 1
    f.close()
    return matrix


def get_blosum_score(seq_1, seq_2, matrix):
    """Compute BLOSUM score."""
    score = 0
    for ch1, ch2 in zip(seq_1, seq_2):
        # a gap or not a gap
        if ch1 == "-" and ch2 == "-":
            continue
        elif ch1 == "-" and ch2 != "-":
            # insertion
            score += INS_PEN
            continue
        elif ch1 != "-" and ch2 == "-":
            score += DEL_PEN
        # this should never happen
        # if ch1 == "-" or ch2 == "-":
        #     score += GAP_SCORE
        #     continue
        else:
            score += matrix[ch1][ch2]
            continue
    return score


def translate_codons(codons_list):
    ### OBSOLETE
    """Translate codon table to ref and que AA sequences."""
    ref_AA_seq = []
    que_AA_seq = []
    for codon in codons_list:
        # extract codon sequences for reference and query
        ref_nuc_seq = codon["ref_codon"]
        que_nuc_seq = codon["que_codon"]
        # if there are any gaps -> remove them
        ref_only_let_nul = ref_nuc_seq.replace("-", "")
        # anyways it's possible that we cannot translate the codon, reasons are various
        # maybe there's N somewhere in the codon, or it's frame-shifted
        ref_only_let = (
            GENETIC_CODE.get(ref_only_let_nul, "X") if len(ref_only_let_nul) > 0 else "-"
        )
        # in case if N query codons correspond to a single reference codon
        # there is guarantee that ref_sequence has only one codon
        # for query -> it's unknown
        full_codons = len(que_nuc_seq) // 3  # get number of codons
        # fill extra codons (that are not in ref but are in query) with gaps:
        ref_AA = list("-" * (full_codons - 1) + ref_only_let)
        in_frame = len(que_nuc_seq) % 3 == 0
        if not in_frame:
            # in this case we don't know exact AA sequence in query
            # X -> frame-shifted codon
            que_AA = list("X" * full_codons)
        else:
            # possible N query codons: split query sequence into pieces of 3 characters
            que_codons = parts(que_nuc_seq, 3)
            que_AA = [GENETIC_CODE.get(c, "X") for c in que_codons]
        ref_AA_seq.extend(ref_AA)
        que_AA_seq.extend(que_AA)
    return ref_AA_seq, que_AA_seq


def compute_score(codon_data):
    """Compute BLOSUM score per exon."""
    blosum_matrix = make_blosum_matrix()
    # go exon_by_exon
    # extract all possible exon nums and sort them
    # exon_nums = sorted(set(x["q_exon_num"] for x in codon_data))
    prot_sequences = {}
    # get exons containing only split codons
    # split_codons = [c for c in codon_data if c["split_"] != 0]
    # split_exons_ = [c["q_exon_num"] - 1 for c in split_codons]
    # split_exons = set(split_exons_).difference(exon_nums)
    exon_pid, exon_blosum = {}, {}

    target_AAs_all = []
    exon_score = {}
    # check stops in split codons
    # exon_number: sequence -> one for reference (tar), another for query (que)
    exon_to_seq_tar = defaultdict(str)
    exon_to_seq_que = defaultdict(str)
    for codon in codon_data:
        exon_num = codon["q_exon_num"]
        split = codon["split_"]  # if split codon - it's a bit tricky
        if split == 0:
            # not a split codon: the codon lies in some exon entirely
            exon_to_seq_tar[exon_num] += codon["ref_codon"]
            exon_to_seq_que[exon_num] += codon["que_codon"]
        else:
            # codon is split between 2 exons (I hope not 3)
            tar_before = codon["ref_codon"][:split]
            que_before = codon["que_codon"][:split]
            tar_after = codon["ref_codon"][split:]
            que_after = codon["que_codon"][split:]
            # split this codon and add different parts to different exons
            exon_to_seq_tar[exon_num - 1] += tar_before
            exon_to_seq_tar[exon_num] += tar_after
            exon_to_seq_que[exon_num - 1] += que_before
            exon_to_seq_que[exon_num] += que_after

    for exon_num in exon_to_seq_tar.keys():
        # exon sequences extracted, now compute scores per exon
        tar_exon = exon_to_seq_tar[exon_num]
        que_exon = exon_to_seq_que[exon_num]
        perc_id = compute_percent_id(tar_exon, que_exon)
        exon_pid[exon_num] = perc_id  # pid: nucleotide %ID

        # for BLOSUM we need AA sequence
        # for AA sequence we need complete codons
        codons_in_exon = [
            x for x in codon_data if x["q_exon_num"] == exon_num and x["split_"] == 0
        ]
        codons_in_exon_w_splice = [x for x in codon_data if x["q_exon_num"] == exon_num]

        if len(codons_in_exon) == 0:
            # possible that there are no full codons in the exon
            # like here: AT ----- G CGA AAA
            # Exon 1 consists on 2 letters: AT
            exon_blosum[exon_num] = 50  # 50 if default value then
            continue

        ref_codons = [x["ref_codon"].replace("-", "") for x in codons_in_exon]
        # que_codons = [x["que_codon"].replace("-", "") for x in codons_in_exon]
        __ref_codons = [x["ref_codon"] for x in codons_in_exon]
        __que_codons = [x["que_codon"] for x in codons_in_exon]

        if len(ref_codons) == 0:  # exon containing only split codons
            exon_score[exon_num] = 50  # default value in this case
            continue

        ref_aa, que_aa = translate_codons(codons_in_exon)
        ref_aa_s_, que_aa_s_ = translate_codons(codons_in_exon_w_splice)
        # returns lists
        ref_aa_s = "".join(ref_aa_s_)
        que_aa_s = "".join(que_aa_s_)
        prot_sequences[exon_num] = {"ref": ref_aa_s, "que": que_aa_s}

        # different rule for query, it's possible that value is absent
        target_AAs_all.extend(ref_aa)  # keep the entire target AA seq for sanity checks
        # get absolute and maximal BLOSUM scores
        verbose("Computing BLOSUM score for: ")
        verbose(str(ref_aa))
        verbose(str(que_aa))
        # blosum score: the blosum score of ref VS query
        # max: reference VS reference
        # it is not ACTUALLY max, but in 99% cases it is
        # in other -> doesn't matter
        blosum_score = get_blosum_score(ref_aa, que_aa, blosum_matrix)
        max_blosum_score = get_blosum_score(ref_aa, ref_aa, blosum_matrix)
        if max_blosum_score > 0:
            rel_blosum_score = (blosum_score / max_blosum_score) * 100
        else:
            # for incredibly short exons everything might happen
            rel_blosum_score = 0
            eprint(f"Warning! Max blosum score < 0!")
        rel_blosum_score = rel_blosum_score if rel_blosum_score >= 0 else 0
        verbose(f"Max BLOSUM score is: {max_blosum_score}")
        verbose(f"Absolute BLOSUM score is: {blosum_score}")
        verbose(f"Relative BLOSUM score is: {rel_blosum_score}")
        exon_blosum[exon_num] = rel_blosum_score
    return exon_pid, exon_blosum, prot_sequences


def split_indexes(indexes):
    """Split indexes list like 1 2 5 in 1 2 and 5."""
    verbose(f"Analyze indexes:\n{indexes}")
    if len(indexes) == 0:
        return [], []
    left, right = [
        indexes[0],
    ], []
    left_now = True
    for i in range(1, len(indexes)):
        prev = indexes[i - 1]
        curr = indexes[i]
        if curr > prev + 1 and left_now:
            left_now = False
        if left_now:
            left.append(curr)
        else:
            right.append(curr)
    return left, right


def extract_query_seq(codons_data):
    """Parse codons data."""
    exon_to_que_seq = defaultdict(str)
    exon_to_ref_seq = defaultdict(str)
    exon_to_indexes = defaultdict(list)

    for codon_data in codons_data:
        exon_num = codon_data["q_exon_num"]
        split = codon_data["split_"]
        # to guarantee that a(n) < a(n + 1)
        indexes = sorted(codon_data["que_coords"])
        if split == 0:
            # regular codon, nothing special
            exon_to_que_seq[exon_num] += codon_data["que_codon"]  # .replace("-", "")
            exon_to_ref_seq[exon_num] += codon_data["ref_codon"]
            exon_to_indexes[exon_num].extend(indexes)
        else:
            # not a regular, split codon
            to_prev_let = codon_data["que_codon"][:split]  # .replace("-", "")
            to_curr_let = codon_data["que_codon"][split:]  # .replace("-", "")
            to_prev_let_r = codon_data["ref_codon"][:split]  # .replace("-", "")
            to_curr_let_r = codon_data["ref_codon"][split:]  # .replace("-", "")
            if len(to_prev_let) == 0 and len(to_curr_let) == 0:
                # since I removed replace -> this branch is impossible
                continue
            exon_to_que_seq[exon_num - 1] += to_prev_let
            exon_to_que_seq[exon_num] += to_curr_let

            exon_to_ref_seq[exon_num - 1] += to_prev_let_r
            exon_to_ref_seq[exon_num] += to_curr_let_r

            left_ind, right_ind = split_indexes(indexes)
            if len(left_ind) > 0 and len(right_ind) > 0:
                exon_to_indexes[exon_num - 1].extend(left_ind)
                exon_to_indexes[exon_num].extend(right_ind)
                continue
            # if we are here, we couldn't split this codon correctly
            if len(to_curr_let.replace("-", "")) > 0:
                exon_to_indexes[exon_num].extend(left_ind)
            else:
                exon_to_indexes[exon_num - 1].extend(left_ind)

    empty_exons = set()
    exon_nums = list(exon_to_que_seq.keys())
    for ex_num in exon_nums:
        # if exon in query contains only gaps, we need to keep it
        # to avoid ghost-exons in the future
        seq = exon_to_que_seq[ex_num].replace("-", "")
        if len(seq) > 0:
            continue
        empty_exons.add(ex_num)
        exon_to_indexes[ex_num] = [0]
        exon_to_que_seq[ex_num] = "-"
    return exon_to_que_seq, exon_to_ref_seq, empty_exons, exon_to_indexes


def get_exon_num_corr(codons_data):
    """Get correspondence between exon numbers in Q and T."""
    ans = defaultdict(set)
    for codon in codons_data:
        ans[codon["q_exon_num"]].add(codon["t_exon_num"])
    if not ans[0]:
        ans[0] = {0}
    return ans


def check_codons_for_aa_sat(codons_data):
    exon_codon_marks = defaultdict(list)
    for codon in codons_data:
        t_exon_num = codon["t_exon_num"]
        ref_codon_seq = codon["ref_codon"]
        que_codon_seq = codon["que_codon"]
        no_gaps = ("-" not in ref_codon_seq) and ("-" not in que_codon_seq)
        no_ins = len(ref_codon_seq) == len(que_codon_seq) == 3
        if no_gaps and no_ins:
            exon_codon_marks[t_exon_num].append(True)
        else:
            exon_codon_marks[t_exon_num].append(False)
    exon_verdict = {}
    for exon_num, codon_marks in exon_codon_marks.items():
        if all(x is True for x in codon_marks):
            exon_verdict[exon_num] = True
        else:
            exon_verdict[exon_num] = False
    return exon_verdict


def check_codon(codon):
    """Mask codon if FS or stop."""
    if codon in Constants.STOP_CODONS:
        # mask stop
        return Constants.XXX_CODON
    elif codon == Constants.GAP_CODON:
        # just a gap
        return codon
    elif codon.count("-") > 0:
        # FS
        return Constants.XXX_CODON
    else:  # normal codon
        return codon


def extract_codon_data(codon_table, excl_exons=None):
    """Get codon alignment."""
    t_codons = []
    q_codons = []
    _excl_exons = excl_exons if excl_exons is not None else set()
    prev_exon_was_del = None

    for codon in codon_table:
        ref_codon = codon["ref_codon"]
        que_codon = codon["que_codon"]
        # needed if we filter D/M exons out
        t_exon_num = codon["t_exon_num"]

        # check whether we need to delete the codon or not
        this_exon_to_del = t_exon_num in _excl_exons
        prev_exon_to_del_this_codon_split = prev_exon_was_del and codon["split_"] > 0

        if this_exon_to_del or prev_exon_to_del_this_codon_split:
            # this codon to be deleted in query
            if len(ref_codon) == 3:
                ref_codon = check_codon(ref_codon)
                t_codons.append(ref_codon)
                q_codons.append(Constants.GAP_CODON)
            elif len(que_codon) % 3 == 0:
                # Frame-preserving ins
                ref_subcodons = [check_codon(x) for x in parts(ref_codon, 3)]
                que_subcodons = [Constants.GAP_CODON for _ in range(len(ref_subcodons))]
                t_codons.extend(ref_subcodons)
                q_codons.extend(que_subcodons)
            elif len(ref_codon) < 3:
                t_codons.append(Constants.XXX_CODON)
                q_codons.append(Constants.GAP_CODON)
            elif ref_codon == "-" or ref_codon == "--":
                # special case to be captured Apr 2023
                t_codons.append(Constants.GAP_CODON)
                q_codons.append(Constants.XXX_CODON)
            elif len(ref_codon) < 3:
                t_codons.append(Constants.XXX_CODON)
                q_codons.append(Constants.XXX_CODON)
            else:
                # something strange
                ref_int = check_codon(ref_codon[-3:])
                t_codons.append(ref_int)
                q_codons.append(Constants.GAP_CODON)

            if this_exon_to_del:
                prev_exon_was_del = True
            else:
                # meaning, prev_exon_to_del_this_codon_split is True
                prev_exon_was_del = False
            continue

        if t_exon_num in _excl_exons:
            prev_exon_was_del = True
            if len(ref_codon) == 3:
                # TODO: check what is that
                ref_codon = check_codon(ref_codon)
            continue

        elif prev_exon_was_del and codon["split_"] > 0:
            # split codon from prev exon
            prev_exon_was_del = False
            continue
        # this codon is added for sure
        prev_exon_was_del = False

        # ordinary branch
        if len(ref_codon) == len(que_codon) == 3:
            # the simplest case
            ref_codon = check_codon(ref_codon)
            que_codon = check_codon(que_codon)
            t_codons.append(ref_codon)
            q_codons.append(que_codon)
            continue
        elif len(que_codon) % 3 == 0:
            # frame-preserving insertion
            ref_subcodons = [check_codon(x) for x in parts(ref_codon, 3)]
            que_subcodons = [check_codon(x) for x in parts(que_codon, 3)]
            t_codons.extend(ref_subcodons)
            q_codons.extend(que_subcodons)
            continue
        elif len(ref_codon) < 3:
            t_codons.append(Constants.XXX_CODON)
            q_codons.append(Constants.XXX_CODON)
        else:
            ref_int = check_codon(ref_codon[-3:])
            que_int = check_codon(que_codon[-3:])
            fs_ref = Constants.GAP_CODON
            fs_que = Constants.XXX_CODON
            t_codons.append(fs_ref)
            q_codons.append(fs_que)
            t_codons.append(ref_int)
            q_codons.append(que_int)
    return {"ref": t_codons, "que": q_codons}


def process_cesar_out(cesar_raw_out, query_loci, inverts):
    """Extract data from raw cesar output."""
    raw_data = parts(cesar_raw_out.split("\n"), n=4)[:-1]
    exon_queries, exon_refs, percIDs, blosums, query_coords, prot_seqs, codon_seqs = (
        {},
        {},
        {},
        {},
        {},
        {},
        {},
    )
    exon_num_corr = {}  # query exon num -> target exon nums
    aa_sat_seq = {}  # exon num -> no dels
    chain_id_to_codon_table = {}  # chain id -> raw codon table

    for part in raw_data:
        # parse output
        target_raw = part[1]
        chain_id = int(part[2][1:])
        query_raw = part[3]
        verbose(f"Processing query ID {chain_id}")

        codons_data = parse_cesar_out(target_raw, query_raw, v=VERBOSE)
        chain_id_to_codon_table[chain_id] = codons_data
        # extract protein sequences also here, just to do it in one place
        part_pIDs, part_blosums, prot_seqs_part = compute_score(codons_data)
        prot_seqs[chain_id] = prot_seqs_part
        # convert codon list to {"ref": [t_codons], "que": [q_codons]}
        codon_seqs_part = extract_codon_data(codons_data)
        codon_seqs[chain_id] = codon_seqs_part

        (
            exon_query_seqs,
            exon_ref_seqs,
            empty_q_exons,
            exon_query_inds,
        ) = extract_query_seq(codons_data)
        exon_num_corr[chain_id] = get_exon_num_corr(codons_data)
        aa_sat_part = check_codons_for_aa_sat(codons_data)
        aa_sat_seq[chain_id] = aa_sat_part
        abs_coords = {}

        # get abs coordinates
        directed = inverts[chain_id]
        locus = query_loci[chain_id]
        q_chrom, region = locus.split(":")
        abs_query_start = int(region.split("-")[0])
        abs_query_end = int(region.split("-")[1])
        for exon_num, indexes in exon_query_inds.items():
            if exon_num in empty_q_exons:
                abs_coords[exon_num] = Constants.UNDEF_REGION
                continue
            elif len(indexes) == 0:
                # 0 indexes: definitely cannot find the exon
                abs_coords[exon_num] = Constants.UNDEF_REGION
                continue
            rel_start, rel_len = indexes[0], len(indexes)
            exon_abs_start = (
                abs_query_start + rel_start if directed else abs_query_end - rel_start
            )
            exon_abs_end = (
                abs_query_start + rel_start + rel_len
                if directed
                else abs_query_end - rel_start - rel_len
            )
            exon_grange = f"{q_chrom}:{exon_abs_start}-{exon_abs_end}"
            abs_coords[exon_num] = exon_grange
        # add to the global dicts
        exon_queries[chain_id] = exon_query_seqs
        exon_refs[chain_id] = exon_ref_seqs
        percIDs[chain_id] = part_pIDs
        blosums[chain_id] = part_blosums
        query_coords[chain_id] = abs_coords

    ret = (
        exon_queries,
        exon_refs,
        percIDs,
        blosums,
        query_coords,
        exon_num_corr,
        prot_seqs,
        codon_seqs,
        aa_sat_seq,
        chain_id_to_codon_table,
    )
    return ret


def __check_fragm_coords_intersect_in_q(query_coords):
    """Check whether coordinates of predicted exons intersect in the query."""
    exon_ranges_w_chrom = list(query_coords.values())
    exon_ranges = [x.split(":")[1].split("-") for x in exon_ranges_w_chrom]
    exon_starts = [int(x[0]) for x in exon_ranges]
    exon_ends = [int(x[1]) for x in exon_ranges]
    exon_intervals = list(zip(exon_starts, exon_ends))
    exon_intervals_sorted = sorted(exon_intervals, key=lambda x: x[0])
    intersecting_intervals = _check_seq_of_intervals_intersect(exon_intervals_sorted)
    return intersecting_intervals


def process_cesar_out__fragments(cesar_raw_out, fragm_data, query_loci, inverts):
    """Process CESAR output for assembled from fragments gene."""
    exon_queries, exon_refs, percIDs, blosums, query_coords, = (
        {},
        {},
        {},
        {},
        {},
    )
    prot_seqs, codon_seqs, exon_num_corr, aa_sat_seq = {}, {}, {}, {}
    cesar_raw_lines = cesar_raw_out.split("\n")
    target_seq_raw = cesar_raw_lines[1]
    query_seq_raw = cesar_raw_lines[3]
    chain_id_to_codon_table = {}  # chain id -> raw codon table

    codons_data = parse_cesar_out(target_seq_raw, query_seq_raw, v=VERBOSE)
    chain_id_to_codon_table[Constants.FRAGMENT_CHAIN_ID] = codons_data
    # extract protein sequences also here, just to do it in one place
    part_pIDs, part_blosums, prot_seqs_part = compute_score(codons_data)
    prot_seqs[Constants.FRAGMENT_CHAIN_ID] = prot_seqs_part
    codon_seqs_part = extract_codon_data(codons_data)
    codon_seqs[Constants.FRAGMENT_CHAIN_ID] = codon_seqs_part
    exon_query_seqs, exon_ref_seqs, empty_q_exons, exon_query_inds = extract_query_seq(
        codons_data
    )
    exon_num_corr[Constants.FRAGMENT_CHAIN_ID] = get_exon_num_corr(codons_data)
    aa_sat_part = check_codons_for_aa_sat(codons_data)
    aa_sat_seq[Constants.FRAGMENT_CHAIN_ID] = aa_sat_part
    abs_coords = {}

    # extract coordinates of different exons!
    q_lens = [x[7] for x in fragm_data]
    verbose(f"QLens: {q_lens}")
    # print(q_lens)
    q_limits = [
        0,
    ]
    for q_len in q_lens:
        new_el = q_len + q_limits[-1]
        q_limits.append(new_el)
    verbose(f"q limits: {q_limits}")

    for exon_num, indexes in exon_query_inds.items():
        verbose(exon_num)
        if exon_num in empty_q_exons:
            verbose("In empty Q exons")
            abs_coords[exon_num] = Constants.UNDEF_REGION
            continue
        elif len(indexes) == 0:
            verbose("Exon not mapped")
            abs_coords[exon_num] = Constants.UNDEF_REGION
            continue
        rel_start_not_corr, rel_len = indexes[0], len(indexes)
        borders_below_start = [x for x in q_limits if x <= rel_start_not_corr]
        num_borders_before = len(borders_below_start)
        chain_index_in_list = num_borders_before - 1
        closest_border = max(borders_below_start)
        resp_chain = fragm_data[chain_index_in_list]
        chain_id = resp_chain[0]
        q_size_ = resp_chain[8]
        directed = inverts[chain_id]
        verbose(f"Direction: {directed}")
        locus = query_loci[chain_id]
        rel_start = rel_start_not_corr - closest_border
        verbose(f"Rel starts not corr: {rel_start_not_corr}, corr: {rel_start}")
        verbose(f"closest: {closest_border}")

        q_chrom, region = locus.split(":")
        verbose(f"Chain_ID: {chain_id} Q chrom: {q_chrom}")
        abs_query_start = int(region.split("-")[0])
        abs_query_end = int(region.split("-")[1])

        exon_abs_start = (
            abs_query_start + rel_start if directed else abs_query_end - rel_start
        )
        exon_abs_end = (
            abs_query_start + rel_start + rel_len
            if directed
            else abs_query_end - rel_start - rel_len
        )
        # check that everything in borders
        left_side_ok = exon_abs_start > 0 and exon_abs_end > 0
        right_size_ok = exon_abs_start <= q_size_ and exon_abs_end <= q_size_
        if not left_side_ok or not right_size_ok:
            # out of borders: this is possible
            # let's just skip this exon
            exon_grange = Constants.UNDEF_REGION
        else:
            # this is completely OK
            exon_grange = f"{q_chrom}:{exon_abs_start}-{exon_abs_end}"
        abs_coords[exon_num] = exon_grange

    # check that abs coords do not intersect
    coords_intersect = __check_fragm_coords_intersect_in_q(abs_coords)
    if coords_intersect:  # error, cannot go further
        print("Error! Cannot stitch fragments properly, exon intervals intersect after merge")
        print(f"Intersecting intervals are: {coords_intersect}")
        print("Abort")
        sys.exit(ERR_CODE_FRAGM_ERR)
    # add to the global dicts
    exon_queries[Constants.FRAGMENT_CHAIN_ID] = exon_query_seqs
    exon_refs[Constants.FRAGMENT_CHAIN_ID] = exon_ref_seqs
    percIDs[Constants.FRAGMENT_CHAIN_ID] = part_pIDs
    blosums[Constants.FRAGMENT_CHAIN_ID] = part_blosums
    query_coords[Constants.FRAGMENT_CHAIN_ID] = abs_coords
    # return the same values as process_cesar_out does
    ret = (
        exon_queries,
        exon_refs,
        percIDs,
        blosums,
        query_coords,
        exon_num_corr,
        prot_seqs,
        codon_seqs,
        aa_sat_seq,
        chain_id_to_codon_table,
    )
    return ret


def merge_regions(reg_list):
    """Merge regions list."""
    all_numbers = flatten(reg_list)
    return min(all_numbers), max(all_numbers)


def check_region(exon_exp_reg, act_reg, t_nums):
    """Check if exon found in the expected region."""
    q_regions = [exon_exp_reg.get(n) for n in t_nums if exon_exp_reg.get(n)]
    if not q_regions:
        # target exons are not covered with chain
        # expected region does not exist actually
        return "EXCL", (0, 0)
    q_region = merge_regions(q_regions)
    they_intersect = intersect_ranges(act_reg, q_region) > 0
    if they_intersect:
        return "INC", q_region
    else:
        return "EXCL", q_region


def arrange_output(
    gene,
    exon_seqs,
    query_exon_sequences,
    p_ids,
    p_bl,
    all_query_coords,
    chain_exon_gap=None,
    chain_exon_class=None,
    chain_exon_exp_reg=None,
    ch_q_to_t_num=None,
    missed=None,
    is_paral=False,
):
    """Arrange final fasta file with additional data."""
    header_template = (
        ">{0} | {1} | {2} | {3} | {4:.2f} | {5:.2f} | {6} | {7} "
        "| exp:{8}-{9} | {10} | {11} | query_exon\n"
    )
    chains = list(query_exon_sequences.keys())  # it is output for these chains
    extra_fields = True if chain_exon_gap else None
    output = ""
    chain_to_excl = {}

    for chain in chains:
        # output chain-by-chain
        query_seqs = query_exon_sequences[chain]
        reference_seqs = exon_seqs[chain]
        query_pids = p_ids[chain]
        query_blosums = p_bl[chain]
        q_to_t_num = ch_q_to_t_num[chain]
        query_coords = all_query_coords[chain]
        exons_gap_data = chain_exon_gap[chain] if chain_exon_gap else {}
        exon_class_data = chain_exon_class[chain] if chain_exon_class else {}
        exon_exp_reg = chain_exon_exp_reg[chain] if chain_exon_exp_reg else None
        exon_nums = query_seqs.keys()
        exons_missed = missed[chain] if missed else []
        exon_to_inc = {}

        # and exon-by-exon
        for exon_num in exon_nums:
            # collect data for gene loss detector!
            reference_seq = reference_seqs.get(exon_num, "N")
            query_seq = query_seqs.get(exon_num, "N")
            query_pid = query_pids.get(exon_num, 0.0)
            query_blo = query_blosums.get(exon_num, 0.0)
            coord = query_coords.get(exon_num, Constants.UNDEF_REGION)

            # extract marks for the exon
            if exon_num in exons_missed or not extra_fields:
                is_gap = "N/A"
                exon_class = "N/A"
                exp_reg, reg_data = ("N/A", "N/A"), "N/A"
            else:
                is_gap = "GAP" if exons_gap_data.get(exon_num, None) else "OK"
                exon_class = exon_class_data.get(exon_num, "N/A")
                # exp_reg = exon_exp_reg[exon_num]
                act_reg = [int(x) for x in coord.split(":")[1].split("-")]
                # reg_data = "INC" if intersect_ranges(act_reg, exp_reg) > 0 else "EXCL"
                t_nums = q_to_t_num[exon_num]
                reg_data, exp_reg = check_region(exon_exp_reg, act_reg, t_nums)
                for n_ in t_nums:
                    exon_to_inc[n_] = True if reg_data == "INC" else False

            # add reference seq for this exon
            ref_datum = (
                f">{gene} | {exon_num} | {chain} | reference_exon\n{reference_seq}\n"
            )
            output += ref_datum
            # and then merge all this stuff altogether
            header = header_template.format(
                gene,
                exon_num,
                chain,
                coord,
                query_pid,
                query_blo,
                is_gap,
                exon_class,
                exp_reg[0],
                exp_reg[1],
                reg_data,
                is_paral,
            )
            output += header
            output += "{0}\n".format(query_seq)
            chain_to_excl[chain] = exon_to_inc
    return output, chain_to_excl


def aa_eq_len_check(ref_sequences, query_exon_sequences):
    """Check that exon lens satisfy A+ criteria."""
    chain_eq_len_sat = {}
    for chain, q_seqs in query_exon_sequences.items():
        chain_eq_len_sat[chain] = {k: False for k in q_seqs.keys()}
        for exon_num, q_seq in q_seqs.items():
            r_seq = ref_sequences.get(exon_num, "")
            if len(r_seq) == len(q_seq):
                chain_eq_len_sat[chain][exon_num] = True
    return chain_eq_len_sat


def invert_complement(seq):
    """Make inverted-complement sequence."""
    reverse = seq[::-1]
    reverse_complement = "".join(
        [COMPLEMENT_BASE.get(c) if COMPLEMENT_BASE.get(c) else "N" for c in reverse]
    )
    return reverse_complement


def save_prot(prot_seq, prot_out):
    """Save protein sequences."""
    f = open(prot_out, "w") if prot_out != "stdout" else sys.stdout
    for proj_name, seq_collection in prot_seq.items():
        exons = sorted(seq_collection.keys())

        ref_aa_seq = seq_collection["ref"]
        que_aa_seq = seq_collection["que"]
        # do not allow 0-length sequences, otherwise it will ruin everything
        que_aa_seq = que_aa_seq if len(que_aa_seq) > 0 else "X"
        f.write(f">{proj_name} | PROT | REFERENCE\n")
        f.write(f"{ref_aa_seq}\n")
        f.write(f">{proj_name} | PROT | QUERY\n")
        f.write(f"{que_aa_seq}\n")
    f.close() if prot_out != "stdout" else None


def save_codons(gene_name, cds_data_all, output):
    """Save codon alignments."""
    if output is None:
        return
    f = open(output, "w") if output != "stdout" else sys.stdout
    for key, cds_data in cds_data_all.items():
        proj_name = f"{gene_name}.{key}"
        ref_codons = cds_data["ref"]
        que_codons = cds_data["que"]
        ref_seq = " ".join(ref_codons)
        que_seq = " ".join(que_codons)
        f.write(f">{proj_name} | CODON | REFERENCE\n")
        f.write(f"{ref_seq}\n")
        f.write(f">{proj_name} | CODON | QUERY\n")
        f.write(f"{que_seq}\n")
    f.close() if output != "stdout" else None


def get_a_plus(chain_exon_class, aa_cesar_sat, aa_block_sat_chain, aa_ex_len):
    """Mark A+ exons as A+."""
    chains = chain_exon_class.keys()
    to_mark_a_p = []
    for chain in chains:
        c_exon_classes = chain_exon_class[chain]
        c_cesar_sat = aa_cesar_sat[chain]
        c_sat_chain = aa_block_sat_chain[chain]
        c_sat_len = aa_ex_len[chain]
        exons = c_exon_classes.keys()
        for exon in exons:
            e_class = c_exon_classes.get(exon)
            e_cesar_sat = c_cesar_sat.get(exon)
            e_chain_sat = c_sat_chain.get(exon)
            e_len_sat = c_sat_len.get(exon)
            bools = (e_cesar_sat, e_chain_sat, e_len_sat)
            if e_class == "A" and all(x is True for x in bools):
                mark = (chain, exon)
                to_mark_a_p.append(mark)
            else:
                pass
    for mark in to_mark_a_p:
        chain = mark[0]
        exon = mark[1]
        chain_exon_class[chain][exon] = "A+"
    return chain_exon_class


def analyse_ref_ss(s_sites):
    """Create list of ss with ignored SSM."""
    del s_sites[0]
    del s_sites[-1]
    intron_elems = parts(s_sites, 2)
    masked_ss = set()
    for num, intron in enumerate(intron_elems, 1):
        don_ = intron[0].lower()
        acc_ = intron[1].lower()
        if acc_ in Constants.ACCEPTOR_SITE and don_ in Constants.DONOR_SITE:
            # too usual and boring
            continue
        verbose(f"Intron num {num} has don: {don_} acc: {acc_}")
        mss_1 = (num, 1)
        mss_2 = (num + 1, 0)
        masked_ss.add(mss_1)
        masked_ss.add(mss_2)
    return masked_ss


def append_u12(u12_base, gene, ref_ss_data):
    """Add U12 data to ref ss data."""
    if u12_base is None:
        return
    f = open(u12_base, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        if trans != gene:
            continue
        exon_num = int(line_data[1])
        side = 0 if line_data[2] == "A" else 1
        u12_site = (exon_num, side)
        ref_ss_data.add(u12_site)
    f.close()


def merge_dicts(dicts_list):
    """Merge a list of dicts."""
    ret = {}
    for d in dicts_list:
        ret.update(d)
    return ret


def intersect_lists(lists):
    """Perform intersection operation on several lists."""
    # noinspection PyTypeChecker
    intersection = reduce(and_, [set(x) for x in lists])
    return list(intersection)


def get_relative_coordinates(exon_exp_region, search_locus, directed, max_len=5000):
    """Get relative coordinates of expected regions.

    Formatting for LASTZ-CESAR optimizer.
    """
    include_regions = list(exon_exp_region.values())
    _, start_end_str = search_locus.split(":")
    start_str, end_str = start_end_str.split("-")
    start, end = int(start_str), int(end_str)
    if directed:
        # just subtract absolute start coord
        include_regions_rel = [(x[0] - start, x[1] - start) for x in include_regions]
    else:
        # a bit more complicated
        include_regions_rel = [(end - x[0], end - x[1]) for x in include_regions]
    include_regions_len_filt = [
        x for x in include_regions_rel if abs(x[0] - x[1]) < max_len
    ]
    return include_regions_len_filt


def _check_seq_of_intervals_intersect(intervals):
    intervals_num = len(intervals)
    for i in range(intervals_num - 1):
        # (start, end)
        curr_one = intervals[i]
        next_one = intervals[i + 1]
        # sorted by beginning
        # if start of the next < end of the curr
        # -> they intersect
        if next_one[0] < curr_one[1]:
            return curr_one[0], curr_one[1]
    return None  # nothing suspicious found


def parse_precomp_orth_loci(transcript_name, path):
    """Read precomputed orthologous loci from a file."""
    ret_1 = {}  # chain to search locus
    ret_2 = {}  # chain to subchain locus

    f = open(path, "r")
    # sample file line:
    # #ORTHLOC	ENST00000262455	1169	JH567521:462931-522892	JH567521:462931-522892
    # suffix - transcript - chain - search locus - subch locus
    for line in f:
        ld = line.rstrip().split("\t")
        if ld[1] != transcript_name:
            continue
        chain_id = int(ld[2])
        ret_1[chain_id] = ld[3]
        ret_2[chain_id] = ld[4]
    f.close()

    return ret_1, ret_2


def redo_codon_sequences(codon_tables, del_mis_exons):
    """Rebuild codon alignments: now excluding deleted and missing exons."""
    codon_ret = {}
    for chain_id, codon_table in codon_tables.items():
        excl_exons = del_mis_exons.get(str(chain_id), set())
        codon_seqs_upd = extract_codon_data(codon_table, excl_exons=excl_exons)
        codon_ret[chain_id] = codon_seqs_upd
    return codon_ret


def extract_prot_sequences_from_codon(gene_name, codon_s):
    """Extract protein sequences from codon"""
    ret = {}
    for chain_id, cds_data in codon_s.items():
        proj_name = f"{gene_name}.{chain_id}"
        ref_codons = cds_data["ref"]
        que_codons = cds_data["que"]
        ref_aa_seq = "".join([GENETIC_CODE.get(x, "X") for x in ref_codons])
        que_aa_seq = "".join([GENETIC_CODE.get(x, "X") for x in que_codons])
        ret[proj_name] = {"ref": ref_aa_seq, "que": que_aa_seq}
    return ret


def realign_exons(args):
    """Entry point."""
    memlim = float(args["memlim"]) if args["memlim"] != "Auto" else None
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash
    # read gene-related data
    bed_data = read_bed(
        args["gene"], args["bdb_bed_file"]
    )  # extract gene data from bed file
    # bed_exons_num = len(bed_data["blocks"])
    # parse gene bed-track: get exon coordinates, sequences and splice sites
    exon_coordinates, exon_sequences, s_sites, exon_flanks = get_exons(
        bed_data, args["tDB"]
    )
    # read chain IDs list:
    chains = (
        [int(x) for x in args["chains"].split(",") if x != ""]
        if args["chains"] != "region"
        else []
    )
    # get path to the chain file (most likely just arg)
    chain_file = (
        find_chain_file(args["ref"], args["bdb_chain_file"])
        if args["chains"] != "region"
        else args["bdb_chain_file"]
    )

    # check if there are stop codons in reference -> we either mask them or halt execution
    exon_sequences, sec_codons = check_ref_exons(exon_sequences, args["mask_stops"])
    # CESAR require some formatting of the reference exon sequences:
    prepared_exons = prepare_exons_for_cesar(exon_sequences)

    # read chain-related data
    query_sequences, query_loci, inverts = {}, {}, {}
    gene_range = "{0}:{1}-{2}".format(
        bed_data["chrom"], bed_data["chromStart"], bed_data["chromEnd"]
    )
    chain_exon_gap, chain_exon_class, chain_exon_exp_reg, chain_missed = {}, {}, {}, {}
    aa_block_sat_chain = {}  # one of dicts to mark exceptionally good exon predictions
    fragments_data = []  # required for fragmented genomes
    chains_in_input = True

    # no need to extract chain blocks and compute where is the orthologous locus
    # it's already precomputed
    chain_to_precomp_search_loci = {}
    chain_to_precomp_subch_loci = {}
    if args["precomputed_orth_loci"]:
        plc_ = parse_precomp_orth_loci(args["gene"], args["precomputed_orth_loci"])
        chain_to_precomp_search_loci = plc_[0]
        chain_to_precomp_subch_loci = plc_[1]

    # if chains and args["fragments"]:
    # the normal branch: call CESAR vs 1+ query sequences
    for chain_id in chains:  # in region more this part is skipped
        verbose(f"\nLoading chain {chain_id}")  # only one place where I need chain data
        # extract chain and coordinates of locus; extract sequence from query genome
        chain_str = get_chain(chain_file, chain_id)
        verbose(f"Chain {chain_id} extracted")
        chain_header = chain_str.split("\n")[0].split()
        # most likely we need only the chain part that intersects the gene
        # and skip the rest:
        verbose("Cutting the chain...")
        
        if args["precomputed_orth_loci"]:
            search_locus = chain_to_precomp_search_loci[chain_id]
            subch_locus = chain_to_precomp_subch_loci[chain_id]
            # TODO: chain_data contains non-necessary information
            # it's already present in the chain_header variable
            # can be optimised (historical reasons)
            _t_strand = True if chain_header[2] == "+" else False
            _q_strand = True if chain_header[7] == "+" else False
            _t_size = int(chain_header[3])
            _q_size = int(chain_header[8])
            chain_data = (_t_strand, _t_size, _q_strand, _q_size),
        else:
            search_locus, subch_locus, chain_data = chain_cut(
                chain_str, gene_range, args["gene_flank"], args["extra_flank"]
            )

        # chain data: t_strand, t_size, q_strand, q_size
        chain_query_strand = chain_data[2]
        chain_query_size = chain_data[3]
        verbose("Chain cut is done.")

        # this call of make_query_seq is actually for extracting
        # query sequence for CESAR:
        verbose("Extracting query sequence...")
        query_seq, directed = make_query_seq(
            chain_id, search_locus, args["qDB"], chain_query_strand, bed_data["strand"]
        )
        verbose("Query sequence extracted")
        q_seq_len = len(query_seq)
        # this is extended query seq (larger locus) for assembly gaps search only!
        # We do not call CESAR for this _query_seq_ext sequence
        _query_seq_ext, directed = make_query_seq(
            chain_id, subch_locus, args["qDB"], chain_query_strand, bed_data["strand"]
        )

        if args["ic"]:  # invert complement required for some reason
            query_seq = invert_complement(query_seq)
            _query_seq_ext = invert_complement(_query_seq_ext)
        if len(query_seq) > args["query_len_limit"] > 0:
            # query length limit exceeded:
            verbose(f"Skipping chain {chain_id} - too long")
            continue

        # extract gaps and block coordinates
        gap_coordinates = find_gaps(
            _query_seq_ext, subch_locus, args["gap_size"], directed
        )
        # blocks are [target_start, target_end, query_start, query_end]
        # TODO: can be optimised here
        # and also can be written to log - if same chain & bed - same output
        subchain_blocks_raw = extract_subchain(chain_str, subch_locus)
        # swap blocks in correct orientation and fill interblock ranges
        subchain_blocks = orient_blocks(subchain_blocks_raw, chain_data)
        # intersect exon: chain blocks and chain blocks: gaps, get exons not covered by chain
        block_intersection_out = intersect_exons_blocks_gaps(
            exon_coordinates,
            subchain_blocks,
            gap_coordinates,
            args["exon_flank"],
            args["uhq_flank"],
        )
        # parse block_intersection_out -> there are many different data:
        exon_blocks = block_intersection_out[0]
        # flanked_exon_blocks = block_intersection_out[1]
        blocks_gaps = block_intersection_out[2]
        missing_exons = block_intersection_out[3]
        # exon_flank_coordinates = block_intersection_out[4]
        margin_cases = block_intersection_out[5]
        aa_block_sat = block_intersection_out[6]
        verbose(f"AA sat: {aa_block_sat}")

        # classify exons, get expected regions
        exon_class, exon_exp_region = classify_predict_exons(
            exon_blocks, subchain_blocks, margin_cases
        )
        # check whether any exon intersects assembly gap in the corresponding region
        exon_gap = find_exons_gaps(
            exon_coordinates, exon_blocks, subchain_blocks, blocks_gaps, gap_coordinates
        )
        # possibly there are multiple chains
        # save data for this particular chain:
        chain_exon_gap[chain_id] = exon_gap
        chain_exon_class[chain_id] = exon_class
        chain_exon_exp_reg[chain_id] = exon_exp_region
        chain_missed[chain_id] = missing_exons
        query_sequences[chain_id] = query_seq
        query_loci[chain_id] = search_locus
        inverts[chain_id] = directed
        aa_block_sat_chain[chain_id] = aa_block_sat

        # some features that we need in case of fragmented gene
        t_start = int(chain_header[5])
        t_end = int(chain_header[6])
        q_chrom, q_start_end_str = search_locus.split(":")
        q_start_end_tup = q_start_end_str.split("-")
        q_start = int(q_start_end_tup[0])
        q_end = int(q_start_end_tup[1])
        q_strand = True if chain_header[9] == "+" else False

        # chain_qSize -> query chromosome/scaffold length
        # q_seq_len -> query sequence (that comes to CESAR) length
        fragment_data = (
            chain_id,
            q_chrom,
            q_strand,
            q_start,
            q_end,
            t_start,
            t_end,
            q_seq_len,
            chain_query_size,
        )
        fragments_data.append(fragment_data)

    if not chains:
        # it is possible in the case of "region" mode
        # possible if CESAR wrapper is used as a standalone script
        # a region given directly
        verbose("Working in the region mode")
        chains_in_input = False
        region_chrom, region_start_end = chain_file.replace(",", "").split(":")
        region_start, region_end = [int(x) for x in region_start_end.split("-")]
        if region_start < region_end:
            region_strand = True
        else:
            region_strand = False
            region_start, region_end = region_end, region_start
        search_locus = "{}:{}-{}".format(region_chrom, region_start, region_end)
        query_seq, directed = make_query_seq(
            "-1", search_locus, args["qDB"], region_strand, bed_data["strand"]
        )
        # mimic the chain parsing result:
        chain_exon_gap = None
        chain_exon_class[-1] = {}
        chain_exon_exp_reg[-1] = {}
        chain_missed[-1] = {}
        query_sequences[-1] = query_seq
        query_loci[-1] = search_locus
        inverts[-1] = directed

    if chains_in_input and args["fragments"]:
        # sort chains, get proper chain_id sorting
        if bed_data["strand"]:  # gene is + -> sort directly
            fragments_data = sorted(fragments_data, key=lambda x: x[5])
        else:  # gene is - -> reverse sort of chains
            fragments_data = sorted(fragments_data, key=lambda x: x[6], reverse=True)

        # merge query feat dictionaries
        exon_gap = merge_dicts(chain_exon_gap.values())
        exon_class = merge_dicts(chain_exon_class.values())
        exon_exp_region = merge_dicts(chain_exon_exp_reg.values())
        aa_block_sat = merge_dicts(aa_block_sat_chain.values())
        missing_exons = intersect_lists(chain_missed.values())

        query_seq_chunks = []
        for elem in fragments_data:
            # stitch query seq in a proper order; elem[0] -> chain_id
            query_seq_chunks.append(query_sequences[elem[0]])
        query_seq = "".join(query_seq_chunks)

        # remove chain_id data from dicts
        chain_ids = list(chain_exon_gap.keys())
        for chain_id in chain_ids:
            del chain_exon_gap[chain_id]
            del chain_exon_class[chain_id]
            del chain_exon_exp_reg[chain_id]
            del aa_block_sat_chain[chain_id]
            del chain_missed[chain_id]
            del query_sequences[chain_id]
        # load new values
        chain_exon_gap[Constants.FRAGMENT_CHAIN_ID] = exon_gap
        chain_exon_class[Constants.FRAGMENT_CHAIN_ID] = exon_class
        chain_exon_exp_reg[Constants.FRAGMENT_CHAIN_ID] = exon_exp_region
        aa_block_sat_chain[Constants.FRAGMENT_CHAIN_ID] = aa_block_sat
        chain_missed[Constants.FRAGMENT_CHAIN_ID] = missing_exons
        query_sequences[Constants.FRAGMENT_CHAIN_ID] = query_seq
        inverts[Constants.FRAGMENT_CHAIN_ID] = None

    # some queries might be skipped -> we can eventually skip all of them
    # which means that there is nothing to call CESAR on
    # then it's better to halt the execution
    die("No queries left") if len(query_sequences.keys()) == 0 else None
    # predict the amount of memory
    qlength_max = max([len(v) for v in query_sequences.values()])
    memory = memory_check(bed_data["block_sizes"], qlength_max, args["estimate_memory"])
    verbose(f"\nExpecting a memory consumption of: {memory} GB")
    # check whether some reference splice sites are non-canonical
    # doesn't apply to single-exon genes
    ref_ss_data = analyse_ref_ss(s_sites) if len(exon_sequences) != 1 else None
    # there are two sources of U12 introns data:
    # 1) U12 file provided at the very beginning
    # 2) If splice site in reference is non-canonical -> we also treat this as U12
    #    even if this splice site is not in the U12 data
    append_u12(args["u12"], args["gene"], ref_ss_data)
    cesar_in_data, cesar_in_filename, is_temp = make_cesar_in(prepared_exons,
                                                              query_sequences,
                                                              ref_ss_data,
                                                              args["cesar_input_save_to"],
                                                              args["temp_dir"])
    # run cesar itself
    if args.get("cesar_binary"):
        cesar_bin = args.get("cesar_binary")
    else:
        cesar_bin = DEFAULT_CESAR

    if not args["cesar_output"]:
        cesar_raw_out = run_cesar(
            cesar_in_filename,
            cesar_in_data,
            memory,
            memlim,
            cesar_bin
        )
    else:  # very specific case, load already saved CESAR output
        with open(args["cesar_output"], "r") as f:
            cesar_raw_out = f.read()
    if cesar_in_filename:
        # TODO: can be potentially simplified
        os.remove(cesar_in_filename) if cesar_in_filename else None  # wipe temp if temp
    # save raw CESAR output and close if required
    save(cesar_raw_out, args["raw_output"], t0) if args["raw_output"] else None
    # process the output, extract different features per exon
    if args["fragments"]:
        # a bit more complicated parsing of fragmented output
        proc_out = process_cesar_out__fragments(
            cesar_raw_out, fragments_data, query_loci, inverts
        )
    else:  # not fragmented: use classic procedure
        proc_out = process_cesar_out(cesar_raw_out, query_loci, inverts)
    query_exon_sequences = proc_out[0]  # sequences of predicted exons in query
    ref_exon_sequences_ali = proc_out[1]  # reference sequences -> aligned
    pIDs = proc_out[2]  # nucleotide %IDs
    pBl = proc_out[3]  # BLOSUM scores for protein sequences
    query_coords = proc_out[4]  # genomic coordinates in the query
    exon_num_corr = proc_out[5]  # in case of intron del: ref/que correspondence
    # TODO: prot_s variable is unused -> probably no need to compute it upstream?
    prot_s = proc_out[6]  # protein sequences in query
    codon_s = proc_out[7]  # dict containing sequences of ref and query codons
    aa_cesar_sat = proc_out[8]  # says whether an exon has outstanding quality
    # raw codon table \ superset of "codon_s" basically
    # after changes in TOGA1.1 is needed again
    codon_tables = proc_out[9]  
    aa_eq_len = aa_eq_len_check(exon_sequences, query_exon_sequences)

    if chains:
        # if and only if there are chains, it's possible to extract
        # another check for exceptional exons
        chain_exon_class = get_a_plus(
            chain_exon_class, aa_cesar_sat, aa_block_sat_chain, aa_eq_len
        )

    # time to arrange all these data altogether
    final_output, chain_ex_inc = arrange_output(
        args["gene"],
        ref_exon_sequences_ali,
        query_exon_sequences,
        pIDs,
        pBl,
        query_coords,
        chain_exon_gap,
        chain_exon_class,
        chain_exon_exp_reg,
        exon_num_corr,
        chain_missed,
        args["paral"],
    )
    exon_to_len = {k + 1: len(v) for k, v in exon_sequences.items()}
    verbose(f"Exon lens are: {exon_to_len}")
    # this is for inact mutations check:
    chain_to_exon_to_properties = (
        chain_exon_class,
        chain_exon_gap,
        pIDs,
        pBl,
        chain_missed,
        chain_ex_inc,
        exon_to_len,
    )
    verbose(f"Chain to exon to properties = {chain_to_exon_to_properties}")
    if args["check_loss"]:  # call inact mutations scanner,
        loss_report, del_mis_exons = inact_mut_check(
            cesar_raw_out,
            v=VERBOSE,
            gene=args["gene"],
            ex_prop=chain_to_exon_to_properties,
            ref_ss=ref_ss_data,
            sec_codons=sec_codons,
            no_fpi=args["no_fpi"],
            alt_f_del=args["alt_frame_del"],
            mask_all_first_10p=args["mask_all_first_10p"],
        )
    else:  # do not call inact mut scanner
        loss_report = None
        del_mis_exons = None
    
    # del_mis_exons contains:
    # chain_id(string) = [0-based exon nums]
    need_correct_codon_and_prot = del_mis_exons is not None and len(del_mis_exons.keys()) > 0
    if need_correct_codon_and_prot:
        # if exists, need to filter codon alignment accordingly
        # chain id is numeric in "codon_table"
        codon_s = redo_codon_sequences(codon_tables, del_mis_exons)

    prot_s = extract_prot_sequences_from_codon(args["gene"], codon_s)

    # save protein/codon ali and text output
    save_prot(prot_s, args["prot_out"])
    save_codons(args["gene"], codon_s, args["codon_out"])
    save(final_output, args["output"], t0, loss_report)
    sys.exit(0)


if __name__ == "__main__":
    t0 = dt.now()
    cmd_args = vars(parse_args())
    realign_exons(cmd_args)
    sys.exit(0)
