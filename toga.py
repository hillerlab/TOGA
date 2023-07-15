#!/usr/bin/env python3
"""Master script for the TOGA pipeline.

Perform all operations from the beginning to the end.
If you need to call TOGA: most likely this is what you need.
"""
import argparse
import sys
import os
import subprocess
import time
from datetime import datetime as dt
import json
import shutil
from math import ceil
from collections import defaultdict
import logging
from twobitreader import TwoBitFile
from modules.common import parts
from modules.filter_bed import prepare_bed_file
from modules.bed_hdf5_index import bed_hdf5_index
from modules.chain_bst_index import chain_bst_index
from modules.merge_chains_output import merge_chains_output
from modules.make_pr_pseudogenes_anno import create_ppgene_track
from modules.merge_cesar_output import merge_cesar_output
from modules.gene_losses_summary import gene_losses_summary
from modules.orthology_type_map import orthology_type_map
from modules.classify_chains import classify_chains
from modules.get_transcripts_quality import classify_transcripts
from modules.make_query_isoforms import get_query_isoforms_data
from modules.collect_prefefined_glp_classes import collect_predefined_glp_cases
from modules.collect_prefefined_glp_classes import add_transcripts_to_missing
from modules.stitch_fragments import stitch_scaffolds
from modules.common import read_isoforms_file
from modules.common import to_log
from modules.common import setup_logger
from version import __version__


__author__ = "Bogdan M. Kirilenko"
__email__ = "kirilenkobm [at] google mail"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

U12_FILE_COLS = 3
U12_AD_FIELD = {"A", "D"}
ISOFORMS_FILE_COLS = 2
NF_DIR_NAME = "nextflow_logs"
NEXTFLOW = "nextflow"
CESAR_PUSH_INTERVAL = 30  # CESAR jobs push interval
ITER_DURATION = 60  # CESAR jobs check interval
MEMLIM_ARG = "--memlim"
FRAGM_ARG = "--fragments"
LOCATION = os.path.dirname(__file__)
CESAR_RUNNER = os.path.abspath(
    os.path.join(LOCATION, "cesar_runner.py")
)  # script that will run jobs
CESAR_RUNNER_TMP = "{0} {1} {2} --check_loss {3} --rejected_log {4}"
CESAR_PRECOMPUTED_REGIONS_DIRNAME = "cesar_precomputed_regions"
CESAR_PRECOMPUTED_MEMORY_DIRNAME = "cesar_precomputed_memory"
CESAR_PRECOMPUTED_ORTHO_LOCI_DIRNAME = "cesar_precomputed_orthologous_loci"

CESAR_PRECOMPUTED_MEMORY_DATA = "cesar_precomputed_memory.tsv"
CESAR_PRECOMPUTED_REGIONS_DATA = "cesar_precomputed_regions.tsv"
CESAR_PRECOMPUTED_ORTHO_LOCI_DATA = "cesar_precomputed_orthologous_loci.tsv"

NUM_CESAR_MEM_PRECOMP_JUBS = 500


TEMP_CHAIN_CLASS = "temp_chain_trans_class"
MODULES_DIR = "modules"
RUNNING = "RUNNING"
CRASHED = "CRASHED"
TEMP = "temp"

# automatically enable flush
# print = functools.partial(print, flush=True)


class Toga:
    """TOGA manager class."""

    def __init__(self, args):
        """Initiate toga class."""
        self.t0 = dt.now()
        # define project name
        if args.project_name:
            self.project_name = args.project_name
        elif args.project_dir:
            _dirname = os.path.dirname(args.project_dir)
            self.project_name = os.path.basename(_dirname)
        else:
            self.project_name = self.__gen_project_name()
        # create project dir
        self.wd = (
            os.path.abspath(args.project_dir)
            if args.project_dir
            else os.path.join(os.getcwd(), self.project_name)
        )
        os.mkdir(self.wd) if not os.path.isdir(self.wd) else None

        # manage logfiles
        _log_filename = self.t0.strftime("%Y.%m.%d_at_%H:%M")
        self.log_file = os.path.join(self.wd, f"toga_{_log_filename}.log")
        self.log_dir = os.path.join(self.wd, "temp_logs")  # temp file to collect logs from processes
        os.mkdir(self.log_dir) if not os.path.isdir(self.log_dir) else None
        setup_logger(self.log_file)

        # check if all files TOGA needs are here
        self.temp_files = []  # remove at the end, list of temp files
        to_log("#### Initiating TOGA class ####")
        self.para = args.para
        self.__check_args_correctness(args)
        self.__modules_addr()
        self.__check_dependencies()
        self.__check_completeness()
        self.toga_exe_path = os.path.dirname(__file__)
        self.version = self.__get_version()
        self.para_bigmem = args.para_bigmem
        self.nextflow_dir = self.__get_nf_dir(args.nextflow_dir)
        self.nextflow_config_dir = args.nextflow_config_dir

        self.temp_wd = os.path.join(self.wd, TEMP)
        self.project_name = self.project_name.replace("/", "")
        os.mkdir(self.temp_wd) if not os.path.isdir(self.temp_wd) else None

        # TODO: maybe simplify the pipeline and exclude bigmem config?
        if args.cesar_bigmem_config:
            self.nextflow_bigmem_config = os.path.abspath(args.cesar_bigmem_config)
        else:  # if none: we cannot call os.path.abspath method
            self.nextflow_bigmem_config = None
        self.__check_nf_config()

        # to avoid crash on filesystem without locks:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        chain_basename = os.path.basename(args.chain_input)
        # dir to collect log files with rejected reference genes:
        self.rejected_dir = os.path.join(self.temp_wd, "rejected")
        os.mkdir(self.rejected_dir) if not os.path.isdir(self.rejected_dir) else None

        # filter chain in this folder
        g_ali_basename = "genome_alignment"
        self.chain_file = os.path.join(self.temp_wd, f"{g_ali_basename}.chain")
        # there is an assumption that chain file has .chain extension
        # chain indexing was a bit problematic: (i) bsddb3 fits perfectly but is very
        # painful to install, (ii) sqlite is also fine but might be dysfunctional on some
        # cluster file systems, so we create chain_ID: (start_byte, offset) dictionary for
        # instant extraction of a particular chain from the chain file
        # we save these dictionaries into two files: a text file (tsv) and binary file with BST
        # depending on the case we will use both (for maximal performance)
        self.chain_index_file = os.path.join(self.temp_wd, f"{g_ali_basename}.bst")
        self.chain_index_txt_file = os.path.join(
            self.temp_wd, f"{g_ali_basename}.chain_ID_position"
        )

        # make the command, prepare the chain file
        if not os.path.isfile(args.chain_input):
            chain_filter_cmd = None
            self.die(f"Error! File {args.chain_input} doesn't exist!")
        elif chain_basename.endswith(".gz"):  # version for gz
            chain_filter_cmd = (
                f"gzip -dc {args.chain_input} | "
                f"{self.CHAIN_SCORE_FILTER} stdin "
                f"{args.min_score} > {self.chain_file}"
                # Tried to replace C binary with AWK, something to think about
                # f"awk -f {self.CHAIN_SCORE_FILTER_AWK} {args.min_score} "
                # f"> {self.chain_file}"
            )
        elif args.no_chain_filter:  # it is .chain and score filter is not required
            chain_filter_cmd = f"rsync -a {args.chain_input} {self.chain_file}"
        else:  # it is .chain | score filter required
            chain_filter_cmd = (
                f"{self.CHAIN_SCORE_FILTER} {args.chain_input} "
                f"{args.min_score} > {self.chain_file}"
            )

        # filter chains with score < threshold
        self.__call_proc(
            chain_filter_cmd, "Please check if you use a proper chain file."
        )

        # bed define bed files addresses
        self.ref_bed = os.path.join(self.temp_wd, "toga_filt_ref_annot.bed")
        self.index_bed_file = os.path.join(self.temp_wd, "toga_filt_ref_annot.hdf5")

        # filter bed file
        bed_filt_rejected_file = "BED_FILTER_REJECTED.txt"
        bed_filt_rejected = os.path.join(self.rejected_dir, bed_filt_rejected_file)
        # keeping UTRs!
        prepare_bed_file(
            args.bed_input,
            self.ref_bed,
            save_rejected=bed_filt_rejected,
            only_chrom=args.limit_to_ref_chrom,
        )

        # mics things
        self.isoforms_arg = args.isoforms if args.isoforms else None
        self.isoforms = None  # will be assigned after completeness check
        self.chain_jobs = args.chain_jobs_num
        self.cesar_binary = (
            self.DEFAULT_CESAR if not args.cesar_binary else args.cesar_binary
        )
        self.opt_cesar_binary = os.path.abspath(
            os.path.join(LOCATION, "cesar_input_optimiser.py")
        )
        self.cesar_is_opt = args.using_optimized_cesar
        self.output_opt_cesar_regions = args.output_opt_cesar_regions
        self.time_log = args.time_marks
        self.stop_at_chain_class = args.stop_at_chain_class
        self.rejected_log = os.path.join(self.wd, "genes_rejection_reason.tsv")

        self.keep_temp = True if args.keep_temp else False
        # define to call CESAR or not to call
        self.t_2bit = self.__find_two_bit(args.tDB)
        self.q_2bit = self.__find_two_bit(args.qDB)

        self.hq_orth_threshold = 0.95
        self.cesar_jobs_num = args.cesar_jobs_num
        self.cesar_buckets = args.cesar_buckets
        self.cesar_mem_limit = args.cesar_mem_limit
        self.cesar_chain_limit = args.cesar_chain_limit
        self.uhq_flank = args.uhq_flank
        self.mask_stops = args.mask_stops
        self.no_fpi = args.no_fpi
        self.o2o_only = args.o2o_only
        self.annotate_paralogs = args.annotate_paralogs
        self.keep_nf_logs = args.do_not_del_nf_logs
        self.exec_cesar_parts_sequentially = args.cesar_exec_seq
        self.ld_model_arg = args.ld_model
        self.mask_all_first_10p = args.mask_all_first_10p

        self.cesar_ok_merged = (
            None  # Flag: indicates whether any cesar job BATCHES crashed
        )
        self.crashed_cesar_jobs = []  # List of individual cesar JOBS that crashed
        self.cesar_crashed_batches_log = os.path.join(
            self.temp_wd, "_cesar_crashed_job_batches.txt"
        )
        self.cesar_crashed_jobs_log = os.path.join(
            self.temp_wd, "_cesar_crashed_jobs.txt"
        )
        self.fragmented_genome = False if args.disable_fragments_joining else True
        self.orth_score_threshold = args.orth_score_threshold
        if self.orth_score_threshold < 0.0 or args.orth_score_threshold > 1.0:
            self.die(
                "orth_score_threshold parameter must be in range [0..1], got "
                f"{self.orth_score_threshold}; Abort"
            )

        self.chain_results_df = os.path.join(self.temp_wd, "chain_results_df.tsv")
        self.nucl_fasta = os.path.join(self.wd, "nucleotide.fasta")
        self.prot_fasta = os.path.join(self.wd, "prot.fasta")
        self.codon_fasta = os.path.join(self.wd, "codon.fasta")
        self.meta_data = os.path.join(self.temp_wd, "exons_meta_data.tsv")
        self.intermediate_bed = os.path.join(self.temp_wd, "intermediate.bed")
        self.orthology_type = os.path.join(self.wd, "orthology_classification.tsv")
        self.trash_exons = os.path.join(self.temp_wd, "trash_exons.bed")
        self.gene_loss_data = os.path.join(self.temp_wd, "inact_mut_data")
        self.query_annotation = os.path.join(self.wd, "query_annotation.bed")
        self.loss_summ = os.path.join(self.wd, "loss_summ_data.tsv")
        # directory to store intermediate files with technically unprocessable transcripts:
        self.technical_cesar_err = os.path.join(self.temp_wd, "technical_cesar_err")
        # unprocessed transcripts to be considered Missing:
        self.technical_cesar_err_merged = os.path.join(
            self.temp_wd, "technical_cesar_err.txt"
        )

        self.bed_fragm_exons_data = os.path.join(
            self.temp_wd, "bed_fragments_to_exons.tsv"
        )
        self.precomp_mem_cesar = os.path.join(
            self.temp_wd, CESAR_PRECOMPUTED_MEMORY_DATA
        )
        self.precomp_reg_dir = None
        self.cesar_mem_was_precomputed = False
        self.u12_arg = args.u12
        self.u12 = None  # assign after U12 file check

        # genes to be classified as missing
        self._transcripts_not_intersected = []
        self._transcripts_not_classified = []
        self.predefined_glp_cesar_split = os.path.join(
            self.temp_wd, "predefined_glp_cesar_split.tsv"
        )

        self.__check_param_files()

        # create symlinks to 2bits: let user know what 2bits were used
        self.t_2bit_link = os.path.join(self.wd, "t2bit.link")
        self.q_2bit_link = os.path.join(self.wd, "q2bit.link")
        self.__make_symlink(self.t_2bit, self.t_2bit_link)
        self.__make_symlink(self.q_2bit, self.q_2bit_link)

        # dump input parameters, object state
        self.toga_params_file = os.path.join(self.temp_wd, "toga_init_state.json")
        self.toga_args_file = os.path.join(self.wd, "project_args.json")
        self.version_file = os.path.join(self.wd, "version.txt")

        with open(self.toga_params_file, "w") as f:
            # default=string is a workaround to serialize datetime object
            json.dump(self.__dict__, f, default=str)
        with open(self.toga_args_file, "w") as f:
            json.dump(vars(args), f, default=str)
        with open(self.version_file, "w") as f:
            f.write(self.version)

        to_log(f"Saving output to {self.wd}")
        to_log(f"Arguments stored in {self.toga_args_file}")

    @staticmethod
    def __make_symlink(src, dest):
        """Create a symlink.

        os.symlink expects that dest doesn't exist.
        Need to make a couple of checks before calling it.
        """
        if os.path.islink(dest):
            return
        elif os.path.isfile(dest):
            return
        os.symlink(src, dest)

    @staticmethod
    def __gen_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name

    def __check_args_correctness(self, args):
        """Check that arguments are correct.

        Error exit if any argument is wrong.
        """
        if args.cesar_buckets:
            # if set, need to check that it could be split into numbers
            comma_sep = args.cesar_buckets.split(",")
            all_numeric = [x.isnumeric() for x in comma_sep]
            if any(x is False for x in all_numeric):
                # there is some non-numeric value
                err_msg = (
                    f"Error! --cesar_buckets value {args.cesar_buckets} is incorrect\n"
                    f"Expected comma-separated list of integers"
                )
                self.die(err_msg)
        if not os.path.isfile(args.chain_input):
            self.die(f"Error! Chain file {args.chain_input} does not exist!")
        if not os.path.isfile(args.bed_input):
            self.die(f"Error! Bed file {args.bed_input} does not exist!")
        return

    def __check_param_files(self):
        """Check that all parameter files exist."""
        files_to_check = [
            self.u12,
            self.t_2bit,
            self.q_2bit,
            self.cesar_binary,
            self.ref_bed,
            self.chain_file,
            self.isoforms_arg,
        ]
        for item in files_to_check:
            if not item:
                # this file just not given
                continue
            elif not os.path.isfile(item):
                self.die(f"Error! File {item} not found!")

        # sanity checks: check that bed file chroms match reference 2bit
        with open(self.ref_bed, "r") as f:
            lines = [line.rstrip().split("\t") for line in f]
            t_in_bed = set(x[3] for x in lines)
            chroms_in_bed = set(x[0] for x in lines)
            # 2bit check function accepts a dict chrom: size
            # from bed12 file we cannot infer sequence length
            # None is just a placeholder that indicated that we don't need
            # to compare chrom lengths with 2bit
            chrom_sizes_in_bed = {x: None for x in chroms_in_bed}
        self.__check_isoforms_file(t_in_bed)
        self.__check_u12_file(t_in_bed)
        self.__check_2bit_file(self.t_2bit, chrom_sizes_in_bed, self.ref_bed)
        # need to check that chain chroms and their sizes match 2bit file data
        with open(self.chain_file, "r") as f:
            header_lines = [x.rstrip().split() for x in f if x.startswith("chain")]
            t_chrom_to_size = {x[2]: int(x[3]) for x in header_lines}
            q_chrom_to_size = {x[7]: int(x[8]) for x in header_lines}
        # q_chrom_in_chain = set(x.split()[7] for x in f if x.startswith("chain"))
        f.close()
        self.__check_2bit_file(self.t_2bit, t_chrom_to_size, self.chain_file)
        self.__check_2bit_file(self.q_2bit, q_chrom_to_size, self.chain_file)

    def __get_nf_dir(self, nf_dir_arg):
        """Define nextflow directory."""
        if nf_dir_arg is None:
            default_dir = os.path.join(self.LOCATION, NF_DIR_NAME)
            os.mkdir(default_dir) if not os.path.isdir(default_dir) else None
            return default_dir
        else:
            os.mkdir(nf_dir_arg) if not os.path.isdir(nf_dir_arg) else None
            return nf_dir_arg

    def __check_2bit_file(self, two_bit_file, chroms_sizes, chrom_file):
        """Check that 2bit file is readable."""
        try:  # try to catch EOFError: if 2bitreader cannot read file
            two_bit_reader = TwoBitFile(two_bit_file)
            # check what sequences are in the file:
            twobit_seq_to_size = two_bit_reader.sequence_sizes()
            twobit_sequences = set(twobit_seq_to_size.keys())
            to_log(f"Found {len(twobit_sequences)} sequences in {two_bit_file}")
        except EOFError as err:  # this is a file but twobit reader couldn't read it
            twobit_seq_to_size = None  # to suppress linter
            twobit_sequences = None
            to_log(str(err))
            to_log(f"Error! Twobitreader cannot open {two_bit_file}")
            self.die("Abort")
        # another check: that bed or chain chromosomes intersect 2bit file sequences
        check_chroms = set(chroms_sizes.keys())  # chroms in the input file
        intersection = twobit_sequences.intersection(check_chroms)
        chroms_not_in_2bit = check_chroms.difference(twobit_sequences)

        if len(chroms_not_in_2bit) > 0:
            # err = (
            #     f"Error! 2bit file: {two_bit_file}; chain/bed file: {chrom_file}; "
            #     f"Different sets of chromosomes!"
            # )
            missing_top_100 = list(chroms_not_in_2bit)[:100]
            missing_str = "\n".join(missing_top_100)
            err = (
                f"Error! 2bit file: {two_bit_file}; chain/bed file: {chrom_file}; "
                f"Some chromosomes present in the chain/bed file are not found in the "
                f"Two bit file. First <=100: {missing_str}"
            )
            self.die(err)
        # check that sizes also match
        for chrom in intersection:
            twobit_seq_len = twobit_seq_to_size[chrom]
            comp_file_seq_len = chroms_sizes[chrom]
            # if None: this is from bed file: cannot compare
            if comp_file_seq_len is None:
                continue
            if twobit_seq_len == comp_file_seq_len:
                continue
            # got different sequence length in chain and 2bit files
            # which means these chains come from something different
            err = (
                f"Error! 2bit file: {two_bit_file}; chain_file: {chrom_file} "
                f"Chromosome: {chrom}; Sizes don't match! "
                f"Size in twobit: {twobit_seq_len}; size in chain: {comp_file_seq_len}"
            )
            self.die(err)
        return

    def __check_u12_file(self, t_in_bed):
        """Sanity check for U12 file."""
        if not self.u12_arg:
            # just not provided: nothing to check
            return
        # U12 file provided
        self.u12 = os.path.join(self.temp_wd, "u12_data.txt")
        filt_lines = []
        f = open(self.u12_arg, "r")
        for num, line in enumerate(f, 1):
            line_data = line.rstrip().split("\t")
            if len(line_data) != U12_FILE_COLS:
                err_msg = (
                    f"Error! U12 file {self.u12} line {num} is corrupted, 3 fields expected; "
                    f"Got {len(line_data)}; please note that a tab-separated file expected"
                )
                self.die(err_msg)
            trans_id = line_data[0]
            if trans_id not in t_in_bed:
                # transcript doesn't appear in the bed file: skip it
                continue
            exon_num = line_data[1]
            if not exon_num.isnumeric():
                err_msg = (
                    f"Error! U12 file {self.u12} line {num} is corrupted, "
                    f"field 2 value is {exon_num}; This field must "
                    f"contain a numeric value (exon number)."
                )
                self.die(err_msg)
            acc_don = line_data[2]
            if acc_don not in U12_AD_FIELD:
                err_msg = (
                    f"Error! U12 file {self.u12} line {num} is corrupted, field 3 value is {acc_don}"
                    f"; This field could have either A or D value."
                )
                self.die(err_msg)
            filt_lines.append(line)  # save this line
        f.close()
        # another check: what if there are no lines after filter?
        if len(filt_lines) == 0:
            err_msg = (
                f"Error! No lines left in the {self.u12_arg} file after filter."
                f"Please check that transcript IDs in this file and bed {self.ref_bed} are consistent"
            )
            self.die(err_msg)
        with open(self.u12, "w") as f:
            f.write("".join(filt_lines))

    def __check_isoforms_file(self, t_in_bed):
        """Sanity checks for isoforms file."""
        if not self.isoforms_arg:
            to_log("Continue without isoforms file: not provided")
            return  # not provided: nothing to check
        # isoforms file provided: need to check correctness and completeness
        # then check isoforms file itself
        _, isoform_to_gene, header = read_isoforms_file(self.isoforms_arg)
        header_maybe_gene = header[
            0
        ]  # header is optional, if not the case: first field is a gene
        header_maybe_trans = header[1]  # and the second is the isoform
        # save filtered isoforms file here:  (without unused transcripts)
        self.isoforms = os.path.join(self.temp_wd, "isoforms.tsv")
        # this set contains isoforms found in the isoforms file
        t_in_i = set(isoform_to_gene.keys())
        # there are transcripts that appear in bed but not in the isoforms file
        # if this set is non-empty: raise an error
        u_in_b = t_in_bed.difference(t_in_i)

        if len(u_in_b) != 0:  # isoforms file is incomplete
            extra_t_list = "\n".join(
                list(u_in_b)[:100]
            )  # show first 100 (or maybe show all?)
            err_msg = (
                f"Error! There are {len(u_in_b)} transcripts in the bed "
                f"file absent in the isoforms file! "
                f"There are the transcripts (first 100):\n{extra_t_list}"
            )
            self.die(err_msg)

        t_in_both = t_in_bed.intersection(t_in_i)  # isoforms data that we save
        # if header absent: second field found in the bed file
        # then we don't need to write the original header
        # if present -> let keep it
        # there is not absolutely correct: MAYBE there is no header at all, but
        # the first line of the isoforms file is not in the bed file so we still will write it
        skip_header = header_maybe_trans in t_in_bed

        # write isoforms file
        f = open(self.isoforms, "w")
        if not skip_header:
            f.write(f"{header_maybe_gene}\t{header_maybe_trans}\n")
        to_log(f"Writing isoforms data for {len(t_in_both)} transcripts.")
        for trans in t_in_both:
            gene = isoform_to_gene[trans]
            f.write(f"{gene}\t{trans}\n")

    def __log_blank(self):
        """Add blank line to log."""
        with open(self.log_file, "a") as f:
            f.write("\n\n")

    def die(self, msg, rc=1):
        """Show msg in stderr, exit with the rc given."""
        self.__log_blank()
        to_log(msg)
        to_log(f"Program finished with exit code {rc}\n")
        # for t_file in self.temp_files:  # remove temp files if required
        #     os.remove(t_file) if os.path.isfile(t_file) and not self.keep_temp else None
        #     shutil.rmtree(t_file) if os.path.isdir(t_file) and not self.keep_temp else None
        self.__mark_crashed()
        sys.exit(rc)

    def __modules_addr(self):
        """Define addresses of modules."""
        self.LOCATION = os.path.dirname(__file__)  # folder containing pipeline scripts
        self.CONFIGURE = os.path.join(self.LOCATION, "configure.sh")
        self.CHAIN_SCORE_FILTER = os.path.join(
            self.LOCATION, "modules", "chain_score_filter"
        )
        self.CHAIN_SCORE_FILTER_AWK = os.path.join(
            self.LOCATION, "modules", "chain_score_filter.awk"
        )
        self.CHAIN_COORDS_CONVERT_LIB = os.path.join(
            self.LOCATION, "modules", "chain_coords_converter_slib.so"
        )
        self.EXTRACT_SUBCHAIN_LIB = os.path.join(
            self.LOCATION, "modules", "extract_subchain_slib.so"
        )
        self.CHAIN_FILTER_BY_ID = os.path.join(
            self.LOCATION, "modules", "chain_filter_by_id"
        )
        self.CHAIN_BDB_INDEX = os.path.join(
            self.LOCATION, MODULES_DIR, "chain_bst_index.py"
        )
        self.CHAIN_INDEX_SLIB = os.path.join(
            self.LOCATION, MODULES_DIR, "chain_bst_lib.so"
        )
        self.BED_BDB_INDEX = os.path.join(
            self.LOCATION, MODULES_DIR, "bed_hdf5_index.py"
        )
        self.SPLIT_CHAIN_JOBS = os.path.join(self.LOCATION, "split_chain_jobs.py")
        self.MERGE_CHAINS_OUTPUT = os.path.join(
            self.LOCATION, MODULES_DIR, "merge_chains_output.py"
        )
        self.CLASSIFY_CHAINS = os.path.join(
            self.LOCATION, MODULES_DIR, "classify_chains.py"
        )
        self.SPLIT_EXON_REALIGN_JOBS = os.path.join(
            self.LOCATION, "split_exon_realign_jobs.py"
        )
        self.MERGE_CESAR_OUTPUT = os.path.join(
            self.LOCATION, MODULES_DIR, "merge_cesar_output.py"
        )
        self.TRANSCRIPT_QUALITY = os.path.join(
            self.LOCATION, MODULES_DIR, "get_transcripts_quality.py"
        )
        self.GENE_LOSS_SUMMARY = os.path.join(
            self.LOCATION, MODULES_DIR, "gene_losses_summary.py"
        )
        self.ORTHOLOGY_TYPE_MAP = os.path.join(
            self.LOCATION, MODULES_DIR, "orthology_type_map.py"
        )
        self.PRECOMPUTE_OPT_CESAR_DATA = os.path.join(
            self.LOCATION, MODULES_DIR, "precompute_regions_for_opt_cesar.py"
        )
        self.MODEL_TRAINER = os.path.join(self.LOCATION, "train_model.py")
        self.DEFAULT_CESAR = os.path.join(self.LOCATION, "CESAR2.0", "cesar")
        self.nextflow_rel_ = os.path.join(self.LOCATION, "execute_joblist.nf")
        self.NF_EXECUTE = os.path.abspath(self.nextflow_rel_)

    def __check_dependencies(self):
        """Check all dependencies."""
        c_not_compiled = any(
            os.path.isfile(f) is False
            for f in [
                self.CHAIN_SCORE_FILTER,
                self.CHAIN_COORDS_CONVERT_LIB,
                self.CHAIN_FILTER_BY_ID,
                self.EXTRACT_SUBCHAIN_LIB,
                self.CHAIN_INDEX_SLIB,
            ]
        )
        if c_not_compiled:
            to_log("Warning! C code is not compiled, trying to compile...")
        imports_not_found = False
        try:
            import twobitreader
            import networkx
            import pandas
            import xgboost
            import joblib
            import h5py
        except ImportError:
            to_log("Warning! Some of the required packages are not installed.")
            imports_not_found = True

        not_nf = shutil.which(NEXTFLOW) is None
        if not self.para and not_nf:
            msg = (
                "Error! Cannot fild nextflow executable. Please make sure you "
                "have a nextflow binary in a directory listed in your $PATH"
            )
            self.die(msg)

        not_all_found = any([c_not_compiled, imports_not_found])
        self.__call_proc(
            self.CONFIGURE, "Could not call configure.sh!"
        ) if not_all_found else None

    def __check_completeness(self):
        """Check if all modules are presented."""
        files_must_be = [
            self.CONFIGURE,
            self.CHAIN_BDB_INDEX,
            self.BED_BDB_INDEX,
            self.SPLIT_CHAIN_JOBS,
            self.MERGE_CHAINS_OUTPUT,
            self.CLASSIFY_CHAINS,
            self.SPLIT_EXON_REALIGN_JOBS,
            self.MERGE_CESAR_OUTPUT,
            self.GENE_LOSS_SUMMARY,
            self.ORTHOLOGY_TYPE_MAP,
        ]
        for _file in files_must_be:
            if os.path.isfile(_file):
                continue
            self.die(f"Error! File {_file} not found!")

    def __check_nf_config(self):
        """Check that nextflow configure files are here."""
        if self.nextflow_bigmem_config and not os.path.isfile(
            self.nextflow_bigmem_config
        ):
            # bigmem config is special for now | sanity check -> defined and no file -> crash
            self.die(f"Error! File {self.nextflow_bigmem_config} not found!")

        if self.nextflow_config_dir is None:
            # no nextflow config provided -> using local executor
            self.cesar_config_template = None
            self.nf_chain_extr_config_file = None
            self.local_executor = True
            return
        # check conflict: if we set para; nf args should not be set
        if self.para and self.nextflow_config_dir:
            self.die(
                "Conflict: --para and --nf_dir should not be used at the same time"
            )
        # check that required config files are here
        if not os.path.isdir(self.nextflow_config_dir):
            self.die(
                f"Error! Nextflow config dir {self.nextflow_config_dir} does not exist!"
            )
        err_msg = (
            "Please note these two files are expected in the nextflow config directory:\n"
            "1) call_cesar_config_template.nf"
            "2) extract_chain_features_config.nf"
        )
        # check CESAR config template first
        nf_cesar_config_temp = os.path.join(
            self.nextflow_config_dir, "call_cesar_config_template.nf"
        )
        if not os.path.isfile(nf_cesar_config_temp):
            self.die(f"Error! File {nf_cesar_config_temp} not found!\n{err_msg}")
        # we need the content of this file
        with open(nf_cesar_config_temp, "r") as f:
            self.cesar_config_template = f.read()
        # check chain extract features config; we need abspath to this file
        self.nf_chain_extr_config_file = os.path.abspath(
            os.path.join(self.nextflow_config_dir, "extract_chain_features_config.nf")
        )
        if not os.path.isfile(self.nf_chain_extr_config_file):
            self.die(
                f"Error! File {self.nf_chain_extr_config_file} not found!\n{err_msg}"
            )
        self.local_executor = False

    def __call_proc(self, cmd, extra_msg=None):
        """Call a subprocess and catch errors."""
        to_log(f"Calling cmd:\n{cmd}\n")
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            to_log(extra_msg) if extra_msg else None
            self.die(f"Error! Process:\n{cmd}\ndied! Abort.")
        to_log("Command finished with exit code 0.")

    def __find_two_bit(self, db):
        """Find a 2bit file."""
        if os.path.isfile(db):
            return os.path.abspath(db)
        # For now here is a hillerlab-oriented solution
        # you can write your own template for 2bit files location
        with_alias = f"/projects/hillerlab/genome/gbdb-HL/{db}/{db}.2bit"
        if os.path.isfile(with_alias):
            return with_alias
        elif os.path.islink(with_alias):
            return os.path.abspath(os.readlink(with_alias))
        self.die(f"Two bit file {db} not found! Abort")

    def run(self):
        """Run toga. Method to be called."""
        # 0) preparation:
        # define the project name and mkdir for it
        # move chain file filtered
        # define initial values
        # make indexed files for the chain
        self.__mark_start()
        to_log("\n\n#### STEP 0: making chain and bed file indexes\n")
        self.__make_indexed_chain()
        self.__make_indexed_bed()
        self.__time_mark("Made indexes")

        # 1) make joblist for chain features extraction
        to_log("\n\n#### STEP 1: Generate extract chain features jobs\n")
        self.__split_chain_jobs()
        self.__time_mark("Split chain jobs")

        # 2) extract chain features: parallel process
        to_log("\n\n#### STEP 2: Extract chain features: parallel step\n")
        self.__extract_chain_features()
        self.__time_mark("Chain jobs done")
        self.__collapse_logs("chain_runner_")

        # 3) create chain features dataset
        to_log("\n\n#### STEP 3: Merge step 2 output\n")
        self.__merge_chains_output()
        self.__time_mark("Chains output merged")

        # 4) classify chains as orthologous, paralogous, etc using xgboost
        to_log("\n\n#### STEP 4: Classify chains using gradient boosting model\n")
        self.__classify_chains()
        self.__time_mark("Chains classified")

        # 5) create cluster jobs for CESAR2.0
        to_log("\n\n#### STEP 5: Generate CESAR jobs")
        # experimental feature, not publically available:
        # self.__precompute_data_for_opt_cesar()
        self.__split_cesar_jobs()
        self.__time_mark("Split cesar jobs done")

        # 6) Create bed track for processed pseudogenes
        to_log("\n\n#### STEP 6: Create processed pseudogenes track\n")
        self.__get_proc_pseudogenes_track()

        # 7) call CESAR jobs: parallel step
        to_log("\n\n### STEP 7: Execute CESAR jobs: parallel step\n")
        self.__run_cesar_jobs()
        self.__time_mark("Cesar jobs done")
        self.__check_cesar_completeness()

        # 8) parse CESAR output, create bed / fasta files
        to_log("\n\n#### STEP 8: Merge STEP 7 output\n")
        self.__merge_cesar_output()
        self.__time_mark("Merged cesar output")

        # 9) classify projections/genes as lost/intact
        # also measure projections confidence levels
        to_log("\n\n#### STEP 9: Gene loss pipeline classification\n")
        self.__transcript_quality()  # maybe remove -> not used anywhere
        self.__gene_loss_summary()
        self.__time_mark("Got gene loss summary")

        # TODO: missing genes! no chain chrom???
        # 10) classify genes as one2one, one2many, etc orthologs
        to_log("\n\n#### STEP 10: Create orthology relationships table\n")
        self.__orthology_type_map()

        # 11) merge logs containing information about skipped genes,transcripts, etc.
        to_log("\n\n#### STEP 11: Cleanup: merge parallel steps output files")
        self.__merge_split_files()
        shutil.rmtree(self.log_dir)
        self.__check_crashed_cesar_jobs()
        # Everything is done

        self.__time_mark("Everything is done")
        if not self.cesar_ok_merged:
            cesar_not_ok_message = (
                f"PLEASE NOTE:\nCESAR RESULTS ARE LIKELY INCOMPLETE"
                f"Please look at: {self.cesar_crashed_batches_log}"
            )
            to_log(cesar_not_ok_message)
        self.__left_done_mark()

    def __collapse_logs(self, prefix):
        """Merge logfiles starting with prefix into a single log."""
        log_filenames_with_prefix = [x for x in os.listdir(self.log_dir) if x.startswith(prefix)]
        log_f = open(self.log_file, "a")
        for log_filename in log_filenames_with_prefix:
            full_path = os.path.join(self.log_dir, log_filename)
            clipped_filename = log_filename.split(".")[0]  # remove .log
            in_f = open(full_path, "r")
            for line in in_f:
                log_f.write(f"{clipped_filename}: {line}")
            in_f.close()
        log_f.close()

    def __mark_start(self):
        """Indicate that TOGA process have started."""
        p_ = os.path.join(self.wd, RUNNING)
        f = open(p_, "w")
        now_ = str(dt.now())
        f.write(f"TOGA process started at {now_}\n")
        f.close()

    def __mark_crashed(self):
        """Indicate that TOGA process died."""
        running_f = os.path.join(self.wd, RUNNING)
        crashed_f = os.path.join(self.wd, CRASHED)
        os.remove(running_f) if os.path.isfile(running_f) else None
        f = open(crashed_f, "w")
        now_ = str(dt.now())
        f.write(f"TOGA CRASHED AT {now_}\n")
        f.close()

    def __make_indexed_chain(self):
        """Make chain index file."""
        # make *.bb file
        to_log("Started chain indexing...")
        chain_bst_index(
            self.chain_file, self.chain_index_file, txt_index=self.chain_index_txt_file
        )
        self.temp_files.append(self.chain_index_file)
        self.temp_files.append(self.chain_file)
        self.temp_files.append(self.chain_index_txt_file)

    def __time_mark(self, msg):
        """Left time mark."""
        if self.time_log is None:
            return
        t = dt.now() - self.t0
        with open(self.time_log, "a") as f:
            f.write(f"{msg} at {t}\n")

    def __make_indexed_bed(self):
        """Create gene_ID: bed line bdb indexed file."""
        to_log("Started bed file indexing...")
        bed_hdf5_index(self.ref_bed, self.index_bed_file)
        self.temp_files.append(self.index_bed_file)

    @staticmethod
    def __get_fst_col(path):
        """Just extract first file column."""
        with open(path, "r") as f:
            ret = [line.rstrip().split("\t")[0] for line in f]
        return ret

    def __split_chain_jobs(self):
        """Wrap split_jobs.py script."""
        # define arguments
        # save split jobs
        self.ch_cl_jobs = os.path.join(self.temp_wd, "chain_classification_jobs")
        # for raw results of this stage
        self.chain_class_results = os.path.join(
            self.temp_wd, "chain_classification_results"
        )
        self.chain_cl_jobs_combined = os.path.join(
            self.temp_wd, "chain_class_jobs_combined"
        )
        rejected_filename = "SPLIT_CHAIN_REJ.txt"
        rejected_path = os.path.join(self.rejected_dir, rejected_filename)
        self.temp_files.append(self.ch_cl_jobs)
        self.temp_files.append(self.chain_class_results)
        self.temp_files.append(self.chain_cl_jobs_combined)

        split_jobs_cmd = (
            f"{self.SPLIT_CHAIN_JOBS} "
            f"{self.chain_file} "
            f"{self.ref_bed} "
            f"{self.index_bed_file} "
            f"--log_file {self.log_file} "
            f"--parallel_logs_dir {self.log_dir} "
            f"--jobs_num {self.chain_jobs} "
            f"--jobs {self.ch_cl_jobs} "
            f"--jobs_file {self.chain_cl_jobs_combined} "
            f"--results_dir {self.chain_class_results} "
            f"--rejected {rejected_path}"
        )

        self.__call_proc(split_jobs_cmd, "Could not split chain jobs!")
        # collect transcripts not intersected at all here
        self._transcripts_not_intersected = self.__get_fst_col(rejected_path)

    def __extract_chain_features(self):
        """Execute extract chain features jobs."""
        # get timestamp to name the project and create a dir for that
        #  time() returns something like: 1595861493.8344169
        timestamp = str(time.time()).split(".")[0]
        project_name = f"{self.project_name}_chain_feats_at_{timestamp}"
        project_path = os.path.join(self.nextflow_dir, project_name)
        to_log(f"Extracting chain features, NF project name: {project_name}")
        to_log(f"Project path: {project_path}")

        if self.para:  # run jobs with para, skip nextflow
            # Hillerlab-specific script to schedule slurm jobs
            cmd = f'para make {project_name} {self.chain_cl_jobs_combined} -q="short"'
            to_log(f"Calling parallel process with:\n{cmd}")
            rc = subprocess.call(cmd, shell=True)
        else:  # calling jobs with nextflow
            cmd = (
                f"nextflow {self.NF_EXECUTE} "
                f"--joblist {self.chain_cl_jobs_combined}"
            )
            if not self.local_executor:
                # not local executor -> provided config files
                # need abspath for nextflow execution
                cmd += f" -c {self.nf_chain_extr_config_file}"
            to_log(f"Calling parallel process with:\n{cmd}")
            os.mkdir(project_path) if not os.path.isdir(project_path) else None
            rc = subprocess.call(cmd, shell=True, cwd=project_path)

        if rc != 0:  # if process (para or nf) died: terminate execution
            self.die(f"Error! Process {cmd} died")
        if not self.keep_nf_logs and not self.para:
            # remove nextflow intermediate files
            # if para: this dir doesn't exist
            shutil.rmtree(project_path) if os.path.isdir(project_path) else None

    def __merge_chains_output(self):
        """Call parse results."""
        # define where to save intermediate table
        merge_chains_output(
            self.ref_bed, self.isoforms, self.chain_class_results, self.chain_results_df
        )
        # .append(self.chain_results_df)  -> UCSC plugin needs that

    def __classify_chains(self):
        """Run decision tree."""
        # define input and output."""
        to_log("Classifying chains")
        self.transcript_to_chain_classes = os.path.join(self.temp_wd, "trans_to_chain_classes.tsv")
        self.pred_scores = os.path.join(self.wd, "orthology_scores.tsv")
        self.se_model = os.path.join(self.LOCATION, "models", "se_model.dat")
        self.me_model = os.path.join(self.LOCATION, "models", "me_model.dat")
        self.ld_model = os.path.join(
            self.LOCATION, "long_distance_model", "long_dist_model.dat"
        )
        cl_rej_log = os.path.join(self.rejected_dir, "classify_chains_rejected.txt")
        ld_arg_ = self.ld_model if self.ld_model_arg else None

        if not os.path.isfile(self.se_model) or not os.path.isfile(self.me_model):
            self.__call_proc(self.MODEL_TRAINER, "Models not found, training...")

        classify_chains(
            self.chain_results_df,
            self.transcript_to_chain_classes,
            self.se_model,
            self.me_model,
            rejected=cl_rej_log,
            raw_out=self.pred_scores,
            annot_threshold=self.orth_score_threshold,
            ld_model=ld_arg_,
        )
        # extract not classified transcripts
        # first column in the rejected log
        self._transcripts_not_classified = self.__get_fst_col(cl_rej_log)

        if self.stop_at_chain_class:
            self.die("User requested to halt TOGA after chain features extraction", rc=0)

    def __get_proc_pseudogenes_track(self):
        """Create annotation of processed genes in query."""
        to_log("Creating processed pseudogenes track.")
        proc_pgenes_track = os.path.join(self.wd, "proc_pseudogenes.bed")
        create_ppgene_track(
            self.pred_scores, self.chain_file, self.index_bed_file, proc_pgenes_track
        )

    @staticmethod
    def __split_file(src, dst_dir, pieces_num):
        """Split file src into pieces_num pieces and save them to dst_dir."""
        f = open(src, "r")
        # skip header
        lines = list(f.readlines())[1:]
        f.close()
        lines_num = len(lines)
        piece_size = ceil(lines_num / pieces_num)
        paths_to_pieces = []
        for num, piece in enumerate(parts(lines, piece_size), 1):
            f_name = f"part_{num}"
            f_path = os.path.join(dst_dir, f_name)
            paths_to_pieces.append(f_path)
            f = open(f_path, "w")
            f.write("".join(piece))
            f.close()
        return paths_to_pieces
    
    def __get_transcript_to_strand(self):
        """Get """
        ret = {}
        f = open(self.ref_bed, 'r')
        for line in f:
            ld = line.rstrip().split("\t")
            trans = ld[3]
            direction = ld[5]
            ret[trans] = direction
        f.close()
        return ret

    def __get_chain_to_qstrand(self):
        ret = {}
        f = open(self.chain_file, "r")
        for line in f:
            if not line.startswith("chain"):
                continue
            fields = line.rstrip().split()
            chain_id = int(fields[-1])
            q_strand = fields[9]
            ret[chain_id] = q_strand
        f.close()
        return ret

    def __fold_exon_data(self, exons_data, out_bed):
        """Convert exon data into bed12."""
        projection_to_exons = defaultdict(list)
        # 1: make projection:
        f = open(exons_data, "r")
        for line in f:
            ld = line.rstrip().split("\t")
            transcript, _chain, _start, _end = ld
            chain = int(_chain)
            start = int(_start)
            end = int(_end)
            region = (start, end)
            projection_to_exons[(transcript, chain)].append(region)
        # 2 - get search loci for each projection
        projection_to_search_loc = {}

        # CESAR_PRECOMPUTED_ORTHO_LOCI_DATA = "cesar_precomputed_orthologous_loci.tsv"
        f = open(self.precomp_query_loci_path, "r")
        for elem in f:
            elem_data = elem.split("\t")
            transcript = elem_data[1]
            chain = int(elem_data[2])
            projection = f"{transcript}.{chain}"
            search_locus = elem_data[3]
            chrom, s_e = search_locus.split(":")
            s_e_split = s_e.split("-")
            start, end = int(s_e_split[0]), int(s_e_split[1])
            projection_to_search_loc[(transcript, chain)] = (chrom, start, end)
        f.close()
        trans_to_strand = self.__get_transcript_to_strand()
        chain_to_qstrand = self.__get_chain_to_qstrand()
        # 3 - save bed 12

        # WARNING DOING BED4
        f = open(out_bed, "w")
        f2 = open(f"{out_bed}4", "w")

        for (transcript, chain), exons in projection_to_exons.items():
            projection = f"{transcript}.{chain}"
            trans_strand = trans_to_strand[transcript]
            chain_strand = chain_to_qstrand[chain]
            to_invert = trans_strand != chain_strand

            exons_sort = sorted(exons, key=lambda x: x[0])
            if to_invert:
                exons_sort = exons_sort[::-1]

            # crash here
            search_locus = projection_to_search_loc[(transcript, chain)]
            chrom = search_locus[0]
            search_start = search_locus[1]
            search_end = search_locus[2]

            if to_invert:
                bed_start = search_end - exons_sort[0][1]
                bed_end = search_end - exons_sort[-1][0]
            else:
                bed_start = exons_sort[0][0] + search_start
                bed_end = exons_sort[-1][1] + search_start

            block_sizes, block_starts = [], []
            for exon in exons_sort:
                abs_start_in_s, abs_end_in_s = exon

                if to_invert:
                    abs_start = search_end - abs_end_in_s
                    abs_end = search_end - abs_start_in_s
                else:
                    abs_start = abs_start_in_s + search_start
                    abs_end = abs_end_in_s + search_start

                f2.write(f"{chrom}\t{abs_start}\t{abs_end}\t{projection}\n")

                rel_start = abs_start - bed_start
                block_size = abs_end - abs_start
                block_sizes.append(block_size)
                block_starts.append(rel_start)
            block_sizes_field = ",".join(map(str, block_sizes))
            block_starts_field = ",".join(map(str, block_starts))
            all_fields = (
                chrom,
                bed_start,
                bed_end,
                projection,
                0,
                0,
                bed_start,
                bed_end,
                "0,0,0",
                len(exons),
                block_sizes_field,
                block_starts_field,
            )
            f.write("\t".join(map(str, all_fields)))
            f.write("\n")
        f.close()
        f2.close()

    def __precompute_data_for_opt_cesar(self):
        """Precompute memory data for optimised CESAR.

        LASTZ optimised CESAR requires different (lesser) amounts of memory.
        To compute them we need to run another cluster-dependent step.
        """
        if not self.cesar_is_opt:
            return  # in case we use standard CESAR this procedure is not needed

        # need gene: chains
        # artificial bed file with all possible exons per gene, maybe the longest if intersection... or not?
        # what if there is a giant exon that does not align but smaller versions do?
        to_log("Computing memory requirements for optimized CESAR")
        to_log("Warning! This is an experimental feature")

        mem_dir = os.path.join(self.wd, CESAR_PRECOMPUTED_MEMORY_DIRNAME)
        self.precomp_reg_dir = os.path.join(self.wd, CESAR_PRECOMPUTED_REGIONS_DIRNAME)

        chain_class_temp_dir = os.path.join(self.wd, TEMP_CHAIN_CLASS)
        self.precomp_query_loci_dir = os.path.join(
            self.wd, CESAR_PRECOMPUTED_ORTHO_LOCI_DIRNAME
        )
        self.precomp_query_loci_path = os.path.join(
            self.temp_wd, CESAR_PRECOMPUTED_ORTHO_LOCI_DATA
        )

        os.mkdir(mem_dir) if not os.path.isdir(mem_dir) else None
        os.mkdir(chain_class_temp_dir) if not os.path.isdir(
            chain_class_temp_dir
        ) else None
        os.mkdir(self.precomp_reg_dir) if not os.path.isdir(
            self.precomp_reg_dir
        ) else None
        os.mkdir(self.precomp_query_loci_dir) if not os.path.isdir(
            self.precomp_query_loci_dir
        ) else None

        self.temp_files.append(mem_dir)
        self.temp_files.append(chain_class_temp_dir)
        self.temp_files.append(self.precomp_reg_dir)
        self.temp_files.append(self.precomp_query_loci_dir)

        # split file containing gene: orthologous/paralogous chains into 100 pieces
        class_pieces = self.__split_file(
            self.transcript_to_chain_classes, chain_class_temp_dir, NUM_CESAR_MEM_PRECOMP_JUBS
        )
        precompute_jobs = []
        # for nextflow abspath is better:
        precomp_abspath = os.path.abspath(self.PRECOMPUTE_OPT_CESAR_DATA)
        # cesar_bin_abspath = os.path.abspath(self.cesar_binary)

        # create joblist: one job per piece
        for num, class_piece in enumerate(class_pieces, 1):
            out_m_path = os.path.join(mem_dir, f"part_{num}")
            out_reg_path = os.path.join(self.precomp_reg_dir, f"part_{num}")
            out_ol_path = os.path.join(self.precomp_query_loci_dir, f"path_{num}")

            cmd = (
                f"{precomp_abspath} {class_piece} {self.wd} {self.opt_cesar_binary} "
                f"{self.t_2bit} {self.q_2bit} {out_m_path} --ro {out_reg_path} --ol {out_ol_path}"
            )
            precompute_jobs.append(cmd)
        # save joblist
        precomp_mem_joblist = os.path.join(
            self.wd, "_temp_cesar_precompute_memory_joblist.txt"
        )
        with open(precomp_mem_joblist, "w") as f:
            f.write("\n".join(precompute_jobs))
            f.write("\n")
        # now push these jobs to cluster
        timestamp = str(time.time()).split(".")[0]
        project_name = f"{self.project_name}_cesar_precompute_{timestamp}"
        project_path = os.path.join(self.nextflow_dir, project_name)
        to_log(f"Precomputing cesar regions, NF project name: {project_name}")
        to_log(f"Project path: {project_path}")

        if self.para:  # run jobs with para, skip nextflow
            cmd = f'para make {project_name} {precomp_mem_joblist} -q="shortmed"'
            to_log(f"Executing parallel step with {cmd}")
            rc = subprocess.call(cmd, shell=True)
        else:  # calling jobs with nextflow
            cmd = f"nextflow {self.NF_EXECUTE} " f"--joblist {precomp_mem_joblist}"
            if not self.local_executor:
                # not local executor -> provided config files
                # need abspath for nextflow execution
                # use the same nf config as for chain features, OK
                cmd += f" -c {self.nf_chain_extr_config_file}"
            to_log(f"Executing parallel step with {cmd}")
            os.mkdir(project_path) if not os.path.isdir(project_path) else None
            rc = subprocess.call(cmd, shell=True, cwd=project_path)

        if rc != 0:  # if process (para or nf) died: terminate execution
            self.die(f"Error! Process {cmd} died")
        if not self.keep_nf_logs and not self.para:
            # remove nextflow intermediate files
            # if para: this dir doesn't exist
            shutil.rmtree(project_path) if os.path.isdir(project_path) else None
        # merge files and quit
        self.__merge_dir(mem_dir, self.precomp_mem_cesar)
        self.cesar_mem_was_precomputed = True

        if self.output_opt_cesar_regions:
            save_to = os.path.join(self.wd, "precomp_regions_for_cesar.bed")
            exon_regions_df = os.path.join(self.wd, "precomp_regions_exons.tsv")
            self.__merge_dir(self.precomp_reg_dir, exon_regions_df)
            self.__merge_dir(self.precomp_query_loci_dir, self.precomp_query_loci_path)
            self.__fold_exon_data(exon_regions_df, save_to)
            to_log(f"Bed file containing precomputed regions is saved to: {save_to}")
            sys.exit(0)

    def __split_cesar_jobs(self):
        """Call split_exon_realign_jobs.py."""
        # TODO: check whether it is reachable
        if not self.t_2bit or not self.q_2bit:
            self.die(
                "There are no 2 bit files provided, cannot go ahead and call CESAR.",
                rc=0,
            )

        if self.fragmented_genome:
            to_log("Detecting fragmented transcripts")
            # need to stitch fragments together
            gene_fragments = stitch_scaffolds(
                self.chain_file, self.pred_scores, self.ref_bed, True
            )
            # for now split exon realign jobs is a subprocess, not a callable function
            # TODO: fix this
            fragm_dict_file = os.path.join(self.temp_wd, "gene_fragments.txt")
            f = open(fragm_dict_file, "w")
            for k, v in gene_fragments.items():
                v_str = ",".join(map(str, v))
                f.write(f"{k}\t{v_str}\n")
            f.close()
            to_log(f"Fragments data saved to {fragm_dict_file}")
        else:
            # no fragment file: ok
            to_log("Skip fragmented genes detection")
            fragm_dict_file = None

        # if we call CESAR
        to_log("Setting up creating CESAR jobs")
        self.cesar_jobs_dir = os.path.join(self.temp_wd, "cesar_jobs")
        self.cesar_combined = os.path.join(self.temp_wd, "cesar_combined")
        self.cesar_results = os.path.join(self.temp_wd, "cesar_results")
        self.cesar_bigmem_jobs = os.path.join(self.temp_wd, "cesar_bigmem")

        self.temp_files.append(self.cesar_jobs_dir)
        self.temp_files.append(self.cesar_combined)
        self.temp_files.append(self.predefined_glp_cesar_split)
        self.temp_files.append(self.technical_cesar_err)
        os.mkdir(self.technical_cesar_err) if not os.path.isdir(
            self.technical_cesar_err
        ) else None

        # different field names depending on --ml flag
        self.temp_files.append(self.cesar_results)
        self.temp_files.append(self.gene_loss_data)
        skipped_path = os.path.join(self.rejected_dir, "SPLIT_CESAR.txt")
        self.paralogs_log = os.path.join(self.temp_wd, "paralogs.txt")
        cesar_binary_to_use = (
            self.opt_cesar_binary if self.cesar_is_opt else self.cesar_binary
        )

        split_cesar_cmd = (
            f"{self.SPLIT_EXON_REALIGN_JOBS} "
            f"{self.transcript_to_chain_classes} {self.ref_bed} "
            f"{self.index_bed_file} {self.chain_index_file} "
            f"{self.t_2bit} "
            f"{self.q_2bit} "
            f"{self.wd} "
            f"--jobs_dir {self.cesar_jobs_dir} "
            f"--jobs_num {self.cesar_jobs_num} "
            f"--combined {self.cesar_combined} "
            f"--bigmem {self.cesar_bigmem_jobs} "  # TODO: refactor this part
            f"--results {self.cesar_results} "
            f"--buckets {self.cesar_buckets} "
            f"--mem_limit {self.cesar_mem_limit} "
            f"--chains_limit {self.cesar_chain_limit} "
            f"--skipped_genes {skipped_path} "
            f"--rejected_log {self.rejected_dir} "
            f"--cesar_binary {cesar_binary_to_use} "
            f"--paralogs_log {self.paralogs_log} "
            f"--uhq_flank {self.uhq_flank} "
            f"--predefined_glp_class_path {self.predefined_glp_cesar_split} "
            f"--unprocessed_log {self.technical_cesar_err} "
            f"--log_file {self.log_file}"
        )

        if self.annotate_paralogs:  # very rare occasion
            split_cesar_cmd += f" --annotate_paralogs"
        if self.mask_all_first_10p:
            split_cesar_cmd += f" --mask_all_first_10p"
        # split_cesar_cmd = split_cesar_cmd + f" --cesar_binary {self.cesar_binary}" \
        #     if self.cesar_binary else split_cesar_cmd
        split_cesar_cmd = (
            split_cesar_cmd + f" --opt_cesar" if self.cesar_is_opt else split_cesar_cmd
        )
        split_cesar_cmd = (
            split_cesar_cmd + " --mask_stops" if self.mask_stops else split_cesar_cmd
        )
        split_cesar_cmd = (
            split_cesar_cmd + f" --u12 {self.u12}" if self.u12 else split_cesar_cmd
        )
        split_cesar_cmd = (
            split_cesar_cmd + " --o2o_only" if self.o2o_only else split_cesar_cmd
        )
        split_cesar_cmd = (
            split_cesar_cmd + " --no_fpi" if self.no_fpi else split_cesar_cmd
        )
        if self.gene_loss_data:
            split_cesar_cmd += f" --check_loss {self.gene_loss_data}"
        if fragm_dict_file:
            split_cesar_cmd += f" --fragments_data {fragm_dict_file}"
        if self.cesar_mem_was_precomputed:
            split_cesar_cmd += f" --precomp_memory_data {self.precomp_mem_cesar}"
            split_cesar_cmd += f" --precomp_regions_data_dir {self.precomp_reg_dir}"
        self.__call_proc(split_cesar_cmd, "Could not split CESAR jobs!")

    def __get_cesar_jobs_for_bucket(self, comb_file, bucket_req):
        """Extract all cesar jobs belong to the bucket."""
        lines = []
        f = open(comb_file, "r")
        for line in f:
            line_data = line.split()
            if not line_data[0].endswith("cesar_runner.py"):
                # just a sanity check: each line is a command
                # calling cesar_runner.py
                self.die(f"CESAR joblist {comb_file} is corrupted!")
            # cesar job ID contains data of the bucket
            jobs = line_data[1]
            jobs_basename_data = os.path.basename(jobs).split("_")
            bucket = jobs_basename_data[-1]
            if bucket_req != bucket:
                continue
            lines.append(line)
        f.close()
        return "".join(lines)

    def __monitor_jobs(self, processes, project_paths, die_if_sc_1=False):
        """Monitor para/ nextflow jobs."""
        if self.exec_cesar_parts_sequentially is True:
            # no need to monitor jobs in this way
            # if they are called sequentially
            return
        iter_num = 0
        while True:  # Run until all jobs are done (or crashed)
            all_done = True  # default val, re-define if something is not done
            for p in processes:
                # check if each process is still running
                running = p.poll() is None
                if running:
                    all_done = False
            if all_done:
                to_log("### CESAR jobs done ###")
                break
            else:
                to_log(f"Status check iteration {iter_num} waiting {ITER_DURATION * iter_num} seconds.")
                time.sleep(ITER_DURATION)
                iter_num += 1
        if not self.keep_nf_logs:
            # remove nextflow intermediate files
            to_log("Removing nextflow temp files")
            for path in project_paths:
                shutil.rmtree(path) if os.path.isdir(path) else None

        if any(p.returncode != 0 for p in processes) and die_if_sc_1 is True:
            # some para/nextflow job died: critical issue
            # if die_if_sc_1 is True: terminate the program
            err = "Error! Some para/nextflow processes died!"
            self.die(err, 1)

    def __run_cesar_jobs(self):
        """Run CESAR jobs using nextflow.

        At first -> push joblists, there might be a few of them
        At second -> monitor joblists, wait until all are done.
        """
        # for each bucket I create a separate joblist and config file
        # different config files because different memory limits
        project_paths = []  # dirs with logs
        processes = []  # keep subprocess objects here
        timestamp = str(time.time()).split(".")[1]  # for project name
        project_names = []

        # get a list of buckets
        if self.cesar_buckets == "0":
            buckets = [
                0,
            ]  # a single bucket
        else:  # several buckets, each int -> memory limit in gb
            buckets = [int(x) for x in self.cesar_buckets.split(",") if x != ""]
        to_log(f"Pushing {len(buckets)} CESAR job lists")

        # generate buckets, create subprocess instances
        for b in buckets:
            # create config file
            # 0 means that that buckets were not split
            mem_lim = b if b != 0 else self.cesar_mem_limit

            if not self.local_executor:
                # running on cluster, need to create config file
                # for this bucket's memory requirement
                config_string = self.cesar_config_template.replace(
                    "${_MEMORY_}", f"{mem_lim}"
                )
                config_file_path = os.path.join(
                    self.temp_wd, f"cesar_config_{b}_queue.nf"
                )
                config_file_abspath = os.path.abspath(config_file_path)
                with open(config_file_path, "w") as f:
                    f.write(config_string)
                self.temp_files.append(config_file_path)
            else:  # no config dir given: use local executor
                # OK if there is a single bucket
                config_file_abspath = None

            if b != 0:  # extract jobs related to this bucket (if it's not 0)
                bucket_tasks = self.__get_cesar_jobs_for_bucket(
                    self.cesar_combined, str(b)
                )
                if len(bucket_tasks) == 0:
                    to_log(f"There are no jobs in the {b}Gb bucket")
                    continue
                joblist_name = f"cesar_joblist_queue_{b}.txt"
                joblist_path = os.path.join(self.temp_wd, joblist_name)
                with open(joblist_path, "w") as f:
                    f.write(bucket_tasks)
                joblist_abspath = os.path.abspath(joblist_path)
                self.temp_files.append(joblist_path)
            else:  # nothing to extract, there is a single joblist
                joblist_abspath = os.path.abspath(self.cesar_combined)

            # create project directory for logs
            nf_project_name = f"{self.project_name}_cesar_at_{timestamp}_q_{b}"
            project_names.append(nf_project_name)
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            project_paths.append(nf_project_path)

            # create subprocess object
            if not self.para:  # create nextflow cmd
                os.mkdir(nf_project_path) if not os.path.isdir(
                    nf_project_path
                ) else None

                cmd = f"nextflow {self.NF_EXECUTE} " f"--joblist {joblist_abspath}"
                if config_file_abspath:
                    cmd += f" -c {config_file_abspath}"
                p = subprocess.Popen(cmd, shell=True, cwd=nf_project_path)
            else:  # create cmd for para
                memory_mb = b * 1000
                cmd = f'para make {nf_project_name} {joblist_abspath} -q="shortmed"'
                if memory_mb > 0:
                    cmd += f" --memoryMb={memory_mb}"
                p = subprocess.Popen(cmd, shell=True)

            sys.stderr.write(f"Pushed cluster jobs with {cmd}\n")

            # wait for the process if calling processes sequentially
            # or wait a minute before pushing the next batch in parallel
            if self.exec_cesar_parts_sequentially is True:
                p.wait()
            else:
                time.sleep(CESAR_PUSH_INTERVAL)
            processes.append(p)

        # push bigmem jobs
        if self.nextflow_bigmem_config and not self.para:
            # if provided: push bigmem jobs also
            nf_project_name = f"{self.project_name}_cesar_at_{timestamp}_q_bigmem"
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            nf_cmd = (
                f"nextflow {self.NF_EXECUTE} "
                f"--joblist {self.cesar_bigmem_jobs} -c {self.nextflow_bigmem_config}"
            )
            # if bigmem joblist is empty or not exist: do nothing
            is_file = os.path.isfile(self.cesar_bigmem_jobs)
            if is_file:  # if a file: we can open and count lines
                f = open(self.cesar_bigmem_jobs, "r")
                big_lines_num = len(f.readlines())
                f.close()
            else:  # file doesn't exist, equivalent to an empty file in our case
                big_lines_num = 0
            # if it's empty: do nothing
            if big_lines_num == 0:
                pass
            else:
                os.mkdir(nf_project_path) if not os.path.isdir(
                    nf_project_path
                ) else None
                project_paths.append(nf_project_path)
                p = subprocess.Popen(nf_cmd, shell=True, cwd=nf_project_path)
                sys.stderr.write(f"Pushed {big_lines_num} bigmem jobs with {nf_cmd}\n")
                if self.exec_cesar_parts_sequentially:
                    p.wait()
                processes.append(p)
        elif self.para and self.para_bigmem:
            # if requested: push bigmem jobs with para
            is_file = os.path.isfile(self.cesar_bigmem_jobs)
            bm_project_name = f"{self.project_name}_cesar_at_{timestamp}_q_bigmem"
            if is_file:  # if a file: we can open and count lines
                f = open(self.cesar_bigmem_jobs, "r")
                big_lines_num = len(f.readlines())
                f.close()
            else:  # file doesn't exist, equivalent to an empty file in our case
                big_lines_num = 0
            # if it's empty: do nothing
            if big_lines_num == 0:
                pass
            else:  # there ARE cesar bigmem jobs: push them
                memory_mb = 500 * 1000  # TODO: bigmem memory param
                cmd = (
                    f"para make {bm_project_name} {self.cesar_bigmem_jobs} "
                    f'-q="bigmem" --memoryMb={memory_mb}'
                )
                p = subprocess.Popen(cmd, shell=True)
                processes.append(p)

        # monitor jobs
        self.__monitor_jobs(processes, project_paths)

        # print CPU runtime (if para), if not -> quit function
        if not self.para:
            return

        for p_name in project_names:
            cmd = f"para time {p_name}"
            p = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout_, stderr_ = p.communicate()
            stdout = stdout_.decode("utf-8")
            stderr = stderr_.decode("utf-8")
            to_log(f"para time output for {p_name}:")
            to_log(stdout)
            to_log(stderr)
            cmd_cleanup = f"para clean {p_name}"
            p = subprocess.Popen(
                cmd_cleanup, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            _, _ = p.communicate()

    @staticmethod
    def __get_bucket_val(mem_val, buckets):
        for b in buckets:
            if b > mem_val:
                return b
        return None

    @staticmethod
    def __read_chain_arg(chains):
        """Read chain argument.

        Return None if argument is invalid."""
        try:
            chain_ids = [int(x) for x in chains.split(",") if x != ""]
            return chain_ids
        except ValueError:
            return None

    def __rebuild_crashed_jobs(self, crashed_jobs):
        """If TOGA has to re-run CESAR jobs we still need some buckets."""
        bucket_to_jobs = defaultdict(list)
        buckets = [int(x) for x in self.cesar_buckets.split(",") if x != ""]
        if len(buckets) == 0:
            buckets.append(self.cesar_mem_limit)

        for elem in crashed_jobs:
            elem_args = elem.split()
            chains_arg = elem_args[2]
            chain_ids = self.__read_chain_arg(chains_arg)
            if chain_ids is None:
                continue
            try:
                memlim_arg_ind = elem_args.index(MEMLIM_ARG) + 1
                mem_val = float(elem_args[memlim_arg_ind])
            except ValueError:
                mem_val = self.cesar_mem_limit
            bucket_lim = self.__get_bucket_val(mem_val, buckets)

            if FRAGM_ARG in elem_args:
                cmd_copy = elem_args.copy()
                cmd_str = " ".join(cmd_copy)
                bucket_to_jobs[bucket_lim].append(cmd_str)
                continue

            for chain_id in chain_ids:
                cmd_copy = elem_args.copy()
                cmd_copy[2] = str(chain_id)
                cmd_str = " ".join(cmd_copy)
                bucket_to_jobs[bucket_lim].append(cmd_str)
        return bucket_to_jobs

    def __check_cesar_completeness(self):
        """Check that all CESAR jobs were executed, quit otherwise."""
        rejected_logs_filenames = [
            x for x in os.listdir(self.rejected_dir) if x.startswith("cesar")
        ]
        if len(rejected_logs_filenames) == 0:
            # nothing crashed
            to_log("No CESAR jobs crashed")
            return
        # collect crashed jobs
        crashed_jobs = []
        for filename in rejected_logs_filenames:
            path = os.path.join(self.rejected_dir, filename)
            f = open(path, "r")
            crashed_jobs_lines = [x for x in f.readlines() if "CESAR JOB FAILURE" in x]
            crashed_jobs_in_file = [x.split("\t")[0] for x in crashed_jobs_lines]
            crashed_jobs.extend(crashed_jobs_in_file)
        crashed_jobs_num = len(crashed_jobs)
        if crashed_jobs_num == 0:
            to_log("No CESAR jobs crashed")
            return
        elif crashed_jobs_num > 1000:  # this is TOO MUCH, must never happen
            err_msg_ = (
                f"Too many ({crashed_jobs_num}) CESAR jobs died, "
                f"please check your input data and re-run TOGA"
            )
            self.die(err_msg_)

        to_log(f"{crashed_jobs_num} CESAR jobs crashed, trying to run again...")
        bucket_to_jobs = self.__rebuild_crashed_jobs(crashed_jobs)
        temp_jobs_dir = os.path.join(self.wd, "RERUN_CESAR_JOBS")
        os.mkdir(temp_jobs_dir) if not os.path.isdir(temp_jobs_dir) else None
        timestamp = str(time.time()).split(".")[1]  # for project name
        p_objects = []
        project_paths = []
        self.rejected_dir_rerun = os.path.join(self.temp_wd, "cesar_jobs_crashed_again")
        os.mkdir(self.rejected_dir_rerun) if not os.path.isdir(
            self.rejected_dir_rerun
        ) else None
        err_log_files = []

        for bucket, jobs in bucket_to_jobs.items():
            to_log(f"Pushing {len(jobs)} jobs into {bucket} GB queue")
            batch_path = f"_cesar_rerun_batch_{bucket}"
            bucket_batch_file = os.path.join(self.wd, batch_path)
            batch_commands = []
            for num, job in enumerate(jobs, 1):
                job_file = f"rerun_job_{num}_{bucket}"
                job_path = os.path.join(temp_jobs_dir, job_file)
                # job_files.append(job_path)
                f = open(job_path, "w")
                f.write(job)
                f.write("\n")
                f.close()
                out_filename = f"rerun_job_{num}_{bucket}.txt"
                output_path = os.path.join(self.cesar_results, out_filename)
                inact_path = os.path.join(self.gene_loss_data, out_filename)
                rejected = os.path.join(self.rejected_dir_rerun, out_filename)
                err_log_files.append(rejected)
                batch_cmd = CESAR_RUNNER_TMP.format(
                    CESAR_RUNNER, job_path, output_path, inact_path, rejected
                )
                batch_commands.append(batch_cmd)
            f = open(bucket_batch_file, "w")
            f.write("\n".join(batch_commands))
            f.write("\n")
            f.close()

            nf_project_name = (
                f"{self.project_name}_RERUN_cesar_at_{timestamp}_q_{bucket}"
            )
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            project_paths.append(nf_project_path)
            if self.para:
                # push this using para
                memory_mb = bucket * 1000
                cmd = f'para make {nf_project_name} {bucket_batch_file} -q="shortmed"'
                if memory_mb > 0:
                    cmd += f" --memoryMb={memory_mb}"
                p = subprocess.Popen(cmd, shell=True)
            else:
                # push jobs using nextflow
                os.mkdir(nf_project_path) if not os.path.isdir(
                    nf_project_path
                ) else None
                config_file_path = os.path.join(
                    self.wd, f"cesar_config_{bucket}_queue.nf"
                )
                config_file_abspath = os.path.abspath(config_file_path)
                cmd = f"nextflow {self.NF_EXECUTE} " f"--joblist {bucket_batch_file}"
                if os.path.isfile(config_file_abspath):
                    cmd += f" -c {config_file_abspath}"
                p = subprocess.Popen(cmd, shell=True, cwd=nf_project_path)
            p_objects.append(p)
            time.sleep(CESAR_PUSH_INTERVAL)
        self.__monitor_jobs(p_objects, project_paths, die_if_sc_1=True)
        # TODO: maybe some extra sanity check

        # need to check whether anything crashed again
        crashed_twice = []
        for elem in err_log_files:
            f = open(elem, "r")
            lines = [x for x in f.readlines() if len(x) > 1]
            f.close()
            crashed_twice.extend(lines)
        shutil.rmtree(self.rejected_dir_rerun)
        if len(crashed_twice) == 0:
            to_log("All CESAR jobs re-ran succesfully!")
            return
        # OK, some jobs crashed twice
        crashed_log = os.path.join(self.wd, "cesar_jobs_crashed.txt")
        f = open(crashed_log, "w")
        f.write("".join(crashed_twice))
        f.close()
        err_msg = f"Some CESAR jobs crashed twice, please check {crashed_log}; Abort"
        self.die(err_msg, 1)

    def __append_technicall_err_to_predef_class(self, transcripts_path, out_path):
        """Append file with predefined clasifications."""
        if not os.path.isfile(transcripts_path):
            # in this case, we don't have transcripts with tech error
            # can simply quit the function
            return
        with open(transcripts_path, "r") as f:
            transcripts_list = [x.rstrip() for x in f]
        f = open(out_path, "a")
        for elem in transcripts_list:
            f.write(f"TRANSCRIPT\t{elem}\tM\n")
        f.close()

    def __merge_cesar_output(self):
        """Merge CESAR output, save final fasta and bed."""
        to_log("Merging CESAR output to make fasta and bed files.")
        merge_c_stage_skipped = os.path.join(self.rejected_dir, "CESAR_MERGE.txt")
        self.temp_files.append(self.intermediate_bed)

        all_ok = merge_cesar_output(
            self.cesar_results,
            self.intermediate_bed,
            self.nucl_fasta,
            self.meta_data,
            merge_c_stage_skipped,
            self.prot_fasta,
            self.codon_fasta,
            self.trash_exons,
            fragm_data=self.bed_fragm_exons_data,
        )

        # need to merge files containing transcripts that were not processed
        # for some technical reason, such as intersecting fragments in the query
        self.__merge_dir(
            self.technical_cesar_err, self.technical_cesar_err_merged, ignore_empty=True
        )

        self.__append_technicall_err_to_predef_class(
            self.technical_cesar_err_merged, self.predefined_glp_cesar_split
        )

        if len(all_ok) == 0:
            # there are no empty output files -> parsed without errors
            to_log("CESAR results merged")
            self.cesar_ok_merged = True
        else:
            # there are some empty output files
            # MAYBE everything is fine
            f = open(self.cesar_crashed_batches_log, "w")
            for err in all_ok:
                _path = err[0]
                _reason = err[1]
                f.write(f"{_path}\t{_reason}\n")
            f.close()
            # but need to notify user anyway

            warning_message = (
                f"WARNING!\nSOME CESAR JOB BATCHES LIKELY CRASHED\n!"
                f"RESULTS ARE LIKELY INCOMPLETE"
                f"PLEASE SEE {self.cesar_crashed_batches_log} FOR DETAILS"
                f"THERE IS A MINOR CHANCE THAT EVERYTHING IS CORRECT"
            )
            to_log(warning_message)
            self.cesar_ok_merged = False

    def __transcript_quality(self):
        """Call module to get transcript quality."""
        self.trans_quality_file = os.path.join(self.temp_wd, "transcript_quality.tsv")
        classify_transcripts(
            self.meta_data,
            self.pred_scores,
            self.hq_orth_threshold,
            self.trans_quality_file,
        )

    def __gene_loss_summary(self):
        """Call gene loss summary."""
        to_log("Calling gene loss summary")
        # collect already known cases
        predef_glp_classes = collect_predefined_glp_cases(
            self.predefined_glp_cesar_split
        )
        predef_missing = add_transcripts_to_missing(
            self._transcripts_not_intersected, self._transcripts_not_classified
        )
        predef_glp_classes.extend(predef_missing)

        # call classifier itself
        gene_losses_summary(
            self.gene_loss_data,
            self.ref_bed,
            self.intermediate_bed,
            self.query_annotation,
            self.loss_summ,
            iforms_file=self.isoforms,
            paral=self.paralogs_log,
            predefined_class=predef_glp_classes,
        )

    def __orthology_type_map(self):
        """Call orthology_type_map.py"""
        # need to combine projections in genes
        query_isoforms_file = os.path.join(self.wd, "query_isoforms.tsv")
        query_gene_spans = os.path.join(self.wd, "query_gene_spans.bed")
        get_query_isoforms_data(
            self.query_annotation,
            query_isoforms_file,
            save_genes_track=query_gene_spans,
        )
        to_log("Calling orthology types mapping step...")
        skipped_ref_trans = os.path.join(self.wd, "ref_orphan_transcripts.txt")
        orthology_type_map(
            self.ref_bed,
            self.query_annotation,
            self.orthology_type,
            ref_iso=self.isoforms,
            que_iso=query_isoforms_file,
            paralogs_arg=self.paralogs_log,
            loss_data=self.loss_summ,
            save_skipped=skipped_ref_trans,
            orth_scores_arg=self.pred_scores,
        )

    @staticmethod
    def __merge_dir(dir_name, output, ignore_empty=False):
        """Merge all files in a directory into one."""
        files_list = os.listdir(dir_name)
        if len(files_list) == 0 and ignore_empty is False:
            sys.exit(f"Error! {dir_name} is empty")
        elif len(files_list) == 0 and ignore_empty is True:
            # in this case we allow empty directories
            # just remove the directory and return
            shutil.rmtree(dir_name)
            return
        buffer = open(output, "w")
        for filename in files_list:
            path = os.path.join(dir_name, filename)
            with open(path, "r") as f:
                # content = [x for x in f.readlines() if x != "\n"]
                content = f.read()
            # lines = "".join(content)
            # buffer.write(lines)
            buffer.write(content)
        buffer.close()
        # shutil.rmtree(dir_name)

    def __check_crashed_cesar_jobs(self):
        """Check whether any CESAR jobs crashed.

        Previously we checked entire batches of CESAR jobs.
        Now we are seeking for individual jobs.
        Most likely they crashed due to internal CESAR error.
        """
        # grep for crashed jobs in the rejected logs
        if not os.path.isfile(self.rejected_log):
            return  # no log: nothing to do
        f = open(self.rejected_log, "r")
        for line in f:
            if "CESAR" not in line:
                # not related to CESAR
                continue
            if "fragment chains oevrlap" in line:
                # they are not really "crashed"
                # those jobs could not produce a meaningful result
                # only intersecting exons
                continue
            # extract cesar wrapper command from the log
            cesar_cmd = line.split("\t")[0]
            self.crashed_cesar_jobs.append(cesar_cmd)
        f.close()
        # save crashed jobs list if it's non-empty
        if len(self.crashed_cesar_jobs) == 0:
            return  # good, nothing crashed
        to_log(f"{len(self.crashed_cesar_jobs)} CESAR wrapper commands failed")
        f = open(self.cesar_crashed_jobs_log, "w")
        for cmd in self.crashed_cesar_jobs:
            f.write(f"{cmd}\n")
        to_log(f"Failed CESAR wrapper commands were written to: {self.cesar_crashed_jobs_log}")
        f.close()

    def __left_done_mark(self):
        """Write a file confirming that everything is done."""
        mark_file = os.path.join(self.wd, "done.status")
        start_mark = os.path.join(self.wd, RUNNING)
        f = open(mark_file, "w")
        now_ = str(dt.now())
        f.write(f"Done at {now_}\n")
        if not self.cesar_ok_merged:
            f.write("\n:Some CESAR batches produced an empty result, please see:\n")
            f.write(f"{self.cesar_crashed_batches_log}\n")
        if len(self.crashed_cesar_jobs) > 0:
            num_ = len(self.crashed_cesar_jobs)
            f.write(f"Some individual CESAR jobs ({num_} jobs) crashed\n")
            f.write(f"Please see:\n{self.cesar_crashed_jobs_log}\n")
        f.close()
        os.remove(start_mark) if os.path.isfile(start_mark) else None

    def __get_version(self):
        """Get git hash and current branch if possible."""
        cmd_hash = "git rev-parse HEAD"
        cmd_branch = "git rev-parse --abbrev-ref HEAD"
        try:
            git_hash = subprocess.check_output(
                cmd_hash, shell=True, cwd=self.toga_exe_path
            ).decode("utf-8").strip()
            git_branch = subprocess.check_output(
                cmd_branch, shell=True, cwd=self.toga_exe_path
            ).decode("utf-8").strip()
        except subprocess.CalledProcessError:
            git_hash = "unknown"
            git_branch = "unknown"
        version = f"Version {__version__}\nCommit: {git_hash}\nBranch: {git_branch}\n"
        to_log(version)
        return version

    def __merge_split_files(self):
        """Merge intermediate/temp files."""
        # merge rejection logs
        self.__merge_dir(self.rejected_dir, self.rejected_log)
        # save inact mutations data
        inact_mut_file = os.path.join(self.wd, "inact_mut_data.txt")
        self.__merge_dir(self.gene_loss_data, inact_mut_file)
        # save CESAR outputs
        cesar_results_merged = os.path.join(self.temp_wd, "cesar_results.txt")
        self.__merge_dir(self.cesar_results, cesar_results_merged)
        # merge CESAR jobs
        cesar_jobs_merged = os.path.join(self.temp_wd, "cesar_jobs_merged")
        self.__merge_dir(self.cesar_jobs_dir, cesar_jobs_merged)
        # remove chain classification data
        shutil.rmtree(self.ch_cl_jobs) if os.path.isdir(self.ch_cl_jobs) else None
        shutil.rmtree(self.chain_class_results) if os.path.isdir(
            self.chain_class_results
        ) else None


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument(
        "chain_input",
        type=str,
        help="Chain file. Extensions like "
        "FILE.chain or FILE.chain.gz are applicable.",
    )
    app.add_argument(
        "bed_input", type=str, help="Bed file with annotations for the target genome."
    )
    app.add_argument(
        "tDB", default=None, help="Reference genome sequence in 2bit format."
    )
    app.add_argument("qDB", default=None, help="Query genome sequence in 2bit format.")
    # global ops
    app.add_argument(
        "--project_dir",
        "--pd",
        default=None,
        help="Project directory. TOGA will save all intermediate and output files "
        "exactly in this directory. If not specified, use CURRENT_DIR/PROJECT_NAME "
        "as default (see below).",
    )
    app.add_argument(
        "--project_name",
        "--pn",
        default=None,
        help="If you don't like to provide a full path to the project directory with "
        "--project_dir you can use this parameter. In this case TOGA will "
        "create project directory in the current directory as "
        '"CURRENT_DIR/PROJECT_NAME". If not provided, TOGA will try to extract '
        "the project name from chain filename, which is not recommended.",
    )
    app.add_argument(
        "--min_score",
        "--msc",
        type=int,
        default=15000,
        help="Chain score threshold. Exclude chains that have a lower score "
        "from the analysis. Default value is 15000.",
    )
    app.add_argument(
        "--isoforms", "-i", type=str, default="", help="Path to isoforms data file"
    )
    app.add_argument(
        "--keep_temp",
        "--kt",
        action="store_true",
        dest="keep_temp",
        help="Do not remove temp files.",
    )
    app.add_argument(
        "--limit_to_ref_chrom",
        default=None,
        help="Find orthologs " "for a single reference chromosome only",
    )
    app.add_argument(
        "--nextflow_dir",
        "--nd",
        default=None,
        help="Nextflow working directory: from this directory nextflow is "
        "executed, also there all nextflow log files are kept",
    )
    app.add_argument(
        "--nextflow_config_dir",
        "--nc",
        default=None,
        help="Directory containing nextflow configuration files "
        "for cluster, pls see nextflow_config_files/readme.txt "
        "for details.",
    )
    app.add_argument(
        "--do_not_del_nf_logs", "--nfnd", action="store_true", dest="do_not_del_nf_logs"
    )
    app.add_argument(
        "--cesar_bigmem_config",
        "--nb",
        default=None,
        help="File containing "
        "nextflow config for BIGMEM CESAR jobs. If not provided, these "
        "jobs will not run (but list of them saved)/ NOT IMPLEMENTED YET",
    )
    app.add_argument(
        "--para",
        "-p",
        action="store_true",
        dest="para",
        help="Hillerlab feature, use para instead of nextflow to "
        "manage cluster jobs.",
    )
    app.add_argument(
        "--para_bigmem",
        "--pb",
        action="store_true",
        dest="para_bigmem",
        help="Hillerlab feature, push bigmem jobs with para",
    )
    # chain features related
    app.add_argument(
        "--chain_jobs_num",
        "--chn",
        type=int,
        default=100,
        help="Number of cluster jobs for extracting chain features. "
        "Recommended from 150 to 200 jobs.",
    )
    app.add_argument(
        "--no_chain_filter",
        "--ncf",
        action="store_true",
        dest="no_chain_filter",
        help="A flag. Do not filter the chain file (make sure you specified "
        "a .chain but not .gz file in this case)",
    )
    app.add_argument(
        "--orth_score_threshold",
        "--ost",
        default=0.5,
        type=float,
        help="Score threshold to distinguish orthologs from paralogs, default 0.5",
    )
    # CESAR part related
    app.add_argument(
        "--cesar_jobs_num",
        "--cjn",
        type=int,
        default=500,
        help="Number of CESAR cluster jobs.",
    )
    app.add_argument(
        "--cesar_binary", default=None, help="CESAR binary address, cesar as default."
    )
    app.add_argument(
        "--using_optimized_cesar",
        "--uoc",
        action="store_true",
        dest="using_optimized_cesar",
        help="Instead of CESAR, use lastz-based optimized version",
    )
    app.add_argument(
        "--output_opt_cesar_regions",
        "--oocr",
        action="store_true",
        dest="output_opt_cesar_regions",
        help="If optimised CESAR is used, save the included regions, "
        "and quit. The parameter is defined for debugging purposes only.",
    )
    app.add_argument(
        "--mask_stops",
        "--ms",
        action="store_true",
        dest="mask_stops",
        help="Mask stop codons in target sequences. "
        "CESAR cannot process them. Using this "
        "parameter please make sure you know what you are doing.",
    )
    app.add_argument(
        "--cesar_buckets",
        "--cb",
        default="0",
        help="Comma-separated list of integers. Split CESAR jobs in buckets "
        "depending on their memory requirements. "
        "See README.md for explanation.",
    )
    app.add_argument(
        "--cesar_exec_seq",
        "--ces",
        action="store_true",
        dest="cesar_exec_seq",
        help="Execute different CESAR jobs partitions sequentially, "
        "not in parallel.",
    )
    app.add_argument(
        "--cesar_chain_limit",
        type=int,
        default=100,
        help="Skip genes that have more that CESAR_CHAIN_LIMIT orthologous "
        "chains. Recommended values are a 50-100.",
    )
    app.add_argument(
        "--cesar_mem_limit",
        type=int,
        default=16,
        help="Ignore genes requiring > N gig to run CESAR",
    )
    app.add_argument(
        "--time_marks",
        "-t",
        default=None,
        help="File to save timings of different steps.",
    )
    app.add_argument("--u12", default=None, help="U12 introns data")
    app.add_argument(
        "--stop_at_chain_class",
        "--sac",
        action="store_true",
        dest="stop_at_chain_class",
        help="Stop after merging chain features.",
    )
    app.add_argument(
        "--uhq_flank", default=50, type=int, help="Flank size for UHQ exons"
    )
    app.add_argument(
        "--o2o_only",
        "--o2o",
        action="store_true",
        dest="o2o_only",
        help="Process only the genes that have a single orthologous chain",
    )
    app.add_argument(
        "--no_fpi",
        action="store_true",
        dest="no_fpi",
        help="Consider some frame-preserving mutations as inactivating. "
        "See documentation for details.",
    )
    app.add_argument(
        "--disable_fragments_joining",
        "--dfj",
        dest="disable_fragments_joining",
        action="store_true",
        help="Disable assembling query genes from pieces",
    )
    app.add_argument(
        "--ld_model",
        dest="ld_model",
        action="store_true",
        help="Apply extra classifier for molecular distances ~1sps.",
    )
    app.add_argument(
        "--annotate_paralogs",
        "--ap",
        dest="annotate_paralogs",
        action="store_true",
        help="Annotate paralogs instead of orthologs. "
        "(experimental feature for very specific needs)",
    )
    app.add_argument(
        "--mask_all_first_10p",
        "--m_f10p",
        action="store_true",
        dest="mask_all_first_10p",
        help="Automatically mask all inactivating mutations in first 10 percent of "
        "the reading frame, ignoring ATG codons distribution. "
        "(Default mode in V1.0, not recommended to use in later versions)",
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    # some sanity checks
    if args.output_opt_cesar_regions and not args.using_optimized_cesar:
        err_msg = (
            "Error! Please use --output_opt_cesar_regions parameter "
            " with --using_optimized_cesar flag only"
        )
        sys.exit(err_msg)
    return args


def main():
    """Entry point."""
    args = parse_args()
    toga_manager = Toga(args)
    toga_manager.run()


if __name__ == "__main__":
    main()
