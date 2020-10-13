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
import functools
from twobitreader import TwoBitFile
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
# from modules.common import eprint
from modules.stitch_fragments import stitch_scaffolds


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

U12_FILE_COLS = 3
U12_AD_FIELD = {"A", "D"}
ISOFORMS_FILE_COLS = 2
NF_DIR_NAME = "nextflow_logs"
CESAR_PUSH_INTERVAL = 30  # CESAR jobs push interval
ITER_DURATION = 60  # CESAR jobs check interval
# automatically enable flush
print = functools.partial(print, flush=True)


class Toga:
    """TOGA manager class."""
    def __init__(self, args):
        """Initiate toga class."""
        self.t0 = dt.now()
        # check if all files TOGA needs are here
        self.temp_files = []  # remove at the end, list of temp files
        print("#### Initiating TOGA class ####")
        print("Checking dependencies...")
        self.__modules_addr()
        self.__check_dependencies()
        self.__check_completeness()
        self.para = args.para
        self.toga_exe_path = os.path.dirname(__file__)
        self.version = self.__get_version()
        self.para_bigmem = args.para_bigmem
        self.nextflow_dir = self.__get_nf_dir(args.nextflow_dir)
        self.nextflow_config_dir = args.nextflow_config_dir
        if args.cesar_bigmem_config:
            self.nextflow_bigmem_config = os.path.abspath(args.cesar_bigmem_config)
        else:  # if none: we cannot call os.path.abspath method
            self.nextflow_bigmem_config = None
        self.__check_nf_config()
        # to avoid crash on filesystem without locks:
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash

        chain_basename = os.path.basename(args.chain_input)

        # create project dir
        self.project_name = self.__gen_project_name() if not args.project_name \
            else args.project_name
        self.wd = os.path.abspath(args.project_dir) if args.project_dir else  \
            os.path.join(os.getcwd(), self.project_name)
        # for safety; need this to make paths later
        self.project_name = self.project_name.replace("/", "")
        os.mkdir(self.wd) if not os.path.isdir(self.wd) else None
        print(f"Output directory: {self.wd}")

        # dir to collect log files with rejected reference genes:
        self.rejected_dir = os.path.join(self.wd, "rejected")
        os.mkdir(self.rejected_dir) if not os.path.isdir(self.rejected_dir) else None

        # filter chain in this folder
        g_ali_basename = "genome_alignment"
        self.chain_file = os.path.join(self.wd, f"{g_ali_basename}.chain")
        # there is an assumption that chain file has .chain extension
        # chain indexing was a bit problematic: (i) bsddb3 fits perfectly but is very
        # painful to install, (ii) sqlite is also fine but might be dysfunctional on some
        # cluster file systems, so we create chain_ID: (start_byte, offset) dictionary for
        # instant extraction of a particular chain from the chain file
        # we save these dictionaries into two files: a text file (tsv) and binary file with BST
        # depending on the case we will use both (for maximal performance)
        self.chain_index_file = os.path.join(self.wd, f"{g_ali_basename}.bst")
        self.chain_index_txt_file = os.path.join(self.wd, f"{g_ali_basename}.chain_ID_position")

        # make the command, prepare the chain file
        if not os.path.isfile(args.chain_input):
            chain_filter_cmd = None
            self.die(f"Error! File {args.chain_input} doesn't exist!")
        elif chain_basename.endswith(".gz"):  # version for gz
            chain_filter_cmd = f"gzip -dc {args.chain_input} | "\
                               f"{self.CHAIN_SCORE_FILTER} stdin "\
                               f"{args.min_score} > {self.chain_file}"
        elif args.no_chain_filter:  # it is .chain and score filter is not required
            chain_filter_cmd = f"rsync -a {args.chain_input} {self.chain_file}"
        else:  # it is .chain | score filter required
            chain_filter_cmd = f"{self.CHAIN_SCORE_FILTER} {args.chain_input} "\
                               f"{args.min_score} > {self.chain_file}"

        # filter chains with score < threshold
        self.__call_proc(chain_filter_cmd, "Please check if you use a proper chain file.")
        
        # bed define bed files addresses
        self.ref_bed = os.path.join(self.wd, "toga_filt_ref_annot.bed")
        self.index_bed_file = os.path.join(self.wd, "toga_filt_ref_annot.hdf5")

        # filter bed file
        bed_filt_rejected_file = "BED_FILTER_REJECTED.txt"
        bed_filt_rejected = os.path.join(self.rejected_dir,
                                         bed_filt_rejected_file)
        # keeping UTRs!
        prepare_bed_file(args.bed_input,
                         self.ref_bed,
                         save_rejected=bed_filt_rejected,
                         only_chrom=args.limit_to_ref_chrom)

        # mics things
        self.isoforms_arg = args.isoforms if args.isoforms else None
        self.isoforms = None  # will be assigned after completeness check
        self.chain_jobs = args.chain_jobs_num
        self.cesar_binary = self.DEFAULT_CESAR if not args.cesar_binary \
            else args.cesar_binary
        self.time_log = args.time_marks
        self.stop_at_chain_class = args.stop_at_chain_class

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
        self.cesar_fields = args.homology_types
        self.mask_stops = args.mask_stops
        self.no_fpi = args.no_fpi
        self.o2o_only = args.o2o_only
        self.keep_nf_logs = args.do_not_del_nf_logs
        self.cesar_ok_merged = None
        self.fragmented_genome = args.fragmented_genome

        self.chain_results_df = os.path.join(self.wd, "chain_results_df.tsv")
        self.nucl_fasta = os.path.join(self.wd, "nucleotide.fasta")
        self.prot_fasta = os.path.join(self.wd, "prot.fasta")
        self.codon_fasta = os.path.join(self.wd, "codon.fasta")
        self.low_conf_bed = os.path.join(self.wd, "low_confidence.bed")
        self.meta_data = os.path.join(self.wd, "exons_meta_data.tsv")
        self.intermediate_bed = os.path.join(self.wd, "intermediate.bed")
        self.orthology_type = os.path.join(self.wd, "orthology_classification.tsv")
        self.classification_log = os.path.join(self.wd, "o_class.log")
        self.trash_exons = os.path.join(self.wd, "trash_exons.bed")
        self.gene_loss_data = os.path.join(self.wd, "inact_mut_data")
        self.query_annotation = os.path.join(self.wd, "query_annotation.bed")
        self.loss_summ = os.path.join(self.wd, "loss_summ_data.tsv")
        self.bed_fragm_exons_data = os.path.join(self.wd, "bed_fragments_to_exons.tsv")
        self.u12_arg = args.u12
        self.u12 = None  # assign after U12 file check
        print("Checking input files correctness...")
        self.__check_param_files()

        # create symlinks to 2bits: let user know what 2bits were used
        self.t_2bit_link = os.path.join(self.wd, "t2bit.link")
        self.q_2bit_link = os.path.join(self.wd, "q2bit.link")
        self.__make_symlink(self.t_2bit, self.t_2bit_link)
        self.__make_symlink(self.q_2bit, self.q_2bit_link)

        # dump input parameters, object state
        self.toga_params_file = os.path.join(self.wd, "toga_init_state.json")
        self.toga_args_file = os.path.join(self.wd, "project_args.json")
        with open(self.toga_params_file, "w") as f:
            # default=string is a workaround to serialize datetime object
            json.dump(self.__dict__, f, default=str)
        with open(self.toga_args_file, "w") as f:
            json.dump(vars(args), f, default=str)
        print("#### TOGA initiated successfully! ####")

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

    def __gen_project_name(self):
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        print(f"Using automatically generated project name: {project_name}")
        return project_name

    def __check_param_files(self):
        """Check that all parameter files exist."""
        files_to_check = [self.u12,
                          self.t_2bit,
                          self.q_2bit,
                          self.cesar_binary,
                          self.ref_bed,
                          self.chain_file,
                          self.isoforms_arg]
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
            print(f"Detected {len(twobit_sequences)} sequences in {two_bit_file}")
        except EOFError as err:  # this is a file but twobit reader couldn't read it
            twobit_seq_to_size = None  # to suppress linter
            twobit_sequences = None
            print(str(err))
            print(f"twobitreader cannot open {two_bit_file}")
            self.die("Abort")
        # another check: that bed or chain chromosomes intersect 2bit file sequences
        check_chroms = chroms_sizes.keys()
        intersection = twobit_sequences.intersection(check_chroms)
        if len(intersection) == 0:
            err = f"Error! 2bit file: {two_bit_file}; chain/bed file: {chrom_file}; " \
                  f"Different sets of chromosomes!"
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
            err = f"Error! 2bit file: {two_bit_file}; chain_file: {chrom_file} " \
                  f"Chromosome: {chrom}; Sizes don't match! " \
                  f"Size in twobit: {twobit_seq_len}; size in chain: {comp_file_seq_len}"
            self.die(err)
        return

    def __check_u12_file(self, t_in_bed):
        """Sanity check for U12 file."""
        if not self.u12_arg:
            # just not provided: nothing to check
            return
        # U12 file provided
        self.u12 = os.path.join(self.wd, "u12_data.txt")
        filt_lines = []
        f = open(self.u12_arg, "r")
        for num, line in enumerate(f, 1):
            line_data = line.rstrip().split("\t")
            if len(line_data) != U12_FILE_COLS:
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, 3 fields expected; "\
                          f"Got {len(line_data)}; please note that a tab-separated file expected"
                self.die(err_msg)
            trans_id = line_data[0]
            if trans_id not in t_in_bed:
                # transcript doesn't appear in the bed file: skip it
                continue
            exon_num = line_data[1]
            if not exon_num.isnumeric():
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, field 2 value is {exon_num}; "\
                          f"This field must contain a numeric value (exon number)."
                self.die(err_msg)
            acc_don = line_data[2]
            if acc_don not in U12_AD_FIELD:
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, field 3 value is {acc_don}; "\
                          f"This field could have either A or D value."
                self.die(err_msg)
            filt_lines.append(line)  # save this line
        f.close()
        # another check: what if there are no lines after filter?
        if len(filt_lines) == 0:
            err_msg = f"Error! No lines left in the {self.u12_arg} file after filter." \
                      f"Please check that transcript IDs in this file and bed {self.ref_bed} are consistent"
            self.die(err_msg)
        with open(self.u12, "w") as f:
            f.write("".join(filt_lines))
        print("U12 file is correct")

    def __check_isoforms_file(self, t_in_bed):
        """Sanity checks for isoforms file."""
        if not self.isoforms_arg:
            return  # not provided: nothing to check
        # isoforms file provided: need to check correctness and completeness
        # then check isoforms file itself
        f = open(self.isoforms_arg, "r")
        self.isoforms = os.path.join(self.wd, "isoforms.tsv")
        header = f.__next__()  # first line is header
        filt_isoforms_lines = [header, ]  # remove isoforms that don't appear in the bed file
        # also we catch isoforms that are in the bed but not in the isoforms file
        t_in_i = []
        for num, line in enumerate(f, 2):
            line_data = line.rstrip().split("\t")
            if len(line_data) != ISOFORMS_FILE_COLS:
                err_msg = f"Error! Isoforms file {self.isoforms} line {num}: " \
                          f"Expected {ISOFORMS_FILE_COLS} fields, got {len(line_data)}"
                self.die(err_msg)
            transcript = line_data[1]
            if transcript in t_in_bed:
                # this isoforms appears in the bed file: keep it
                filt_isoforms_lines.append(line)
                t_in_i.append(line_data[1])
            else:  # this isoform doesn't appear in the bed: we can skip it
                continue
        f.close()
        # this set contains isoforms found in the isoforms file
        t_in_i = set(t_in_i)
        # there are transcripts that appear in bed but not in the isoforms file
        # if this set is non-empty: raise an error
        u_in_b = t_in_bed.difference(t_in_i)

        if len(u_in_b) != 0:  # isoforms file is incomplete
            extra_t_list = "\n".join(list(u_in_b)[:100])  # show first 100 (or maybe show all?)
            err_msg = f"Error! There are {len(u_in_b)} transcripts in the bed file absent in the isoforms file! " \
                      f"There are the transcripts (first 100):\n{extra_t_list}"
            self.die(err_msg)
        # write isoforms file
        with open(self.isoforms, "w") as f:
            f.write("".join(filt_isoforms_lines))
        print("Isoforms file is OK")

    def die(self, msg, rc=1):
        """Show msg in stderr, exit with the rc given."""
        print(msg)
        print(f"Program finished with exit code {rc}\n")
        # for t_file in self.temp_files:  # remove temp files if required
        #     os.remove(t_file) if os.path.isfile(t_file) and not self.keep_temp else None
        #     shutil.rmtree(t_file) if os.path.isdir(t_file) and not self.keep_temp else None
        sys.exit(rc)

    def __modules_addr(self):
        """Define addresses of modules."""
        self.LOCATION = os.path.dirname(__file__)  # folder containing pipeline scripts
        self.CONFIGURE = os.path.join(self.LOCATION, "configure.sh")
        self.CHAIN_SCORE_FILTER = os.path.join(self.LOCATION, "modules", "chain_score_filter")
        self.CHAIN_COORDS_CONVERT_LIB = os.path.join(self.LOCATION,
                                                     "modules",
                                                     "chain_coords_converter_slib.so")
        self.EXTRACT_SUBCHAIN_LIB = os.path.join(self.LOCATION, "modules", "extract_subchain_slib.so")
        self.CHAIN_FILTER_BY_ID = os.path.join(self.LOCATION, "modules", "chain_filter_by_id")

        self.CHAIN_BDB_INDEX = os.path.join(self.LOCATION, "modules", "chain_bst_index.py")
        self.CHAIN_INDEX_SLIB = os.path.join(self.LOCATION, "modules", "chain_bst_lib.so")
        self.BED_BDB_INDEX = os.path.join(self.LOCATION, "modules", "bed_hdf5_index.py")
        self.SPLIT_CHAIN_JOBS = os.path.join(self.LOCATION, "split_chain_jobs.py")
        self.MERGE_CHAINS_OUTPUT = os.path.join(self.LOCATION, "modules", "merge_chains_output.py")
        self.CLASSIFY_CHAINS = os.path.join(self.LOCATION, "modules", "classify_chains.py")
        self.SPLIT_EXON_REALIGN_JOBS = os.path.join(self.LOCATION, "split_exon_realign_jobs.py")
        self.MERGE_CESAR_OUTPUT = os.path.join(self.LOCATION, "modules", "merge_cesar_output.py")
        self.TRANSCRIPT_QUALITY = os.path.join(self.LOCATION, "modules", "get_transcripts_quality.py")
        self.GENE_LOSS_SUMMARY = os.path.join(self.LOCATION, "modules", "gene_losses_summary.py")
        self.ORTHOLOGY_TYPE_MAP = os.path.join(self.LOCATION, "modules", "orthology_type_map.py")
        self.MODEL_TRAINER = os.path.join(self.LOCATION, "train_model.py")
        self.DEFAULT_CESAR = os.path.join(self.LOCATION, "CESAR2.0", "cesar")
        self.nextflow_rel_ = os.path.join(self.LOCATION, "execute_joblist.nf")
        self.NF_EXECUTE = os.path.abspath(self.nextflow_rel_)

    def __check_dependencies(self):
        """Check all dependencies."""
        print("check if binaries are compiled and libs are installed...")
        c_not_compiled = any(os.path.isfile(f) is False for f in [self.CHAIN_SCORE_FILTER,
                                                                  self.CHAIN_COORDS_CONVERT_LIB,
                                                                  self.CHAIN_FILTER_BY_ID,
                                                                  self.EXTRACT_SUBCHAIN_LIB,
                                                                  self.CHAIN_INDEX_SLIB])
        if c_not_compiled:
            print("Warning! C code is not compiled, trying to compile...")
        imports_not_found = False
        try:
            import twobitreader
            import networkx
            import pandas
            import xgboost
            import joblib
            import h5py
        except ImportError:
            print("Warning! Some of the required packages are not installed.")
            imports_not_found = True

        not_all_found = any([c_not_compiled, imports_not_found])
        self.__call_proc(self.CONFIGURE, "Could not call configure.sh!")\
            if not_all_found else print("All dependencies found")

    def __check_completeness(self):
        """Check if all modules are presented."""
        files_must_be = [self.CONFIGURE,
                         self.CHAIN_BDB_INDEX,
                         self.BED_BDB_INDEX,
                         self.SPLIT_CHAIN_JOBS,
                         self.MERGE_CHAINS_OUTPUT,
                         self.CLASSIFY_CHAINS,
                         self.SPLIT_EXON_REALIGN_JOBS,
                         self.MERGE_CESAR_OUTPUT,
                         self.GENE_LOSS_SUMMARY,
                         self.ORTHOLOGY_TYPE_MAP]
        for _file in files_must_be:
            if os.path.isfile(_file):
                continue
            self.die(f"Error! File {_file} not found!")
    
    def __check_nf_config(self):
        """Check that nextflow configure files are here."""
        if self.nextflow_bigmem_config and not os.path.isfile(self.nextflow_bigmem_config):
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
            self.die("Conflict: --para and --nf_dir should not be used at the same time")
        # check that required config files are here
        if not os.path.isdir(self.nextflow_config_dir):
            self.die(f"Error! Nextflow config dir {self.nextflow_config_dir} does not exist!")
        err_msg = "Please note these two files are expected in the nextflow config directory:\n" \
                  "1) call_cesar_config_template.nf" \
                  "2) extract_chain_features_config.nf"
        # check CESAR config template first
        nf_cesar_config_temp = os.path.join(self.nextflow_config_dir, "call_cesar_config_template.nf")
        if not os.path.isfile(nf_cesar_config_temp):
            self.die(f"Error! File {nf_cesar_config_temp} not found!\n{err_msg}")
        # we need the content of this file
        with open(nf_cesar_config_temp, "r") as f:
            self.cesar_config_template = f.read()
        # check chain extract features config; we need abspath to this file
        self.nf_chain_extr_config_file = os.path.abspath(os.path.join(self.nextflow_config_dir,
                                                                      "extract_chain_features_config.nf"))
        if not os.path.isfile(self.nf_chain_extr_config_file):
            self.die(f"Error! File {self.nf_chain_extr_config_file} not found!\n{err_msg}")
        self.local_executor = False

    def __call_proc(self, cmd, extra_msg=None):
        """Call a subprocess and catch errors."""
        print(f"{cmd} in progress...")
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            print(extra_msg) if extra_msg else None
            self.die(f"Error! Process {cmd} died! Abort.")
        print(f"{cmd} done with code 0")

    def __find_two_bit(self, db):
        """Find a 2bit file."""
        if os.path.isfile(db):
            return db
        # For now here is a hillerlab-oriented solution
        # you can write your own template for 2bit files location
        with_alias = f"/projects/hillerlab/genome/gbdb-HL/{db}/{db}.2bit"
        if os.path.isfile(with_alias):
            return with_alias
        self.die(f"Two bit file {db} not found! Abort")

    def run(self):
        """Run toga. Method to be called."""
        # 0) preparation:
        # define the project name and mkdir for it
        # move chain file filtered
        # define initial values
        # make indexed files for the chain
        print("#### STEP 0: making chain and bed file indexes\n")
        self.__make_indexed_chain()
        self.__make_indexed_bed()
        self.__time_mark("Made indexes")

        # 1) make joblist for chain features extraction
        print("#### STEP 1: Generate extract chain features jobs\n")
        self.__split_chain_jobs()
        self.__time_mark("Split chain jobs")

        # 2) extract chain features: parallel process
        print("#### STEP 2: Extract chain features: parallel step\n")
        self.__extract_chain_features()
        self.__time_mark("Chain jobs done")

        # 3) create chain features dataset
        print("#### STEP 3: Merge step 2 output\n")
        self.__merge_chains_output()
        self.__time_mark("Chains output merged")

        # 4) classify chains as orthologous, paralogous, etc using xgboost
        print("#### STEP 4: Classify chains using gradient boosting model\n")
        self.__classify_chains()
        self.__time_mark("Chains classified")

        # 5) create cluster jobs for CESAR2.0
        print("#### STEP 5: Generate CESAR jobs")
        self.__split_cesar_jobs()
        self.__time_mark("Split cesar jobs done")

        # 6) Create bed track for processed pseudogenes
        print("# STEP 6: RESERVED (PLACEHOLDER)\n")
        # self.__get_proc_pseudogenes_track()

        # 7) call CESAR jobs: parallel step
        print("#### STEP 7: Execute CESAR jobs: parallel step (the most time consuming)\n")
        self.__run_cesar_jobs()
        self.__time_mark("Cesar jobs done")

        # 8) parse CESAR output, create bed / fasta files
        print("#### STEP 8: Merge STEP 7 output\n")
        self.__merge_cesar_output()
        self.__time_mark("Merged cesar output")

        # 9) classify projections/genes as lost/intact
        # also measure projections confidence levels
        print("#### STEP 9: Gene loss pipeline classification\n")
        self.__transcript_quality()
        self.__gene_loss_summary()
        self.__time_mark("Got gene loss summary")

        # 10) classify genes as one2one, one2many, etc orthologs
        print("#### STEP 10: Create orthology relationships table\n")
        self.__orthology_type_map()

        # 11) merge logs containing information about skipped genes,transcripts, etc.
        print("#### STEP 11: Cleanup: merge parallel steps output files")
        self.__merge_split_files()
        # Everything is done

        self.__time_mark("Everything is done")
        if not self.cesar_ok_merged:
            print("PLEASE NOTE: SOME CESAR JOBS CRASHED")
            print("RESULTS ARE LIKELY INCOMPLETE")
        self.die(f"Done! Estimated time: {dt.now() - self.t0}", rc=0)

    def __make_indexed_chain(self):
        """Make chain index file."""
        # make *.bb file
        print("make_indexed in progress...")
        chain_bst_index(self.chain_file,
                        self.chain_index_file,
                        txt_index=self.chain_index_txt_file)
        self.temp_files.append(self.chain_index_file)
        self.temp_files.append(self.chain_file)
        self.temp_files.append(self.chain_index_txt_file)
        print("Indexed")

    def __time_mark(self, msg):
        """Left time mark."""
        if self.time_log is None:
            return
        t = dt.now() - self.t0
        with open(self.time_log, "a") as f:
            f.write(f"{msg} at {t}\n")

    def __make_indexed_bed(self):
        """Create gene_ID: bed line bdb indexed file."""
        print("index_bed in progress...")
        bed_hdf5_index(self.ref_bed, self.index_bed_file)
        self.temp_files.append(self.index_bed_file)
        print("Bed file indexed")

    def __split_chain_jobs(self):
        """Wrap split_jobs.py script."""
        # define arguments
        # save split jobs
        ch_cl_jobs = os.path.join(self.wd, "chain_classification_jobs")
        # for raw results of this stage
        self.chain_class_results = os.path.join(self.wd, "chain_classification_results")
        self.chain_cl_jobs_combined = os.path.join(self.wd, "chain_class_jobs_combined")
        rejected_filename = "SPLIT_CHAIN_REJ.txt"
        rejected_path = os.path.join(self.rejected_dir,
                                     rejected_filename)
        self.temp_files.append(ch_cl_jobs)
        self.temp_files.append(self.chain_class_results)
        self.temp_files.append(self.chain_cl_jobs_combined)

        split_jobs_cmd = f"{self.SPLIT_CHAIN_JOBS} {self.chain_file} " \
                         f"{self.ref_bed} {self.index_bed_file} " \
                         f"--jobs_num {self.chain_jobs} " \
                         f"--jobs {ch_cl_jobs} " \
                         f"--jobs_file {self.chain_cl_jobs_combined} " \
                         f"--results_dir {self.chain_class_results} " \
                         f"--rejected {rejected_path}"

        self.__call_proc(split_jobs_cmd, "Could not split chain jobs!")
    
    def __extract_chain_features(self):
        """Execute extract chain features jobs."""
        # get timestamp to name the project and create a dir for that
        #  time() returns something like: 1595861493.8344169
        timestamp = str(time.time()).split(".")[0]
        project_name = f"{self.project_name}_chain_feats_at_{timestamp}"
        project_path = os.path.join(self.nextflow_dir, project_name)
        print(f"Extract chain features, project name: {project_name}")
        print(f"Project path: {project_path}")

        if self.para:  # run jobs with para, skip nextflow
            cmd = f"para make {project_name} {self.chain_cl_jobs_combined} -q=\"short\""
            print(f"Calling {cmd}")
            rc = subprocess.call(cmd, shell=True)
        else:  # calling jobs with nextflow
            cmd = f"nextflow {self.NF_EXECUTE} " \
                  f"--joblist {self.chain_cl_jobs_combined}"
            if not self.local_executor:
                # not local executor -> provided config files
                # need abspath for nextflow execution
                cmd += f" -c {self.nf_chain_extr_config_file}"
            print(f"Calling {cmd}")
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
        print("Merging chain output...")        
        merge_chains_output(self.ref_bed, self.isoforms,
                            self.chain_class_results, self.chain_results_df)
        # .append(self.chain_results_df)  -> UCSC plugin needs that

    def __classify_chains(self):
        """Run decision tree."""
        # define input and output."""
        print("Decision tree in progress...")
        self.orthologs = os.path.join(self.wd, "trans_to_chain_classes.tsv")
        self.pred_scores = os.path.join(self.wd, "orthology_scores.tsv")
        self.se_model = os.path.join(self.LOCATION, "models", "se_model.dat")
        self.me_model = os.path.join(self.LOCATION, "models", "me_model.dat")
        cl_rej_log = os.path.join(self.rejected_dir, "classify_chains_rejected.txt")
        if not os.path.isfile(self.se_model) or not os.path.isfile(self.me_model):
            self.__call_proc(self.MODEL_TRAINER, "Models not found, training...")
        classify_chains(self.chain_results_df, self.orthologs, self.se_model,
                        self.me_model, rejected=cl_rej_log, raw_out=self.pred_scores)
        if self.stop_at_chain_class:
            self.die("User requested to halt after chain features extraction", rc=0)

    def __get_proc_pseudogenes_track(self):
        """Create annotation of processed genes in query."""
        print("Creating processed pseudogenes track.")
        proc_pgenes_track = os.path.join(self.wd, "proc_pseudogenes.bed")
        create_ppgene_track(self.orthologs, self.chain_file, self.index_bed_file, proc_pgenes_track)

    def __split_cesar_jobs(self):
        """Call split_exon_realign_jobs.py."""
        if not self.t_2bit or not self.q_2bit:
            self.die("There is no 2 bit files provided, cannot go ahead and call CESAR.", rc=0)

        if self.fragmented_genome:
            print("Detect fragmented transcripts")
            # need to stitch fragments together
            gene_fragments = stitch_scaffolds(self.chain_file, self.pred_scores, self.ref_bed, True)
            # for now split exon realign jobs is a subprocess, not a callable function
            # TODO: fix this
            fragm_dict_file = os.path.join(self.wd, "gene_fragments.txt")
            f = open(fragm_dict_file, "w")
            for k, v in gene_fragments.items():
                v_str = ",".join(map(str, v))
                f.write(f"{k}\t{v_str}\n")
            f.close()
            print(f"Detected {len(gene_fragments.keys())} fragmented transcripts")
            print(f"Fragments data saved to {fragm_dict_file}")
        else:
            # no fragment file: ok
            print("Skip fragmented genes detection")
            fragm_dict_file = None
        # if we call CESAR
        cesar_jobs_dir = os.path.join(self.wd, "cesar_jobs")
        self.cesar_combined = os.path.join(self.wd, "cesar_combined")
        self.cesar_results = os.path.join(self.wd, "cesar_results")
        self.cesar_bigmem_jobs = os.path.join(self.wd, "cesar_bigmem")
        self.temp_files.append(cesar_jobs_dir)
        self.temp_files.append(self.cesar_combined)

        # different field names depending on --ml flag
        fields = "ORTH,TRANS"

        self.temp_files.append(self.cesar_results)
        self.temp_files.append(self.gene_loss_data)
        skipped_path = os.path.join(self.rejected_dir, "SPLIT_CESAR.txt")
        self.paralogs_log = os.path.join(self.wd, "paralogs.txt")

        split_cesar_cmd = f"{self.SPLIT_EXON_REALIGN_JOBS} "\
                          f"{self.orthologs} {self.ref_bed} " \
                          f"{self.index_bed_file} {self.chain_index_file} " \
                          f"{self.t_2bit} {self.q_2bit} " \
                          f"--jobs_dir {cesar_jobs_dir} " \
                          f"--jobs_num {self.cesar_jobs_num} " \
                          f"--combined {self.cesar_combined} " \
                          f"--bigmem {self.cesar_bigmem_jobs} " \
                          f"--results {self.cesar_results} " \
                          f"--buckets {self.cesar_buckets} " \
                          f"--mem_limit {self.cesar_mem_limit} " \
                          f"--chains_limit {self.cesar_chain_limit} " \
                          f"--fields {fields} " \
                          f"--skipped_genes {skipped_path} " \
                          f"--rejected_log {self.rejected_dir} " \
                          f"--cesar_binary {self.cesar_binary} " \
                          f"--paralogs_log {self.paralogs_log} " \
                          f"--uhq_flank {self.uhq_flank}"

        split_cesar_cmd = split_cesar_cmd + f" --cesar_binary {self.cesar_binary}" \
            if self.cesar_binary else split_cesar_cmd
        split_cesar_cmd = split_cesar_cmd + " --mask_stops" if self.mask_stops \
            else split_cesar_cmd
        split_cesar_cmd = split_cesar_cmd + f" --u12 {self.u12}" if self.u12 \
            else split_cesar_cmd
        split_cesar_cmd = split_cesar_cmd + " --o2o_only" if self.o2o_only \
            else split_cesar_cmd
        split_cesar_cmd = split_cesar_cmd + " --no_fpi" if self.no_fpi \
            else split_cesar_cmd
        if self.gene_loss_data:
            split_cesar_cmd += f" --check_loss {self.gene_loss_data}"
        if fragm_dict_file:
            split_cesar_cmd += f" --fragments_data {fragm_dict_file}"
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

        # get a list of buckets
        if self.cesar_buckets == "0":
            buckets = [0, ]  # a single bucket
        else:  # several buckets, each int -> memory limit in gb
            buckets = [int(x) for x in self.cesar_buckets.split(",") if x != ""]
        print(f"Pushing {len(buckets)} joblists")

        for b in buckets:
            # create config file
            # 0 means that that buckets were not split
            mem_lim = b if b != 0 else self.cesar_mem_limit

            if not self.local_executor:
                # running on cluster, need to create config file
                # for this bucket's memory requirement
                config_string = self.cesar_config_template.replace("${_MEMORY_}", f"{mem_lim}")
                config_file_path = os.path.join(self.wd, f"cesar_config_{b}_queue.nf")
                config_file_abspath = os.path.abspath(config_file_path)
                with open(config_file_path, "w") as f:
                    f.write(config_string)
                self.temp_files.append(config_file_path)
            else:  # no config dir given: use local executor
                # OK if there is a single bucket
                config_file_abspath = None

            if b != 0:  # extract jobs related to this bucket (if it's not 0)
                bucket_tasks = self.__get_cesar_jobs_for_bucket(self.cesar_combined, str(b))
                if len(bucket_tasks) == 0:
                    print(f"There are no jobs in the {b} bucket")
                    continue
                joblist_name = f"cesar_joblist_queue_{b}.txt"
                joblist_path = os.path.join(self.wd, joblist_name)
                with open(joblist_path, "w") as f:
                    f.write(bucket_tasks)
                joblist_abspath = os.path.abspath(joblist_path)
                self.temp_files.append(joblist_path)
            else:  # nothing to extract, there is a single joblist
                joblist_abspath = os.path.abspath(self.cesar_combined)

            # create project directory for logs
            nf_project_name = f"{self.project_name}_cesar_at_{timestamp}_q_{b}"
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            project_paths.append(nf_project_path)

            if not self.para:  # create nextflow cmd
                os.mkdir(nf_project_path) if not os.path.isdir(nf_project_path) else None

                cmd = f"nextflow {self.NF_EXECUTE} " \
                      f"--joblist {joblist_abspath}"
                if config_file_abspath:
                    cmd += f" -c {config_file_abspath}"
                p = subprocess.Popen(cmd, shell=True, cwd=nf_project_path)
            else:  # create cmd for para
                memory_mb = b * 1000
                cmd = f"para make {nf_project_name} {joblist_abspath} -q=\"day\""
                if memory_mb > 0:
                    cmd += f" --memoryMb={memory_mb}"
                p = subprocess.Popen(cmd, shell=True)
            # create subprocess object
            sys.stderr.write(f"Pushed cluster jobs with {cmd}\n")
            processes.append(p)
            time.sleep(CESAR_PUSH_INTERVAL)

        if self.nextflow_bigmem_config and not self.para:
            # if provided: push bigmem jobs also
            nf_project_name = f"{self.project_name}_cesar_at_{timestamp}_q_bigmem"
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            nf_cmd = f"nextflow {self.NF_EXECUTE} " \
                     f"--joblist {self.cesar_bigmem_jobs} -c {self.nextflow_bigmem_config}"
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
                os.mkdir(nf_project_path) if not os.path.isdir(nf_project_path) else None
                project_paths.append(nf_project_path)
                p = subprocess.Popen(nf_cmd, shell=True, cwd=nf_project_path)
                sys.stderr.write(f"Pushed {big_lines_num} bigmem jobs with {nf_cmd}\n")
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
                cmd = f"para make {bm_project_name} {self.cesar_bigmem_jobs} "\
                      f"-q=\"bigmem\" --memoryMb={memory_mb}"
                p = subprocess.Popen(cmd, shell=True)
                processes.append(p)

        # monitor jobs, iteration counter
        iter_num = 0
        while True:  # Run until all jobs are done (or crashed)
            all_done = True  # default val, re-define if something is not done
            for p in processes:
                # check if each process is still running
                running = p.poll() is None
                if running:
                    all_done = False
            if all_done:
                print("CESAR jobs done")
                break
            else:
                print(f"Iter {iter_num} waiting {ITER_DURATION * iter_num} seconds Not done")
                time.sleep(ITER_DURATION)
                iter_num += 1
        if not self.keep_nf_logs:
            # remove nextflow intermediate files
            print("Removing nextflow temp files")
            for path in project_paths:
                shutil.rmtree(path) if os.path.isdir(path) else None

    def __merge_cesar_output(self):
        """Merge CESAR output, save final fasta and bed."""
        print("Merging CESAR output to make fasta and bed files.")
        merge_c_stage_skipped = os.path.join(self.rejected_dir,
                                             "CESAR_MERGE.txt")
        self.temp_files.append(self.intermediate_bed)

        all_ok = merge_cesar_output(self.cesar_results,
                                    self.intermediate_bed,
                                    self.nucl_fasta,
                                    self.meta_data,
                                    merge_c_stage_skipped,
                                    self.prot_fasta,
                                    self.codon_fasta,
                                    self.trash_exons,
                                    fragm_data=self.bed_fragm_exons_data)
        if all_ok:
            # there are no empty output files
            print("CESAR results merged")
            self.cesar_ok_merged = True
        else:
            # there are some empty output files
            # MAYBE everything is fine
            # but need to notify user anyway
            print("WARNING!\nSOME CESAR JOBS LIKELY CRASHED\n!")
            print("RESULTS ARE LIKELY INCOMPLETE")
            print("PLEASE SEE LOGS FOR DETAILS")
            self.cesar_ok_merged = False

    def __transcript_quality(self):
        """Call module to get transcript quality."""
        self.trans_quality_file = os.path.join(self.wd, "transcript_quality.tsv")
        classify_transcripts(self.meta_data,
                             self.pred_scores,
                             self.hq_orth_threshold,
                             self.trans_quality_file)

    def __gene_loss_summary(self):
        """Call gene loss summary."""
        print("Calling gene loss summary")
        gene_losses_summary(self.gene_loss_data,
                            self.ref_bed,
                            self.intermediate_bed,
                            self.query_annotation,
                            self.loss_summ,
                            iforms=self.isoforms,
                            paral=self.paralogs_log)

    def __orthology_type_map(self):
        """Call orthology_type_map.py"""
        # need to combine projections in genes
        query_isoforms_file = os.path.join(self.wd, "query_isoforms.tsv")
        query_gene_spans = os.path.join(self.wd, "query_gene_spans.bed")
        get_query_isoforms_data(self.query_annotation, query_isoforms_file, save_genes_track=query_gene_spans)
        print("Calling orthology_type_map...")
        skipped_ref_trans = os.path.join(self.wd, "ref_orphan_transcripts.txt")
        orthology_type_map(self.ref_bed,
                           self.query_annotation,
                           self.orthology_type,
                           ref_iso=self.isoforms,
                           que_iso=query_isoforms_file,
                           paralogs_arg=self.paralogs_log,
                           loss_data=self.loss_summ,
                           save_skipped=skipped_ref_trans)

    @staticmethod
    def __merge_dir(dir_name, output):
        """Merge all files in a directory into one."""
        files_list = os.listdir(dir_name)
        if len(files_list) == 0:
            sys.exit(f"Error! {dir_name} is empty")
        buffer = open(output, "w")
        for filename in files_list:
            path = os.path.join(dir_name, filename)
            with open(path, "r") as f:
                content = [x for x in f.readlines() if x != "\n"]
            lines = "".join(content)
            buffer.write(lines)
        buffer.close()
        shutil.rmtree(dir_name)

    def __get_version(self):
        """Get git hash if possible."""
        cmd = "git rev-parse HEAD"
        try:
            version = subprocess.check_output(cmd, shell=True, cwd=self.toga_exe_path).decode("utf-8")
        except subprocess.CalledProcessError:
            version = "unknown"
        return version

    def __merge_split_files(self):
        """Merge intermediate/temp files."""
        # merge rejection logs
        rejected_log = os.path.join(self.wd, "genes_rejection_reason.tsv")
        self.__merge_dir(self.rejected_dir, rejected_log)
        # save inact mutations data
        inact_mut_file = os.path.join(self.wd, "inact_mut_data.txt")
        self.__merge_dir(self.gene_loss_data, inact_mut_file)
        # save CESAR outputs
        cesar_results_merged = os.path.join(self.wd, "cesar_results.txt")
        self.__merge_dir(self.cesar_results, cesar_results_merged)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_input", type=str,
                     help="Chain file. Extensions like "
                          "FILE.chain or FILE.chain.gz are applicable.")
    app.add_argument("bed_input", type=str,
                     help="Bed file with annotations for the target genome.")
    app.add_argument("tDB", default=None,
                     help="Reference genome sequence in 2bit format.")
    app.add_argument("qDB", default=None,
                     help="Query genome sequence in 2bit format.")
    # global ops
    app.add_argument("--project_dir", "--pd", default=None,
                     help="Project directory. TOGA will save all intermediate and output files "
                          "exactly in this directory. If not specified, use CURRENT_DIR/PROJECT_NAME "
                          "as default (see below).")
    app.add_argument("--project_name", "--pn", default=None,
                     help="If you don't like to provide a full path to the project directory with "
                          "--project_dir you can use this parameter. In this case TOGA will "
                          "create project directory in the current directory as "
                          "\"CURRENT_DIR/PROJECT_NAME\". If not provided, TOGA will try to extract "
                          "the project name from chain filename, which is not recommended.")
    app.add_argument("--min_score", "--msc", type=int, default=15000,
                     help="Chain score threshold. Exclude chains that have a lower score "
                     "from the analysis. Default value is 15000.")
    app.add_argument("--isoforms", "-i", type=str, default="",
                     help="Path to isoforms data file")
    app.add_argument("--keep_temp", "--kt", action="store_true",
                     dest="keep_temp", help="Do not remove temp files.")
    app.add_argument("--limit_to_ref_chrom", default=None, help="Find orthologs "
                     "for a single reference chromosome only")
    # app.add_argument("--limit_to_query_chrom", default=None, help="Annotate "
    #                  "a particular query scaffold/chromosome only")
    # nextflow related
    app.add_argument("--nextflow_dir", "--nd", default=None,
                     help="Nextflow working directory: from this directory nextflow is "
                          "executed, also there all nextflow log files are kept")
    app.add_argument("--nextflow_config_dir", "--nc", default=None,
                     help="Directory containing nextflow configuration files "
                          "for cluster, pls see nextflow_config_files/readme.txt "
                          "for details.")
    app.add_argument("--do_not_del_nf_logs", "--nfnd", action="store_true",
                     dest="do_not_del_nf_logs")
    app.add_argument("--cesar_bigmem_config", "--nb", default=None, help="File containing "
                     "nextflow config for BIGMEM CESAR jobs. If not provided, these "
                     "jobs will not run (but list of them saved)/ NOT IMPLEMENTED YET")
    app.add_argument("--para", "-p", action="store_true", dest="para",
                     help="Hillerlab feature, use para instead of nextflow to "
                          "manage cluster jobs.")
    app.add_argument("--para_bigmem", "--pb", action="store_true", dest="para_bigmem",
                     help="Hillerlab feature, push bigmem jobs with para")
    # chain features related
    app.add_argument("--chain_jobs_num", "--chn", type=int, default=100,
                     help="Number of cluster jobs for extracting chain features. "
                          "Recommended from 150 to 200 jobs.")
    app.add_argument("--no_chain_filter", "--ncf", action="store_true",
                     dest="no_chain_filter",
                     help="A flag. Do not filter the chain file (make sure you specified "
                          "a .chain but not .gz file in this case)")
    # CESAR part related
    app.add_argument("--cesar_jobs_num", "--cjn", type=int, default=500,
                     help="Number of CESAR cluster jobs.")
    app.add_argument("--cesar_binary", default=None,
                     help="CESAR binary address, cesar as default.")
    app.add_argument("--mask_stops", "--ms", action="store_true", dest="mask_stops",
                     help="Mask stop codons in target sequences. "
                          "CESAR cannot process them. Using this "
                          "parameter please make sure you know what you are doing.")
    app.add_argument("--cesar_buckets", "--cb", default="0",
                     help="Comma-separated list of integers. Split CESAR jobs in buckets "
                          "depending on their memory requirements. "
                          "See README.md for explanation.")
    app.add_argument("--cesar_chain_limit", type=int, default=100,
                     help="Skip genes that have more that CESAR_CHAIN_LIMIT orthologous "
                          "chains. Recommended values are a 50-100.")
    app.add_argument("--cesar_mem_limit", type=int, default=16,
                     help="Ignore genes requiring > N gig to run CESAR")
    app.add_argument("--homology_types", "--ht", default="ORTH,TRANS",
                     help="TOGA classifies chain in the following categories: "
                          "ORTH, PARA ans TRANS. This parameter controls, "
                          "for which groups TOGA will seek for genes. For "
                          "example, ORTH,TRANS means that TOGA will seek for genes "
                          "throw chains which classified as ORTH and TRANS")
    app.add_argument("--time_marks", "-t", default=None,
                     help="File to save timings of different steps.")
    app.add_argument("--u12", default=None, help="U12 introns data")
    app.add_argument("--stop_at_chain_class", "--sac", action="store_true",
                     dest="stop_at_chain_class",
                     help="Stop after merging chain features.")
    app.add_argument("--uhq_flank", default=50, type=int, help="Flank size for UHQ exons")
    app.add_argument("--o2o_only", "--o2o", action="store_true", dest="o2o_only",
                     help="Process only the genes that have a single orthologous chain")
    app.add_argument("--no_fpi", action="store_true", dest="no_fpi",
                     help="Consider some frame-preserving mutations as inactivating. "
                          "See documentation for details.")
    app.add_argument("--fragmented_genome", "-f", dest="fragmented_genome", action="store_true",
                     help="Annotation of a fragmented genome: need to assemble query genes "
                          "from pieces")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def main():
    """Entry point."""
    args = parse_args()
    toga_manager = Toga(args)
    toga_manager.run()


if __name__ == "__main__":
    main()
