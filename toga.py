#!/usr/bin/env python3
"""Master script for the TOGA pipeline.

Perform all operations from the beginning to the end.
If you need to call TOGA: most likely this is what you need.
"""
import argparse
import sys
import os
import shutil
import subprocess
import time
from datetime import datetime as dt
import json
from twobitreader import TwoBitFile
from modules.filter_bed import prepare_bed_file
from modules.chain_bdb_index import chain_bdb_index
from modules.bed_bdb_index import bed_bdb_index
from modules.merge_chains_output import merge_chains_output
from modules.merge_cesar_output import merge_cesar_output
from modules.gene_losses_summary import gene_losses_summary
from modules.orthology_type_map import orthology_type_map
from modules.classify_chains import classify_chains
from modules.get_transcripts_quality import classify_transcripts
from modules.make_query_isoforms import get_query_isoforms_data


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

U12_FILE_COLS = 3
U12_AD_FIELD = {"A", "D"}
ISOFORMS_FILE_COLS = 2
NF_DIR_NAME = "nextflow_logs"
CESAR_PUSH_INTERVAL = 30  # CESAR jobs push interval
ITER_DURATION = 120  # CESAR jobs check interval


class Toga:
    """TOGA manager class."""
    def __init__(self, args):
        """Initiate toga class."""
        self.t0 = dt.now()
        # check if all files TOGA needs are here
        self.temp_files = []  # remove at the end, list of temp files
        self.__modules_addr()
        self.__check_dependencies()
        self.__check_completeness()
        self.nextflow_dir = self.__get_nf_dir(args.nextflow_dir)
        self.nextflow_config_dir = args.nextflow_config_dir
        self.__check_nf_config()

        eprint("mkdir_and_move_chain in progress...")
        chain_basename = os.path.basename(args.chain_initial)

        # create project dir
        self.project_name = chain_basename.split(".")[1] if not args.project_name \
            else args.project_name
        self.wd = args.project_folder if args.project_folder else  \
            os.path.join(os.getcwd(), self.project_name)
        os.mkdir(self.wd) if not os.path.isdir(self.wd) else None

        # dir to collect log files with rejected reference genes:
        self.rejected_dir = os.path.join(self.wd, "rejected")
        os.mkdir(self.rejected_dir) if not os.path.isdir(self.rejected_dir) else None

        # filter chain in this folder
        chain_filename = chain_basename.replace(".gz", "")  # cut gz
        self.chain = os.path.join(self.wd, chain_filename)
        # there is an assumption that people will usually use .chain file
        # extension for chain files
        index_file = os.path.basename(self.chain).replace(".chain", ".bdb")
        self.chain_index_file = os.path.join(self.wd, index_file)

        # make the command, prepare the chain file
        if chain_basename.endswith(".gz"):  # version for gz
            chain_filter_cmd = f"gzip -dc {args.chain_initial} | "\
                               f"{self.CHAIN_SCORE_FILTER} stdin "\
                               f"{args.min_score} > {self.chain}"
        elif args.no_chain_filter:  # it is .chain and score filter is not required
            chain_filter_cmd = f"rsync -a {args.chain_initial} {self.chain}"
        else:  # it is .chain | score filter required
            chain_filter_cmd = f"{self.CHAIN_SCORE_FILTER} {args.chain_initial} "\
                               f"{args.min_score} > {self.chain}"

        # filter chains with score < threshold
        self.__call_proc(chain_filter_cmd, "Please check if you use a proper chain file.")
        
        # bed define bed files addresses
        self.ref_bed = os.path.join(self.wd, f"toga_filt_ref.{os.path.basename(args.bed_initial)}")
        index_bed_filename = os.path.basename(self.ref_bed).replace(".bed", ".bdb")
        self.index_bed_file = os.path.join(self.wd, index_bed_filename)

        # filter bed file
        bed_filt_rejected_file = "BED_FILTER_REJECTED.txt"
        bed_filt_rejected = os.path.join(self.rejected_dir,
                                         bed_filt_rejected_file)
        # keeping UTRs!
        prepare_bed_file(args.bed_initial,
                         self.ref_bed,
                         save_rejected=bed_filt_rejected,
                         only_chrom=args.limit_to_ref_chrom)

        # mics things
        self.isoforms = args.isoforms if args.isoforms else None
        self.chain_jobs = args.chain_jobs_num
        self.cesar_binary = self.DEFAULT_CESAR if not args.cesar_binary \
            else args.cesar_binary
        self.time_log = args.time_marks
        self.stop_at_chain_class = args.stop_at_chain_class

        self.keep_temp = True if args.keep_temp else False
        # define to call CESAR or not to call
        self.t_2bit = self.__find_two_bit(args.tDB)
        self.q_2bit = self.__find_two_bit(args.qDB)

        self.hq_orth_treshold = 0.95
        self.cesar_jobs_num = args.cesar_jobs_num
        self.cesar_buckets = args.cesar_buckets
        self.cesar_mem_limit = args.cesar_mem_limit
        self.cesar_chain_limit = args.cesar_chain_limit
        self.uhq_flank = args.uhq_flank
        self.cesar_fields = args.homology_types
        self.mask_stops = args.mask_stops
        self.no_fpi = args.no_fpi
        self.o2o_only = args.o2o_only

        self.chain_results_df = os.path.join(self.wd, "chain_results_df.tsv")
        self.nucl_fasta = os.path.join(self.wd, "nucleotide.fasta")
        self.prot_fasta = os.path.join(self.wd, "prot.fasta")
        self.final_bed = os.path.join(self.wd, "query_annotation.bed")
        self.low_conf_bed = os.path.join(self.wd, "low_confidence.bed")
        self.meta_data = os.path.join(self.wd, "exons_meta_data.tsv")
        self.intermediate_bed = os.path.join(self.wd, "intermediate.bed")
        self.orthology_type = os.path.join(self.wd, "orthology_classification.tsv")
        self.classification_log = os.path.join(self.wd, "o_class.log")
        self.trash_exons = os.path.join(self.wd, "trash_exons.bed")
        self.gene_loss_data = os.path.join(self.wd, "inact_mut_data")
        self.query_annotation = os.path.join(self.wd, "query_annotation.bed")
        self.loss_summ = os.path.join(self.wd, "loss_summ_data.tsv")
        self.u12 = args.u12
        self.__check_param_files()

        # dump input parameters, object state
        self.toga_params_file = os.path.join(self.wd, "toga_init_state.json")
        self.toga_args_file = os.path.join(self.wd, "project_args.json")
        with open(self.toga_params_file, "w") as f:
            # default=string is a workaround to serialize datetime object
            json.dump(self.__dict__, f, default=str)
        with open(self.toga_args_file, "w") as f:
            json.dump(vars(args), f, default=str)
        print("TOGA initiated successfully!")

    def __check_param_files(self):
        """Check that all parameter files exist."""
        files_to_check = [self.u12,
                          self.t_2bit,
                          self.q_2bit,
                          self.cesar_binary,
                          self.ref_bed,
                          self.chain,
                          self.isoforms]
        for item in files_to_check:
            if not item:
                # this file just not given
                continue
            elif not os.path.isfile(item):
                self.die(f"Error! File {item} not found!")

        # sanity checks
        self.__check_isoforms_file()
        self.__check_u12_file()
        self.__check_2bit_file(self.t_2bit)
        self.__check_2bit_file(self.q_2bit)

    def __get_nf_dir(self, nf_dir_arg):
        """Define nextflow directory."""
        if nf_dir_arg is None:
            default_dir = os.path.join(self.LOCATION, NF_DIR_NAME)
            os.mkdir(default_dir) if not os.path.isdir(default_dir) else None
            return default_dir
        else:
            os.mkdir(nf_dir_arg) if not os.path.isdir(nf_dir_arg) else None
            return nf_dir_arg

    def __check_2bit_file(self, two_bit_file):
        """Check that 2bit file is readable."""
        # TODO: check that chroms from ref bed are in the 2bit
        try:
            two_bit_reader = TwoBitFile(two_bit_file)
            seq_sizes = two_bit_reader.sequence_sizes()
            print(f"Detected {len(seq_sizes)} sequences in {two_bit_file}")
        except Exception as err:
            print(str(err))
            print(f"twobitreader cannot open {two_bit_file}")
            self.die("Abort")
        return
        
    def __check_u12_file(self):
        """Sanity check for U12 file."""
        if not self.u12:
            # just not provided: nothing to check
            return
        f = open(self.u12, "r")
        for num, line in enumerate(f, 1):
            line_data = line.rstrip().split("\t")
            if len(line_data) != U12_FILE_COLS:
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, 3 fields expected; "\
                          f"Got {len(line_data)}; please note that a tab-separated file expected"
                self.die(err_msg)
            _ = line_data[0]
            exon_num = line_data[1]
            if not exon_num.isnumeric():
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, field 2 value is {exon_num}; "\
                          f"This field must contain a numeric value (exon number)."
            acc_don = line_data[2]
            if acc_don not in U12_AD_FIELD:
                err_msg = f"Error! U12 file {self.u12} line {num} is corrupted, field 3 value is {acc_don}; "\
                          f"This field could have either A or D value."
                self.die(err_msg)
        f.close()
        eprint("U12 file is correct")
    
    def __check_isoforms_file(self):
        """Sanity checks for isoforms file."""
        if not self.isoforms:
            return  # not provided: nothing to check
        f = open(self.isoforms, "r")
        # this is the header, maybe we can check something there
        _ = f.__next__().rstrip().split("\t")
        t_in_i = []
        t_in_bed = []
        for num, line in enumerate(f, 2):
            line_data = line.rstrip().split("\t")
            if len(line_data) != ISOFORMS_FILE_COLS:
                err_msg = f"Error! Isoforms file {self.isoforms} line {num}: " \
                          f"Expected {ISOFORMS_FILE_COLS} fields, got {len(line_data)}"
                self.die(err_msg)
            t_in_i.append(line_data[1])
        f.close()
        # check with bed 12 file
        f = open(self.ref_bed, "r")
        for line in f:
            line_data = line.rstrip().split("\t")
            # sanity checked already
            t_in_bed.append(line_data[3])
        f.close()

        t_in_bed = set(t_in_bed)
        t_in_i = set(t_in_i)
        u_in_b = t_in_bed.difference(t_in_i)
        if len(u_in_b) != 0:
            extra_t_list = "\n".join(list(u_in_b)[:100])
            err_msg = f"Error! There are {len(u_in_b)} transctipts in the bed file absent in the isoforms file! " \
                      f"There are the transctipts (first 100):\n{extra_t_list}"
            self.die(err_msg)
            # eprint(err_msg)
        else:
            eprint("Isoforms file is OK")

    def die(self, msg, rc=1):
        """Show msg in stderr, exit with the rc given."""
        sys.stderr.write(msg + "\n")
        sys.stderr.write("Program finished with exit code {}\n".format(rc))
        for t_file in self.temp_files:  # remove temp files if required
            os.remove(t_file) if os.path.isfile(t_file) and not self.keep_temp else None
            shutil.rmtree(t_file) if os.path.isdir(t_file) and not self.keep_temp else None
        sys.exit(rc)

    def __modules_addr(self):
        """Define addresses of modules."""
        self.LOCATION = os.path.dirname(__file__)  # folder containing pipeline scripts
        self.CONFIGURE = os.path.join(self.LOCATION, "configure.sh")
        self.CHAIN_SCORE_FILTER = os.path.join(self.LOCATION, "modules", "chain_score_filter")
        self.CHAIN_COORDS_CONVERT_LIB = os.path.join(self.LOCATION, "modules",
                                                     "chain_coords_converter_slib.so")
        self.EXTRACT_SUBCHAIN_LIB = os.path.join(self.LOCATION, "modules", "extract_subchain_slib.so")
        self.CHAIN_FILTER_BY_ID = os.path.join(self.LOCATION, "modules", "chain_filter_by_id")

        self.CHAIN_BDB_INDEX = os.path.join(self.LOCATION, "modules", "chain_bdb_index.py")
        self.BED_BDB_INDEX = os.path.join(self.LOCATION, "modules", "bed_bdb_index.py")
        self.SPLIT_CHAIN_JOBS = os.path.join(self.LOCATION, "split_chain_jobs.py")
        self.MERGE_CHAINS_OUTPUT = os.path.join(self.LOCATION, "modules", "merge_chains_output.py")
        self.CLASSIFY_CHAINS = os.path.join(self.LOCATION, "modules", "classify_chains.py")
        self.SPLIT_EXON_REALIGN_JOBS = os.path.join(self.LOCATION, "split_exon_realign_jobs.py")
        self.MERGE_CESAR_OUTPUT = os.path.join(self.LOCATION, "modules", "merge_cesar_output.py")
        self.TRANSCRIPT_QUALITY = os.path.join(self.LOCATION, "modules", "get_transcripts_quality.py")
        self.GENE_LOSS_SUMMARY = os.path.join(self.LOCATION, "modules", "gene_losses_summary.py")
        self.ORTHOLOGY_TYPE_MAP = os.path.join(self.LOCATION, "modules", "orthology_type_map.py")
        self.MODEL_TRAINER = os.path.join(self.LOCATION, "train_model.py")
        self.DEFAULT_CESAR = os.path.join(self.LOCATION, "cesar")
        self.nextlow_rel_ = os.path.join(self.LOCATION, "execute_joblist.nf")
        self.NF_EXECUTE = os.path.abspath(self.nextlow_rel_)

    def __check_dependencies(self):
        """Check all dependencies."""
        eprint("check if binaries are compiled and libs are installed...")
        c_not_compiled = any(os.path.isfile(f) is False for f in [self.CHAIN_SCORE_FILTER,
                                                                  self.CHAIN_COORDS_CONVERT_LIB,
                                                                  self.CHAIN_FILTER_BY_ID,
                                                                  self.EXTRACT_SUBCHAIN_LIB])
        if c_not_compiled:
            eprint("Warning! C code is not compiled, trying to compile...")
        imports_not_found = False
        try:
            import twobitreader
            import networkx
            import pandas
            import xgboost
            import joblib
            import bsddb3
        except ImportError:
            eprint("Warning! Some of the required packages are not installed.")
            imports_not_found = True

        not_all_found = any([c_not_compiled, imports_not_found])
        self.__call_proc(self.CONFIGURE, "Could'd not call configure.sh!")\
            if not_all_found else eprint("All dependencies found")

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
        if self.nextflow_config_dir is None:
            # no nextflow config provided -> using local executor
            self.cesar_config_template = None
            self.nf_chain_extr_config_file = None
            self.local_executor = True
            return
        # check that required config files are here
        if not os.path.isdir(self.nextflow_config_dir):
            self.die(f"Error! Nextflow config dir {self.nextflow_config_dir} doesn't exist!")
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
        eprint("{0} in progress...".format(cmd))
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            eprint(extra_msg) if extra_msg else None
            self.die(f"Error! Process {cmd} died! Abort.")
        eprint(f"{cmd} done with code 0")

    def __find_two_bit(self, db):
        """Find a 2bit file."""
        # TODO: in public release must be simplified
        if os.path.isfile(db):
            return db
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
        self.__make_indexed_chain()
        self.__make_indexed_bed()
        self.__time_mark("Made indexes")

        # 1) split_chain_jobs.py - - make chain: genes joblist
        self.__split_chain_jobs()
        self.__time_mark("Splitted chain jobs")
        # 2) call "ch_cl_jobs" from previous stage as a bunch of cluster job
        self.__extract_chain_features()
        self.__time_mark("Chain jobs done")
        # 3) merge_chains_output.py  --> merge and rearrange data in "results folder"
        self.__merge_chains_output()
        self.__time_mark("Chains output merged")
        # 4) ./classify_chains.py
        # classify chains from raw output according the config file
        self.__classify_chains()
        self.__time_mark("Chains classified")
        # 5)./ split_exon_realign_jobs.py --> create cluster jobs for CESAR2.0
        self.__split_cesar_jobs()
        self.__time_mark("Split cesar jobs done")
        # 6) call combined cesar jobs
        self.__run_cesar_jobs()
        self.__time_mark("Dome cesar jobs")
        # 7) ./parse_cesar_bdb.py
        # extract fasta and bed files from CESAR2.0 output
        self.__merge_cesar_output()
        self.__time_mark("Merged cesar output")
        # 8). get pseudogene track
        self.__transcript_quality()
        self.__gene_loss_summary()
        self.__time_mark("Got gene loss summary")
        # 9) ./orthology_type_map.py
        self.__orthology_type_map()
        # 10) merge rejection reasons
        self.__merge_rejection_logs()
        # Everything is done
        self.__time_mark("Everything is done")
        self.die(f"Done! Estimated time: {dt.now() - self.t0}", rc=0)

    def __make_indexed_chain(self):
        """Make chain index file."""
        # make *.bb file
        eprint("make_indexed in progress...")
        chain_bdb_index(self.chain, self.chain_index_file)
        self.temp_files.append(self.chain_index_file)
        self.temp_files.append(self.chain)
        eprint("Indexed")

    def __time_mark(self, msg):
        """Left time mark."""
        if self.time_log is None:
            return
        t = dt.now() - self.t0
        with open(self.time_log, "a") as f:
            f.write(f"{msg} at {t}\n")

    def __make_indexed_bed(self):
        """Create gene_ID: bed line bdb indexed file."""
        eprint("index_bed in progress...")
        bed_bdb_index(self.ref_bed, self.index_bed_file)
        self.temp_files.append(self.index_bed_file)
        eprint("Indexed")

    def __split_chain_jobs(self):
        """Wrap split_jobs.py script."""
        # define arguments
        # save split jobs
        ch_cl_jobs = os.path.join(self.wd, "chain_classification_jobs")
        # for raw results of this stage
        self.chain_class_results = os.path.join(self.wd, "chain_classification_results")
        self.chain_cl_jobs_combined  = os.path.join(self.wd, "chain_class_jobs_combined")
        rejected_filename = "SPLIT_CHAIN_REJ.txt"
        rejected_path = os.path.join(self.rejected_dir,
                                     rejected_filename)
        self.temp_files.append(ch_cl_jobs)
        self.temp_files.append(self.chain_class_results)
        self.temp_files.append(self.chain_cl_jobs_combined)

        split_jobs_cmd = f"{self.SPLIT_CHAIN_JOBS} {self.chain} " \
                         f"{self.ref_bed} {self.index_bed_file} " \
                         f"--jobs_num {self.chain_jobs} " \
                         f"--jobs {ch_cl_jobs} " \
                         f"--jobs_file {self.chain_cl_jobs_combined} " \
                         f"--results_dir {self.chain_class_results} " \
                         f"--rejected {rejected_path}"

        self.__call_proc(split_jobs_cmd, "Could not split chain jobs!")
    
    def __extract_chain_features(self):
        """Execute extract chain features jobs."""
        nf_cmd = f"nextflow {self.NF_EXECUTE} " \
                 f"--joblist {self.chain_cl_jobs_combined}"
        if not self.local_executor:
            # not local executor -> provided config files
            # need abspath for nextflow execution
            nf_cmd += f" -c {self.nf_chain_extr_config_file}"
        # get timestamp to name the project and create a dir for that
        #  time() returns somrting like: 1595861493.8344169
        tmstmp = str(time.time()).split(".")[0]
        nf_project_name = f"{self.project_name}_chain_feats_at_{tmstmp}"
        nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
        os.mkdir(nf_project_path) if not os.path.isdir(nf_project_path) else None
        rc = subprocess.call(nf_cmd, shell=True, cwd=nf_project_path)
        if rc != 0:
            self.die(f"Error! Process {nf_cmd} died")

    def __merge_chains_output(self):
        """Call parse results."""
        # define where to save intermediate table
        eprint("Merging chain output...")        
        merge_chains_output(self.ref_bed, self.isoforms,
                            self.chain_class_results, self.chain_results_df)
        self.temp_files.append(self.chain_results_df)

    def __classify_chains(self):
        """Run decision tree."""
        # define input and output."""
        eprint("Decision tree in progress...")
        orthologs_file = "trans_to_chain_classes.tsv"
        pred_scores_file = "orthology_scores.tsv"
        self.orthologs = os.path.join(self.wd, orthologs_file)
        self.pred_scores = os.path.join(self.wd, pred_scores_file)
        self.se_model = os.path.join(self.LOCATION, "models", "se_model.dat")
        self.me_model = os.path.join(self.LOCATION, "models", "me_model.dat")
        cl_rej_log = os.path.join(self.rejected_dir, "classify_chains_rejected.txt")
        if not os.path.isfile(self.se_model) or not os.path.isfile(self.me_model):
            self.__call_proc(self.MODEL_TRAINER, "Could not train ML model!")
        classify_chains(self.chain_results_df, self.orthologs, self.se_model,
                        self.me_model, rejected=cl_rej_log, raw_out=self.pred_scores)
        if self.stop_at_chain_class:
            self.die("User requested to halt after chain features extraction", rc=0)

    def __split_cesar_jobs(self):
        """Call split_exon_realign_jobs.py."""
        if not self.t_2bit or not self.q_2bit:
            self.die("There is no 2 bit files provided, cannot go ahead and call CESAR.", rc=0)

        # if we call CESAR
        cesar_jobs_dir = os.path.join(self.wd, "cesar_jobs")
        self.cesar_combined = os.path.join(self.wd, "cesar_combined")
        self.cesar_results = os.path.join(self.wd, "cesar_results")
        self.temp_files.append(cesar_jobs_dir)
        self.temp_files.append(self.cesar_combined)

        # different field names depending on --ml flag
        fields = "ORTH,TRANS"

        self.temp_files.append(self.cesar_results)
        skipped_path = os.path.join(self.rejected_dir, "SPLIT_CESAR.txt")
        self.paralogs_log = os.path.join(self.wd, "paralogs.txt")

        split_cesar_cmd = f"{self.SPLIT_EXON_REALIGN_JOBS} "\
                          f"{self.orthologs} {self.ref_bed} " \
                          f"{self.index_bed_file} {self.chain_index_file} " \
                          f"{self.t_2bit} {self.q_2bit} " \
                          f"--jobs_dir {cesar_jobs_dir} " \
                          f"--jobs_num {self.cesar_jobs_num} " \
                          f"--combined {self.cesar_combined} " \
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
        self.__call_proc(split_cesar_cmd, "Could not split CESAR jobs!")
    
    def __run_cesar_jobs(self):
        """Run CESAR jobs using nextflow.
        
        At first -> push joblists, there might be a few of them
        At second -> monitor joblists, wait until all are done.
        """
        # for each bucket I create a separate joblist and config file
        # different config files because different memory limits
        processes = []  # keep subprocess objects here
        tmstmp = str(time.time()).split(".")[1]  # for project name
        # get a list of buckets
        if self.cesar_buckets == "0":
            buckets = [0, ]  # a single bucket
        else:  # several buckets, each int -> memory limit in gb
            buckets = [int(x) for x in self.cesar_buckets.split(",") if x != ""]
        print(f"Pushing {len(buckets)} joblists")
        # cmd to grep bucket-related commands
        grep_bucket_templ = "cat {0} | grep _{1}.bdb"
        for b in buckets:
            # create config file
            # 0 means that that buckets were not split
            memlim = b if b != 0 else self.cesar_mem_limit
            if not self.local_executor:
                # running on cluster, need to create config file
                # for this bucket's memory requirement
                config_string = self.cesar_config_template.replace("${_MEMORY_}", f"{memlim}")
                config_file_path = os.path.join(self.wd, f"cesar_config_{b}_queue.nf")
                config_file_abspath = os.path.abspath(config_file_path)
                with open(config_file_path, "w") as f:
                    f.write(config_string)
                self.temp_files.append(config_file_path)
            else:  # no config dir given: use local executor
                # OK if there is a single bucket
                config_file_abspath = None
            # extract jobs related to this bucket (if it's not 0)
            if b != 0:
                grep_bucket_cmd = grep_bucket_templ.format(self.cesar_combined, b)
                bucket_tasks = subprocess.check_output(grep_bucket_cmd, shell=True).decode("utf-8")
                joblist_name = f"cesar_joblist_queue_{b}.txt"
                joblist_path = os.path.join(self.wd, joblist_name)
                with open(joblist_path, "w") as f:
                    f.write(bucket_tasks)
                joblist_abspath = os.path.abspath(joblist_path)
                self.temp_files.append(joblist_path)
            else:  # nothing to extract, there is a single joblist
                joblist_abspath = os.path.abspath(self.cesar_combined)
            # create project directory for logs
            nf_project_name = f"{self.project_name}_cesar_at_{tmstmp}_q_{b}"
            nf_project_path = os.path.join(self.nextflow_dir, nf_project_name)
            os.mkdir(nf_project_path) if not os.path.isdir(nf_project_path) else None
            # create subprocess object
            nf_cmd = f"nextflow {self.NF_EXECUTE} " \
                     f"--joblist {joblist_abspath}"
            if config_file_abspath:
                nf_cmd += f" -c {config_file_abspath}"
            p = subprocess.Popen(nf_cmd, shell=True, cwd=nf_project_path)
            sys.stderr.write(f"Pushed cluster jobs with {nf_cmd}")
            processes.append(p)
            time.sleep(CESAR_PUSH_INTERVAL)

        # monitor jobs
        iter_num = 0
        while True:
            # Run until all jobs are done (or crashed)
            all_done = True  # default val, re-define if something is not done
            for p in processes:
                # check if each process is still running
                running = p.poll() is None
                if running:
                    all_done = False
            if all_done:
                print("All jobs succeeded!")
                break
            else:
                print(f"Iter {iter_num} waiting {ITER_DURATION * iter_num} seconds Not done")
                time.sleep(ITER_DURATION)
                iter_num += 1


    def __merge_cesar_output(self):
        """Merge CESAR output, save final fasta and bed."""
        eprint("Merging CESAR output to make final fasta and pre-final bed files.")
        merge_c_stage_skipped = os.path.join(self.rejected_dir,
                                             "CESAR_MERGE.txt")
        self.temp_files.append(self.intermediate_bed)

        merge_cesar_output(self.cesar_results,
                           self.intermediate_bed,
                           self.nucl_fasta,
                           self.meta_data,
                           merge_c_stage_skipped,
                           self.prot_fasta,
                           self.trash_exons)

    def __transcript_quality(self):
        """Call module to get transcript quality."""
        self.trans_quality_file = os.path.join(self.wd, "transcript_quality.tsv")
        classify_transcripts(self.meta_data,
                             self.pred_scores,
                             self.hq_orth_treshold,
                             self.trans_quality_file)

    def __gene_loss_summary(self):
        """Call gene loss summary."""
        eprint("Calling gene loss summary")
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
        eprint("Calling orthology_type_map...")
        skipped_ref_trans = os.path.join(self.wd, "ref_orphan_transcripts.txt")
        orthology_type_map(self.ref_bed,
                           self.query_annotation,
                           self.orthology_type,
                           ref_iso=self.isoforms,
                           que_iso=query_isoforms_file,
                           paralogs_arg=self.paralogs_log,
                           loss_data=self.loss_summ,
                           save_skipped=skipped_ref_trans)

    def __merge_rejection_logs(self):
        """Merge files containing data about rejected transcripts/genes."""
        rejected_log = os.path.join(self.wd,
                                    "genes_rejection_reason.tsv")

        # sanity checks
        if not os.path.isdir(self.rejected_dir):
            sys.exit(f"Error! {self.rejected_dir} is not a directory")
        log_files = os.listdir(self.rejected_dir)
        if len(self.rejected_dir) == 0:
            sys.exit(f"Error! {self.rejected_dir} is empty")

        buffer = open(rejected_log, "w")  # write everything there
        for log_file in log_files:
            log_path = os.path.join(self.rejected_dir, log_file)
            with open(log_path, "r") as f:
                content = [x for x in f.readlines() if x != "\n"]
            buffer.write("".join(content))
        buffer.close()


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_initial", type=str,
                     help="Chain file. Extensions like "
                          "FILE.chain or FILE.chain.gz are applicable.")
    app.add_argument("bed_initial", type=str,
                     help="Bed file with annotations for the target genome.")
    app.add_argument("tDB", default=None,
                     help="Alias of 2bit file for target genome.")
    app.add_argument("qDB", default=None,
                     help="Alias 2bit file for query genome.")
    # global ops
    app.add_argument("--project_folder", default=None,
                     help="Working directory. If not specified, "
                          "use CURRENT_DIR/PROJECT_NAME as default.")
    app.add_argument("--project_name", "--pn", default=None,
                     help="Project name, for example a name of species."
                          "If not set, tries to get the name from the chain file. "
                          "Create folder with this name.")
    app.add_argument("--min_score", "--msc", type=int, default=15000,
                     help="Filter the chains, preserve only those"
                          " with score bigger that this parameter. "
                          "15000 is default.")
    app.add_argument("--isoforms", "-i", type=str, default="",
                     help="Isoforms dictionary for parse_results.")
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
    # chain features related
    app.add_argument("--chain_jobs_num", "--chn", type=int, default=50,
                     help="Number of chain features extractor jobs.")
    app.add_argument("--no_chain_filter", "--ncf", action="store_true",
                     dest="no_chain_filter",
                     help="Do not filter the chain file (be sure you "
                          "specified a .chain but not .gz file in this case")
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
                     help="See ./split_exon_realign_jobs.py --help")
    app.add_argument("--cesar_chain_limit", type=int, default=100,
                     help="Ignore genes having more than N orthologs.")
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
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def main():
    """Entry point."""
    args = parse_args()
    toga_manager = Toga(args)
    toga_manager.run()


if __name__ == "__main__":
    main()
