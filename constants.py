"""Class holding all project-wide constants."""
import os

__author__ = "Bogdan M. Kirilenko"

# TODO: complete and organise the class
# TODO: think about splitting into subclasses


class Constants:
    LOCATION = os.path.dirname(__file__)

    U12_FILE_COLS = 3
    U12_AD_FIELD = {"A", "D"}
    ISOFORMS_FILE_COLS = 2
    NF_DIR_NAME = "nextflow_logs"
    NEXTFLOW = "nextflow"
    CESAR_PUSH_INTERVAL = 30  # CESAR jobs push interval
    ITER_DURATION = 60  # CESAR jobs check interval
    MEMLIM_ARG = "--memlim"
    FRAGM_ARG = "--fragments"

    CESAR_RUNNER = os.path.abspath(os.path.join(LOCATION, "cesar_runner.py"))
    CESAR_RUNNER_TMP = "{0} {1} {2} --check_loss {3} --rejected_log {4}"
    CESAR_PRECOMPUTED_REGIONS_DIRNAME = "cesar_precomputed_regions"
    CESAR_PRECOMPUTED_MEMORY_DIRNAME = "cesar_precomputed_memory"
    CESAR_PRECOMPUTED_ORTHO_LOCI_DIRNAME = "cesar_precomputed_orthologous_loci"

    CESAR_PRECOMPUTED_MEMORY_DATA = "cesar_precomputed_memory.tsv"
    CESAR_PRECOMPUTED_REGIONS_DATA = "cesar_precomputed_regions.tsv"
    CESAR_PRECOMPUTED_ORTHO_LOCI_DATA = "cesar_precomputed_orthologous_loci.tsv"

    NUM_CESAR_MEM_PRECOMP_JOBS = 500
    PARA_STRATEGIES = ["nextflow", "para", "custom"]  # TODO: add snakemake

    TEMP_CHAIN_CLASS = "temp_chain_trans_class"
    MODULES_DIR = "modules"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"

    # lists of features required by single and multi exon models
    SE_MODEL_FEATURES = ["gl_exo", "flank_cov", "exon_perc", "synt_log"]
    ME_MODEL_FEATURES = ["gl_exo", "loc_exo", "flank_cov", "synt_log", "intr_perc"]
