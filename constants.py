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
    GNU_PARALLEL = "parallel"
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
    PARA_STRATEGIES = ["nextflow", "para", "uge", "custom"]  # TODO: add snakemake
    CONFIG_STRATEGIES = ["uge"] # strategies that need JSON configuration
    PARA_SCHEMA_DIR = os.path.abspath(os.path.join(LOCATION, "parallel_config/schema"))
    PARA_TEMPLATES_DIR = os.path.abspath(os.path.join(LOCATION, "parallel_config/templates"))
    TEMP_CHAIN_CLASS = "temp_chain_trans_class"
    MODULES_DIR = "modules"
    RUNNING = "RUNNING"
    CRASHED = "CRASHED"
    TEMP = "temp"

    # lists of features required by single and multi exon models
    SE_MODEL_FEATURES = ["gl_exo", "flank_cov", "exon_perc", "synt_log"]
    ME_MODEL_FEATURES = ["gl_exo", "loc_exo", "flank_cov", "synt_log", "intr_perc"]

    # from CESAR_wrapper.py #
    FRAGMENT_CHAIN_ID = -1
    ORTH_LOC_LINE_SUFFIX = "#ORTHLOC"
    UNDEF_REGION = "None:0-0"

    # Sequence related #
    ATG_CODON = "ATG"
    XXX_CODON = "XXX"
    GAP_CODON = "---"
    NNN_CODON = "NNN"
    STOP_CODONS = {"TAG", "TGA", "TAA"}

    ACCEPTOR_SITE = ("ag",)
    DONOR_SITE = (
        "gt",
        "gc",
    )


class ConstColors:
    BLUE = "0,0,200"
    LIGHT_BLUE = "0,200,255"
    LIGHT_RED = "255,50,50"
    SALMON = "255,160,120"
    GREY = "130,130,130"
    BROWN = "159,129,112"
    BLACK = "10,10,10"


class InactMutClassesConst:
    MISS_EXON = "Missing exon"
    DEL_EXON = "Deleted exon"
    DEL_MISS = {MISS_EXON, DEL_EXON}
    COMPENSATION = "COMPENSATION"
    SSM = "SSM"
    # (ag)acceptor-EXON-donor(gt)
    SSM_D = "SSMD"  # Donor, right, GT,GC
    SSM_A = "SSMA"  # Acceptor, left, AG


# Standalone constants #
COMPLEMENT_BASE = {
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


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "---": "-",
    "NNN": "X",
}
