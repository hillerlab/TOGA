#!/usr/bin/env python3
"""Create CESAR joblist.

According to predicted orthologous chains create CESAR jobs.
Merge them into joblists.
"""
import argparse
import os
import sys
import math
from collections import defaultdict
from datetime import datetime as dt
from re import finditer, IGNORECASE
import ctypes
from twobitreader import TwoBitFile
from modules.common import parts
from modules.common import split_in_n_lists
from modules.common import chain_extract_id
from modules.common import eprint
from modules.common import make_cds_track
from modules.common import die

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# 0 gene; 1 chains; 2 bed_file; 3 bdb_chain_file; 4 tDB; 5 qDB; 6 memlim gig;
LOCATION = os.path.dirname(__file__)
WRAPPER_ABSPATH = os.path.abspath(os.path.join(LOCATION, "CESAR_wrapper.py"))
WRAPPER_TEMPLATE = (
    WRAPPER_ABSPATH
    + " {0} {1} {2} {3} {4} {5} --cesar_binary {6}"
    + " --uhq_flank {7} --memlim {8}"
)
CESAR_RUNNER = os.path.abspath(
    os.path.join(LOCATION, "cesar_runner.py")
)  # script that will run jobs
LONG_LOCI_FIELDS = {
    "GGLOB",
    "TRANS",
}  # chain classes that could lead to very long query loci
BIGMEM_LIM = 1000  # mem limit for bigmem partition
REL_LENGTH_THR = 50
ABS_LENGTH_TRH = 1000000
EXTRA_MEM = 100000  # extra memory "just in case"
BIGMEM_JOBSNUM = 100  # TODO: make a parameter?
REF_LEN_THRESHOLD = 0.05  # if query length < 5% CDS then skip it

ASM_GAP_SIZE = 10
ASM_GAP_PATTERN = r"N{" + str(ASM_GAP_SIZE) + ",}"

M = "M"
L = "L"

ORTHOLOG = "ORTH"
PARALOG = "PARA"
TRANS = "TRANS"
PROJECTION = "PROJECTION"
TRANSCRIPT = "TRANSCRIPT"

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


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("orthologs_file", help="Output of the chain classifier.")
    app.add_argument("bed_file", type=str, help="BED FILE")
    app.add_argument("bdb_bed_file", type=str, help="BDB BED FILE")
    app.add_argument("bdb_chain_file", type=str, help="BDB CHAIN FILE")
    app.add_argument("tDB", type=str, help="target 2 bit")
    app.add_argument("qDB", type=str, help="query 2 bit")

    app.add_argument(
        "--cesar_binary",
        type=str,
        default="cesar",
        help="CESAR2.0 binary address, cesar as default.",
    )
    app.add_argument(
        "--jobs_num",
        type=int,
        default=300,
        help="Total number of cluster jobs, 300 is recommended."
        " Resulting number may slightly vary in case of buckets "
        "param usage due to round issues.",
    )
    app.add_argument(
        "--buckets",
        default="0",
        help=""
        "If you need to split the cluster jobs in different classes"
        " according the memory consumption use this parameter. To do "
        " that write comma-separated list of memory levels. For "
        "example, --buckets 10,30 means that there are two classes of "
        "jobs - consuming 10 and 30 gb. All jobs consuming more than 30gb "
        "are ignored. Job names will be 'cesar_job_[job_number]_[memory_class]' "
        "like cesar_job_9999_30 - meaning all tasks in this file require "
        "no more that 30Gb. --buckets 0 means no separation.",
    )
    app.add_argument(
        "--mask_stops",
        "--ms",
        action="store_true",
        dest="mask_stops",
        help="Mask stop codons in target sequences. CESAR cannot process them."
        "Using this parameter please make sure you know what you are doing.",
    )
    app.add_argument(
        "--chains_limit",
        type=int,
        default=15,
        help="Skip genes with amount of orthologs more than the limit.",
    )
    app.add_argument(
        "--skipped_genes",
        default=None,
        help="If a gene was skipped due to memory of number "
        " of chain limit, save it into a file.",
    )
    app.add_argument(
        "--mem_limit",
        type=float,
        default=50,
        help="Skip genes requiring more than X GB to call CESAR",
    )
    app.add_argument("--jobs_dir", default="cesar_jobs", help="Save jobs in.")
    app.add_argument(
        "--combined", default="cesar_combined", help="Combined cluster jobs."
    )
    app.add_argument("--bigmem", default="cesar_bigmem", help="CESAR bigmem joblist")
    app.add_argument("--results", default="cesar_results", help="Save results to.")
    app.add_argument(
        "--check_loss", default=None, help="Call internal gene loss pipeline"
    )
    app.add_argument("--u12", default=None, help="Add U12 introns data")
    app.add_argument(
        "--rejected_log", default=None, help="Save rejection data in this dir"
    )
    app.add_argument(
        "--paralogs_log",
        default=os.path.join(os.path.dirname(__file__), "paralogs.log"),
        help="Write a list of genes for which only paralogous chains were detected.",
    )
    app.add_argument("--uhq_flank", default=50, type=int, help="UHQ flank size")
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
        "--annotate_paralogs",
        "--ap",
        action="store_true",
        dest="annotate_paralogs",
        help="Annotate paralogs instead of orthologs.",
    )
    app.add_argument(
        "--fragments_data", help="Gene: fragments file for fragmented genomes."
    )
    app.add_argument(
        "--opt_cesar",
        action="store_true",
        dest="opt_cesar",
        help="Using lastz-optimized version of CESAR",
    )
    app.add_argument(
        "--precomp_memory_data",
        default=None,
        help="Memory consumption was already precomputed",
    )
    app.add_argument(
        "--predefined_glp_class_path",
        default=None,
        help="Save preliminary projection classification for: "
        "(i) Projections with too short query region (L or M) and "
        "(ii) Projections with very long query region (M)",
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_u12_data(u12_data_file):
    """Read U12 introns."""
    u12_data = defaultdict(list)
    if not u12_data_file:
        # not provided
        return u12_data
    f = open(u12_data_file, "r")
    f.__next__()
    for line in f:
        line_data = line[:-1].split("\t")
        trans = line_data[0]
        exon_num = int(line_data[1])
        site = line_data[2]
        val = (exon_num, site)
        u12_data[trans].append(val)
    f.close()
    return u12_data


def define_buckets(lim, buckets):
    """Return memory limit in Gig if required. Get classes."""
    if buckets == "0":
        # split was not required
        return lim, {0: []}
    # buckets assigned
    buckets_values = sorted([int(x) for x in buckets.split(",") if x != ""])
    buckets = {x: [] for x in buckets_values}
    lim = buckets_values[-1]
    return lim, buckets


def read_orthologs(orthologs_file, only_o2o=False, annotate_paralogs=False):
    """Read orthologs file."""
    # convert fields param string to list
    # fields = [x.upper() for x in fields_raw.split(",") if x != ""]
    genes_chains = {}
    chain_gene_field = {}
    skipped = []  # genes skipped at this stage
    _no_chains_intersecting = []
    f = open(orthologs_file, "r")  # open the file
    f.__next__()  # skip header
    # first column: transcript identifier
    # then: chain class fields (like column 2 - orthologous chains, 3 - paralogous)
    for line in f:
        # parse line
        line_info = line[:-1].split("\t")
        # "0" is a placeholder meaning "no chains there"
        gene = line_info[0]
        selected, chains = [], {}

        chains[ORTHOLOG] = [x for x in line_info[1].split(",") if x != "0"]
        chains[PARALOG] = [x for x in line_info[2].split(",") if x != "0"]
        chains[TRANS] = [x for x in line_info[3].split(",") if x != "0"]
        # Processed pseudogenes column ignored -> they are processed separately
        all_chains = chains[ORTHOLOG] + chains[PARALOG] + chains[TRANS]

        if len(all_chains) == 0:
            # no way in running CESAR on this gene
            # because there are no chains we could use
            skipped.append((gene, "0", "No chains intersecting the gene"))
            _no_chains_intersecting.append(gene)
            continue

        # user can ask to process only the genes that have a single orthologous chain
        # here we check that this is the case
        not_one2one = len(chains[ORTHOLOG]) == 0 or len(chains[ORTHOLOG]) > 1
        if only_o2o and not_one2one:  # we requested only a single orthologous chain
            skipped.append((gene, "0", "Only one2one requested, this gene didn't pass"))
            continue

        # use orthologous chains by default,
        # if no orthologous chains -> use spanning chains (TRANS)
        # no spanning chains -> use paralogous
        if annotate_paralogs:
            selected_field = PARALOG
        elif len(chains[ORTHOLOG]) > 0:
            selected_field = ORTHOLOG
        elif len(chains[TRANS]) > 0:
            selected_field = TRANS
        else:
            selected_field = PARALOG

        selected = chains[selected_field].copy()
        # mark used field
        for chain in selected:
            key = (chain, gene)
            chain_gene_field[key] = selected_field

        # write to the dict, gene to chains we will use
        genes_chains[gene] = selected

    f.close()
    die(
        "Error! No gene:chains pairs selected! Probably --fields parameter is wrong!"
    ) if len(genes_chains) == 0 else None
    return genes_chains, chain_gene_field, skipped, _no_chains_intersecting


def read_bed(bed):
    """Read bed 12 file.

    For each transcript extract genetic coordinates and exon sizes.
    """
    bed_data = {}
    f = open(bed, "r")
    for line in f:
        cds_track = make_cds_track(line).split("\t")
        bed_info = line[:-1].split("\t")
        chrom = bed_info[0]
        chrom_start = int(bed_info[1])
        chrom_end = int(bed_info[2])
        name = bed_info[3]
        block_sizes = [int(x) for x in cds_track[10].split(",") if x != ""]
        bed_data[name] = (chrom, chrom_start, chrom_end, block_sizes)
    f.close()
    return bed_data


def define_short_q_proj_stat(q_chrom, q_start, q_end, q_2bit):
    """If the query sequence is shorter than 50% CDS,

    we need to check whether projection is really deleted (Lost)
    or is just missing due to assembly gaps.
    """
    query_genome_sequence = TwoBitFile(q_2bit)
    query_chrom = query_genome_sequence[q_chrom]
    query_seq = query_chrom[q_start:q_end].upper()
    # N are two-sided?
    # simply copy pasted solution from CESAR wrapper.py
    # can be further optimised
    gap_ranges = 0
    for match in finditer(ASM_GAP_PATTERN, query_seq, IGNORECASE):
        span_start, span_end = match.span()
        # gap_ranges.append((seq_start + span_start, seq_start + span_end))
        gap_ranges += 1
    if gap_ranges == 0:
        # no assembly gaps: really deleted -> Lost
        return L
    else:
        # there are assembly gaps -> Missing
        return M


def precompute_regions(
    batch, bed_data, bdb_chain_file, chain_gene_field, limit, q_2bit
):
    """Precompute region for each chain: bed pair."""
    eprint("Precompute regions for each gene:chain pair...")
    chain_to_genes, skipped = defaultdict(list), []
    predef_glp = {}  # predefined GLP classification

    # revert the dict, from gene2chain to chain2genes
    for gene, chains_not_sorted in batch.items():
        if len(chains_not_sorted) == 0:
            skipped.append((gene, "no orthologous chains"))
            continue
        chains = sorted(chains_not_sorted, key=lambda x: int(x))
        chains = chains[:limit]

        if len(chains_not_sorted) > limit:
            # skip genes that have > limit orthologous chains
            chains_skipped = chains[limit:]
            skipped.append(
                (
                    gene,
                    ",".join(chains_skipped),
                    f"number of chains ({limit} chains) limit exceeded",
                )
            )
            # add each projection individually
            # further append to GLP classification
            for c_ in chains_skipped:
                proj_id = f"{gene}.{c_}"
                predef_glp[proj_id] = f"{PROJECTION}\t{M}"

        for chain_id in chains:
            chain_to_genes[chain_id].append(gene)
    # read regions themselves
    gene_chain_grange = defaultdict(dict)
    chains_num, iter_num = len(chain_to_genes.keys()), 0

    for chain_id, genes in chain_to_genes.items():
        # extract chain itself
        chain_body = chain_extract_id(bdb_chain_file, chain_id).encode()
        all_gene_ranges = []
        genes_cds_length = []
        for gene in genes:
            # get genomic coordinates for each gene
            gene_data = bed_data.get(gene)
            grange = f"{gene_data[0]}:{gene_data[1]}-{gene_data[2]}"
            cds_length = sum(gene_data[3])
            genes_cds_length.append(cds_length)
            all_gene_ranges.append(grange)

        # we need to get corresponding regions in the query
        # for now we have chain blocks coordinates and gene
        # regions in the reference genome
        # use chain_coords_converter shared library to
        # convert target -> query coordinates via chain
        # first need to convert to C-types
        c_chain = ctypes.c_char_p(chain_body)
        c_shift = ctypes.c_int(2)
        granges_bytes = [s.encode("utf-8") for s in all_gene_ranges]
        granges_num = len(all_gene_ranges)
        c_granges_num = ctypes.c_int(granges_num)
        granges_arr = (ctypes.c_char_p * (granges_num + 1))()
        granges_arr[:-1] = granges_bytes
        granges_arr[granges_num] = None

        # then call the function
        raw_ch_conv_out = ch_lib.chain_coords_converter(
            c_chain, c_shift, c_granges_num, granges_arr
        )
        chain_coords_conv_out = []  # keep lines here
        # convert C output to python-readable type
        for i in range(granges_num + 1):
            chain_coords_conv_out.append(raw_ch_conv_out[i].decode("utf-8"))

        for line in chain_coords_conv_out[1:]:
            # then parse the output
            # line contains information about transcript range in the query
            # and the corresponding locus in the reference
            line_info = line.rstrip().split()
            # line info is: region num, region in reference, region in query
            # one line per one gene, in the same order
            num = int(line_info[0])
            # regions format is chrom:start-end
            q_chrom = line_info[1].split(":")[0]
            q_grange = line_info[1].split(":")[1].split("-")
            q_start, q_end = int(q_grange[0]), int(q_grange[1])
            que_len = q_end - q_start
            t_grange = line_info[2].split(":")[1].split("-")
            t_start, t_end = int(t_grange[0]), int(t_grange[1])
            tar_len = t_end - t_start
            len_delta = abs(tar_len - que_len)
            delta_gene_times = len_delta / tar_len
            gene = genes[num]  # shared lib returns data per gene in the same order
            proj_id = f"{gene}.{chain_id}"
            cds_length = genes_cds_length[num]
            min_query_length = cds_length * REF_LEN_THRESHOLD
            field = chain_gene_field.get((chain_id, gene))
            # check that corresponding region in the query is not too long
            # for instance query locus is 50 times longer than the gene
            # or it's longer than 1M base and also this is a TRANS chain
            high_rel_len = delta_gene_times > REL_LENGTH_THR
            high_abs_len = len_delta > ABS_LENGTH_TRH
            long_loci_field = field in LONG_LOCI_FIELDS
            if (high_rel_len or high_abs_len) and long_loci_field:
                skipped.append((gene, chain_id, "too long query locus"))
                # print(f"TOO LONG: {proj_id}")
                predef_glp[proj_id] = f"{PROJECTION}\t{M}"
                continue
            # in contrast, if query locus is too short (<5% CDS length)
            # then CESAR might not build HMM properly, we skip this
            # hard to imagine in what case such an input will give us any meaningful result
            if que_len < min_query_length:
                # in this case we need to check whether the gene is truly deleted
                # in the corresponding locus or is missing
                # to separate these cases, TOGA checks whether the region contains
                # assembly gaps
                # skipped.append((gene, chain_id, "too short query locus"))
                proj_stat = define_short_q_proj_stat(q_chrom, q_start, q_end, q_2bit)
                predef_glp[proj_id] = f"{PROJECTION}\t{proj_stat}"
                # print(f"Too short query: {proj_id}: {que_len}: {proj_stat}")
                continue
            # for each chain-gene pair save query region length
            # need this for required memory estimation
            gene_chain_grange[gene][chain_id] = que_len

        del raw_ch_conv_out  # not sure if necessary but...
        iter_num += 1  # verbosity
        eprint(f"Chain {iter_num} / {chains_num}", end="\r")
    return gene_chain_grange, skipped, predef_glp


def fill_buckets(buckets, all_jobs):
    """Split jobs in buckets according their memory consumption."""
    if 0 in buckets.keys():  # do not split it
        buckets[0] = list(all_jobs.keys())
        return buckets
    # buckets were set
    memlims = sorted(buckets.keys())
    prev_lim = 0
    for memlim in memlims:
        # buckets and memory limits are pretty much the same
        # if buckets are 5 and 10 then:
        # memlim[5] -> jobs that require <= 5Gb
        # memlim[10] -> jobs that require > 5Gb AND <= 10Gb
        buckets[memlim] = [
            job for job, job_mem in all_jobs.items() if prev_lim < job_mem <= memlim
        ]
        prev_lim = memlim
    # remove empty
    filter_buckets = {k: v for k, v in buckets.items() if len(v) > 0}
    return filter_buckets


def save_jobs(filled_buckets, bucket_jobs_num, jobs_dir):
    """Save cesar calls in the dir assigned."""
    os.mkdir(jobs_dir) if not os.path.isdir(jobs_dir) else None
    file_num, to_combine = 0, []
    for bucket_id, jobs in filled_buckets.items():
        num_of_files = bucket_jobs_num[bucket_id]
        # just in case
        num_of_files = len(jobs) if num_of_files >= len(jobs) else num_of_files
        size_of_file = len(jobs) // num_of_files
        # size_of_file = size_of_file + 1 if len(jobs) % num_of_files != 0 else size_of_file
        jobs_split = parts(jobs, n=size_of_file)
        for part in jobs_split:
            file_num += 1
            file_name = f"cesar_job_{file_num}_{bucket_id}"
            file_path = os.path.abspath(os.path.join(jobs_dir, file_name))
            f = open(file_path, "w")
            f.write("\n".join(part) + "\n")
            f.close()
            to_combine.append(file_path)
    return to_combine


def save_bigmem_jobs(bigmem_joblist, jobs_dir):
    """Save bigmem jobs."""
    # TODO: try to merge with save_jobs() func
    # one bigmem job per joblist, but not more than 100
    # if > 100: something is wrong
    joblist_size = len(bigmem_joblist)
    num_of_parts = joblist_size if joblist_size <= BIGMEM_LIM else BIGMEM_JOBSNUM
    if num_of_parts == 0:
        return None  # no bigmem jobs

    bigmem_parts = split_in_n_lists(bigmem_joblist, num_of_parts)
    bigmem_files_num = len(bigmem_parts)  # in case if num of jobs < BIGMEM_JOBSNUM
    bigmem_paths = []
    if bigmem_files_num == 0:
        return None  # no bigmem jobs at all
    for num, bigmem_part in enumerate(bigmem_parts):
        file_name = f"cesar_job_{num}_bigmem"
        file_path = os.path.abspath(os.path.join(jobs_dir, file_name))
        f = open(file_path, "w")
        f.write("\n".join(bigmem_part) + "\n")
        f.close()
        bigmem_paths.append(file_path)
    return bigmem_paths


def save_combined_joblist(
    to_combine, combined_file, results_dir, inact_mut_dat, rejected_log, name=""
):
    """Save joblist of joblists (combined joblist)."""
    f = open(combined_file, "w")
    for num, comb in enumerate(to_combine, 1):
        basename = os.path.basename(comb).split(".")[0]
        results_path = os.path.abspath(os.path.join(results_dir, basename + ".txt"))
        combined_command = f"{CESAR_RUNNER} {comb} {results_path}"
        if inact_mut_dat:
            loss_data_path = os.path.join(inact_mut_dat, f"{basename}.inact_mut.txt")
            combined_command += f" --check_loss {loss_data_path}"
        if rejected_log:
            log_path = os.path.join(rejected_log, f"{basename}.txt")
            combined_command += f" --rejected_log {log_path}"
        f.write(combined_command + "\n")
    f.close()


def read_fragments_data(in_file):
    """Read gene: fragments file."""
    ret = {}
    f = open(in_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        chain_str = line_data[1]
        # chains = [int(x) for x in line_data.split(",") if x != ""]
        # actually there are strings:
        chains = [x for x in chain_str.split(",") if x != ""]
        ret[gene] = chains
    f.close()
    return ret


def read_precomp_mem(precomp_file):
    """Read precomputed memory if exists."""
    ret = {}
    if precomp_file is None:
        return ret
    f = open(precomp_file, "r")
    for line in f:
        if line == "\n":
            continue
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        mem_raw = float(line_data[1])
        mem = math.ceil(mem_raw) + 1.25
        ret[gene] = mem
    f.close()
    return ret


def main():
    """Entry point."""
    t0 = dt.now()
    args = parse_args()
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash

    # as default we create CESAR jobs for chains with "orth" or "trans" class
    # but user could select another set of chain classes

    # read U12 introns: to create a list of U12-containing genes
    # need it to make subsequent commands
    u12_data = read_u12_data(args.u12)

    # if memory is precomputed: use it
    precomp_mem = read_precomp_mem(args.precomp_memory_data)
    # get lists of orthologous chains per each gene
    # skipped_1 - no chains found -> log them
    predefined_glp_class = {}  # for projections which are M and L without CESAR
    # m_ -> to be added to Missing bucket
    batch, chain_gene_field, skipped_1, m_ = read_orthologs(
        args.orthologs_file, only_o2o=args.o2o_only, annotate_paralogs=args.annotate_paralogs
    )
    for gene in m_:  # classify transcripts with no intersecting chains as missing
        predefined_glp_class[gene] = f"{TRANSCRIPT}\t{M}"
    # split cesar jobs in different buckets (if user requested so)
    # like put all jobs that require < 5Gig in the bucket 1
    # jobs requiring 5 to 15 Gb to bucket 2 and so on
    # CESAR might be very memory-consuming -> so we care about this
    mem_limit, buckets = define_buckets(args.mem_limit, args.buckets)

    # load reference bed file data; coordinates and exon sizes
    bed_data = read_bed(args.bed_file)
    # check if cesar binary exists
    die(
        f"Error! Cannot find cesar executable at {args.cesar_binary}!"
    ) if not os.path.isfile(args.cesar_binary) else None

    # if this is a fragmmented genome: we need to change CESAR commands for
    # split genes
    if args.fragments_data:
        gene_fragments_dict = read_fragments_data(args.fragments_data)
    else:  # better to create empty dict and call dict.get()
        gene_fragments_dict = dict()
    # pre-compute chain : gene : region data
    # collect the second list of skipped genes
    # skipped_2 -> too long corresponding regions in query
    regions, skipped_2, predef_glp = precompute_regions(
        batch,
        bed_data,
        args.bdb_chain_file,
        chain_gene_field,
        args.chains_limit,
        args.qDB,
    )
    predefined_glp_class.update(predef_glp)

    # start making the jobs
    all_jobs = {}
    skipped_3 = []
    bigmem_jobs = []

    for gene in batch.keys():
        u12_this_gene = u12_data.get(gene)
        block_sizes = bed_data[gene][3]

        gene_chains_data = regions.get(gene)
        # check that there is something for this gene
        if not gene_chains_data:
            continue
        elif len(gene_chains_data.keys()) == 0:
            continue

        gene_fragments = gene_fragments_dict.get(gene, False)
        if gene_fragments:
            # this is a fragmented gene, need to change the procedure a bit
            gene_chains_data = {
                k: v for k, v in gene_chains_data.items() if k in gene_fragments
            }
        chains = gene_chains_data.keys()
        if len(chains) == 0:
            continue
        chains_arg = ",".join(chains)  # chain ids -> one of the cmd args

        # if memory is precomputed then use it
        precomp_gig = precomp_mem.get(gene, None)
        if precomp_gig is None:
            # proceed to memory estimation
            # the same procedure as inside CESAR2.0 code
            num_states, r_length = 0, 0

            # required memory depends on numerous params
            # first, we need reference transcript-related parameters
            # query-related parameters will be later
            for block_size in block_sizes:
                # num_states += 6 + 6 * reference->num_codons + 1 + 2 + 2 + 22 + 6;
                #  /* 22 and 6 for acc and donor states */
                num_codons = block_size // 3
                num_states += 6 + 6 * num_codons + 1 + 2 + 2 + 22 + 6
                # r_length += 11 + 6 * fasta.references[i]->length
                # + donors[i]->length + acceptors[i]->length;
                r_length += block_size

            # now compute query sequence-related parameters
            query_lens = [v for v in gene_chains_data.values()]
            if (
                gene_fragments
            ):  # in case of fragmented genome: we stitch queries together
                # so query length = sum of all queries
                q_length_max = sum(query_lens)
            else:  # not fragmented genome: processins queries separately
                # thus we need only the max length
                q_length_max = max(query_lens)
            # and now compute the amount of required memory
            memory = (
                (num_states * 4 * 8)
                + (num_states * q_length_max * 4)
                + (num_states * 304)
                + (2 * q_length_max + r_length) * 8
                + (q_length_max + r_length) * 2 * 1
                + EXTRA_MEM
            )
            gig = math.ceil(memory / 1000000000) + 0.25
        else:
            # memory was precomputed
            gig = precomp_gig

        # gig = compute_amount_of_memory(block_sizes, q_length_max, args.opt_cesar)
        # # 0 gene; 1 chains; 2 bed_file; 3 bdb chain_file; 4 tDB; 5 qDB; 6 output; 7 cesar_bin
        job = WRAPPER_TEMPLATE.format(
            gene,
            chains_arg,
            os.path.abspath(args.bdb_bed_file),
            os.path.abspath(args.bdb_chain_file),
            os.path.abspath(args.tDB),
            os.path.abspath(args.qDB),
            os.path.abspath(args.cesar_binary),
            args.uhq_flank,
            gig,
        )
        # add some flags if required
        job = job + " --mask_stops" if args.mask_stops else job
        job = job + " --check_loss" if args.check_loss else job
        job = job + " --no_fpi" if args.no_fpi else job
        job = job + " --fragments" if gene_fragments else job
        job = job + " --opt_cesar" if args.opt_cesar else job

        # add U12 introns data if this gene has them:
        job = job + f" --u12 {os.path.abspath(args.u12)}" if u12_this_gene else job

        # define whether it's an ordinary or a bigmem job
        # depending on the memory requirements
        if gig <= mem_limit:  # ordinary job
            all_jobs[job] = gig
        elif gig <= BIGMEM_LIM:
            skipped_3.append((gene, ",".join(chains), f"requires {gig}) -> bigmem job"))
            predef_glp[gene] = f"{TRANSCRIPT}\tM"
            bigmem_jobs.append(job)
        else:
            skipped_3.append(
                (
                    gene,
                    ",".join(chains),
                    f"big mem limit ({BIGMEM_LIM} gig) exceeded (needs {gig})",
                )
            )
            predef_glp[gene] = f"{TRANSCRIPT}\tM"

    eprint(f"\nThere are {len(all_jobs.keys())} jobs in total.")
    eprint("Splitting the jobs.")
    # split jobs in buckets | compute proportions
    filled_buckets = fill_buckets(buckets, all_jobs)
    prop_sum = sum([k * len(v) for k, v in filled_buckets.items()])
    # estimate proportion of a bucket in the runtime
    buckets_prop = (
        {k: (k * len(v)) / prop_sum for k, v in filled_buckets.items()}
        if 0 not in filled_buckets.keys()
        else {0: 1.0}
    )
    eprint("Bucket proportions are:")
    eprint("\n".join([f"{k} -> {v}" for k, v in buckets_prop.items()]))
    eprint(f"Also there are {len(bigmem_jobs)} bigmem jobs")
    # get number of jobs for each bucket
    bucket_jobs_num = {k: math.ceil(args.jobs_num * v) for k, v in buckets_prop.items()}
    # save jobs, get comb lines
    to_combine = save_jobs(filled_buckets, bucket_jobs_num, args.jobs_dir)
    # save combined jobs, combined is a file containing paths to separate jobs
    os.mkdir(args.results) if not os.path.isdir(args.results) else None
    os.mkdir(args.check_loss) if args.check_loss and not os.path.isdir(
        args.check_loss
    ) else None

    # save joblist of joblists
    save_combined_joblist(
        to_combine, args.combined, args.results, args.check_loss, args.rejected_log
    )

    # save bigmem jobs, a bit different logic
    bigmem_paths = save_bigmem_jobs(bigmem_jobs, args.jobs_dir)
    if bigmem_paths:
        save_combined_joblist(
            bigmem_paths,
            args.bigmem,
            args.results,
            args.check_loss,
            args.rejected_log,
            name="bigmem",
        )

    # save skipped genes if required
    if args.skipped_genes:
        skipped = skipped_1 + skipped_2 + skipped_3
        f = open(args.skipped_genes, "w")
        # usually we have gene + reason why skipped
        # we split them with tab
        f.write("\n".join(["\t".join(x) for x in skipped]) + "\n")
        f.close()

    if args.predefined_glp_class_path:
        # if we know GLP class for some of the projections: save it
        f = open(args.predefined_glp_class_path, "w")
        for k, v in predefined_glp_class.items():
            f.write(f"{k}\t{v}\n")
        f.close()

    # save IDs of paralogous projections
    # skip if we annotate only paralogs
    if not args.annotate_paralogs:
        f = open(args.paralogs_log, "w")
        for k, v in chain_gene_field.items():
            if v != "PARA":
                continue
            gene_ = f"{k[1]}.{k[0]}\n"
            f.write(gene_)
        f.close()

    eprint(f"Estimated: {dt.now() - t0}")
    sys.exit(0)


if __name__ == "__main__":
    main()
