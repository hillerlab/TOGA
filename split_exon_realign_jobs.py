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
import ctypes
from modules.common import parts
from modules.common import chain_extract_id
from modules.common import eprint
from modules.common import die

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# 0 gene; 1 chains; 2 bed_file; 3 bdb_chain_file; 4 tDB; 5 qDB; 6 memlim gig;
LOCATION = os.path.dirname(__file__)
WRAPPER_ABSPATH = os.path.abspath(os.path.join(LOCATION, "CESAR_wrapper.py"))
WRAPPER_TEMPLATE = WRAPPER_ABSPATH \
                   + " {0} {1} {2} {3} {4} {5} --memlim {6} --cesar_binary {7}" \
                   + " --uhq_flank {8}"
CESAR_RUNNER = os.path.abspath(os.path.join(LOCATION, "cesar_runner.py"))  # script that will run jobs
LONG_LOCI_FIELDS = {"GGLOB", "TRANS"}  # chain classes that could lead to very long query loci
REL_LENGTH_THR = 50
ABS_LENGTH_TRH = 500000
EXTRA_MEM = 100000  # extra memory "just in case"

# connect shared lib; define input and output data types
chain_coords_conv_lib_path = os.path.join(LOCATION,
                                          "modules",
                                          "chain_coords_converter_slib.so")

ch_lib = ctypes.CDLL(chain_coords_conv_lib_path)
ch_lib.chain_coords_converter.argtypes = [ctypes.c_char_p,
                                          ctypes.c_int,
                                          ctypes.c_int,
                                          ctypes.POINTER(ctypes.c_char_p)]
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

    app.add_argument("--cesar_binary", type=str, default="cesar",
                     help="CESAR2.0 binary address, cesar as default.")
    app.add_argument("--jobs_num", type=int, default=300,
                     help="Total number of cluster jobs, 300 is recommended."
                     " Resulting number may slightly vary in case of buckets "
                     "param usage due to round issues.")
    app.add_argument("--buckets", default="0", help=""
                     "If you need to split the cluster jobs in different classes"
                     " according the memory consumption use this parameter. To do "
                     " that write comma-separated list of memory levels. For "
                     "example, --buckets 10,30 means that there are two classes of "
                     "jobs - consuming 10 and 30 gb. All jobs consuming more than 30gb "
                     "are ignored. Job names will be 'cesar_job_[job_number]_[memory_class]' "
                     "like cesar_job_9999_30 - meaning all tasks in this file require "
                     "no more that 30Gb. --buckets 0 means no separation.")
    app.add_argument("--fields", default=None,
                     help="Use those chains that are placed in these fields "
                          " in orthologs file. Comma-separated list. For example "
                          "PERF,GLOK - for perfect and good local chains.")
    app.add_argument("--mask_stops", "--ms", action="store_true", dest="mask_stops",
                     help="Mask stop codons in target sequences. CESAR cannot process them."
                     "Using this parameter please make sure you know what you are doing.")
    app.add_argument("--chains_limit", type=int, default=15,
                     help="Skip genes with amount of orthologs more than the limit.")
    app.add_argument("--skipped_genes", default=None,
                     help="If a gene was skipped due to memory of number "
                          " of chain limit, save it into a file.")
    app.add_argument("--mem_limit", type=float, default=50,
                     help="Skip genes requiring more than X GB to call CESAR")
    app.add_argument("--jobs_dir", default="cesar_jobs", help="Save jobs in.")
    app.add_argument("--combined", default="cesar_combined", help="Combined cluster jobs.")
    app.add_argument("--results", default="cesar_results", help="Save results to.")
    app.add_argument("--check_loss", default=None, help="Call internal gene loss pipeline")
    app.add_argument("--u12", default=None, help="Add U12 introns data")
    app.add_argument("--rejected_log", default=None, help="Save rejection data in this dir")
    app.add_argument("--paralogs_log", default=os.path.join(os.path.dirname(__file__), "paralogs.log"), 
                     help="Write a list of genes for which only paralogous chains were detected.")
    app.add_argument("--uhq_flank", default=50, type=int, help="UHQ flank size")
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


def read_orthologs(orthologs_file, fields_raw, only_o2o=False):
    """Read orthologs file."""
    # convert fields param string to list
    fields = [x.upper() for x in fields_raw.split(",") if x != ""]
    genes_chains = {}
    chain_gene_field = {}
    skipped = []  # genes skipped at this stage
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

        chains["ORTH"] = [x for x in line_info[1].split(",") if x != "0"]
        chains["PARA"] = [x for x in line_info[2].split(",") if x != "0"]
        chains["TRANS"] = [x for x in line_info[3].split(",") if x != "0"]
        # Processed pseudogenes column ignored -> they are processed separately
        all_chains = chains["ORTH"] + chains["PARA"] + chains["TRANS"]

        if len(all_chains) == 0:
            # no way in running CESAR on this gene
            # because there are no chains we could use
            skipped.append((gene, "0", "No chains intersecting the gene"))
            continue

        # user can ask to process only the genes that have a single orthologous chain
        # here we check that this is the case
        not_one2one = len(chains["ORTH"]) == 0 or len(chains["ORTH"]) > 1
        if only_o2o and not_one2one:  # we requested only a single orthologous chain
            skipped.append((gene, "0", "Only one2one requested, this gene didn't pass"))
            continue

        # get those are chosen in FIELDS
        for field in fields:
            # field is most likely "ORTH" or "TRANS"
            field_chains = chains.get(field)
            if not field_chains:
                continue
            selected.extend(field_chains)
            for chain in field_chains:
                key = (chain, gene)
                chain_gene_field[key] = field

        # if a gene has no orthologous chains, then use paralogous
        # if no paralogous -> log this gene
        if not selected:
            # no orthologous chains
            # we try to use paralogous chains
            # of course, log this data
            selected = all_chains.copy()
            keys = [(chain, gene) for chain in selected]
            for key in keys:
                chain_gene_field[key] = "PARALOG"
        # write to the dict, gene to chains we will use
        genes_chains[gene] = selected

    f.close()
    die("Error! No gene:chains pairs selected! Probably --fields parameter is wrong!") \
        if len(genes_chains) == 0 else None
    return genes_chains, chain_gene_field, skipped


def read_bed(bed):
    """Read bed 12 file.
    
    For each transcript extract genetic coordinates and exon sizes.
    """
    bed_data = {}
    f = open(bed, "r")
    for line in f:
        bed_info = line[:-1].split("\t")
        chrom = bed_info[0]
        chrom_start = int(bed_info[1])
        chrom_end = int(bed_info[2])
        name = bed_info[3]
        block_sizes = [int(x) for x in bed_info[10].split(',') if x != '']
        bed_data[name] = (chrom, chrom_start, chrom_end, block_sizes)
    f.close()
    return bed_data


def precompute_regions(batch, bed_data, bdb_chain_file, chain_gene_field, limit):
    """Precompute region for each chain: bed pair."""
    eprint("Precompute regions for each gene:chain pair...")
    chain_to_genes, skipped = defaultdict(list), []
    # revert the dict, from gene2chain to chain2genes
    for gene, chains in batch.items():
        if len(chains) == 0:
            skipped.append((gene, ",".join(chains), "no orthologous chains"))
            continue
        chains_ = sorted(chains, key=lambda x: int(x))
        chains_ = chains_[:limit]
        if len(chains) > limit:
            # skip genes that have > limit orthologous chains
            skipped.append((gene, ",".join(chains_[limit:]),
                            f"number of chains ({limit} chains) limit exceeded"))
        for chain in chains_:
            chain_to_genes[chain].append(gene)
    # read regions themselves
    gene_chain_grange = defaultdict(dict)
    chains_num, iter_num = len(chain_to_genes.keys()), 0

    for chain_id, genes in chain_to_genes.items():
        # extract chain itself
        chain_body = chain_extract_id(bdb_chain_file, chain_id).encode()
        all_gene_ranges = []
        for gene in genes:
            # get genomic coordinates for each gene
            gene_data = bed_data.get(gene)
            grange = f"{gene_data[0]}:{gene_data[1]}-{gene_data[2]}"
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
        raw_ch_conv_out = ch_lib.chain_coords_converter(c_chain,
                                                        c_shift,
                                                        c_granges_num,
                                                        granges_arr)
        chain_coords_conv_out = []  # keep lines here
        # convert C output to python-readable type
        for i in range(granges_num + 1):
            chain_coords_conv_out.append(raw_ch_conv_out[i].decode("utf-8"))

        for line in chain_coords_conv_out[1:]:
            # then parse the output
            line_info = line.rstrip().split()
            # line info is: region num, region in reference, region in query
            # one line per one gene, in the same order
            num = int(line_info[0])
            # regions format is chrom:start-end
            q_grange = line_info[1].split(":")[1].split("-")
            q_start, q_end = int(q_grange[0]), int(q_grange[1])
            que_len = q_end - q_start
            t_grange = line_info[2].split(":")[1].split("-")
            t_start, t_end = int(t_grange[0]), int(t_grange[1])
            tar_len = t_end - t_start
            len_delta = abs(tar_len - que_len)
            delta_gene_times = len_delta / tar_len
            gene = genes[num]
            field = chain_gene_field.get((chain_id, gene))
            # check that corresponding region in the query is not too long
            # if so: skip this
            high_rel_len = delta_gene_times > REL_LENGTH_THR
            high_abs_len = len_delta > ABS_LENGTH_TRH
            long_loci_field = field in LONG_LOCI_FIELDS
            if (high_rel_len or high_abs_len) and long_loci_field:
                skipped.append((gene, chain_id, "too long query locus"))
                continue
            # for each chain-gene pair save query region length
            # need this for required memory estimation
            gene_chain_grange[gene][chain_id] = que_len

        del raw_ch_conv_out  # not sure if necessary but...
        iter_num += 1  # verbosity
        eprint(f"Chain {iter_num} / {chains_num}", end="\r")
    return gene_chain_grange, skipped


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
        buckets[memlim] = [job for job, job_mem in all_jobs.items() if prev_lim < job_mem <= memlim]
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


def main():
    """Entry point."""
    t0 = dt.now()
    args = parse_args()
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise it could crash

    # as default we create CESAR jobs for chains with "orth" or "trans" class
    # but user could select another set of chain classes
    fields = "ORTH,TRANS" if args.fields is None else args.fields

    # read U12 introns: to create a list of U12-containing genes
    # need it to make subsequent commands
    u12_data = read_u12_data(args.u12)

    # get lists of orthologous chains per each gene
    # skipped_1 - no chains found -> log them
    batch, chain_gene_field, skipped_1 = read_orthologs(args.orthologs_file,
                                                        fields,
                                                        only_o2o=args.o2o_only)
    # split cesar jobs in different buckets (if user requested so)
    # like put all jobs that require < 5Gig in the bucket 1
    # jobs requiring 5 to 15 Gb to bucket 2 and so on
    # CESAR might be very memory-consuming -> so we care about this
    mem_limit, buckets = define_buckets(args.mem_limit, args.buckets)

    # load reference bed file data; coordinates and exon sizes
    bed_data = read_bed(args.bed_file)
    # check if cesar binary exists
    die(f"Error! Cannot find cesar executable at {args.cesar_binary}!") if \
        not os.path.isfile(args.cesar_binary) else None

    # pre-compute chain : gene : region data
    # collect the second list of skipped genes
    # skipped_2 -> too long corresponding regions in query
    regions, skipped_2 = precompute_regions(batch,
                                            bed_data,
                                            args.bdb_chain_file,
                                            chain_gene_field,
                                            args.chains_limit)
    
    # start making the jobs
    all_jobs = {}
    skipped_3 = []

    for gene in batch.keys():
        u12_this_gene = u12_data.get(gene)
        block_sizes = bed_data[gene][3]

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

        gene_chains_data = regions.get(gene)
        # check that there is something for this gene
        if not gene_chains_data:
            continue
        elif len(gene_chains_data) == 0:
            continue

        chains = gene_chains_data.keys()
        chains_arg = ",".join(chains)  # chain ids -> one of the cmd args
        
        # now compute query sequence-related parameters
        query_lens = [v for v in gene_chains_data.values()]
        q_length_max = max(query_lens)
        # and now compute the amount of required memory
        memory = (num_states * 4 * 8) + \
                 (num_states * q_length_max * 4) + \
                 (num_states * 304) + \
                 (2 * q_length_max + r_length) * 8 + \
                 (q_length_max + r_length) * 2 * 1 + EXTRA_MEM

        # convert to gigs + 0.25 extra gig
        gig = math.ceil(memory / 1000000000) + 0.25 
        if gig > mem_limit:
            # it is going to consume TOO much memory
            # skip this gene -> save to log
            skipped_3.append((gene, ",".join(chains),
                             f"memory limit ({mem_limit} gig) exceeded (needs {gig})"))
            continue

        # # 0 gene; 1 chains; 2 bed_file; 3 bdb chain_file; 4 tDB; 5 qDB; 6 output; 7 cesar_bin
        job = WRAPPER_TEMPLATE.format(gene, chains_arg,
                                      os.path.abspath(args.bdb_bed_file),
                                      os.path.abspath(args.bdb_chain_file),
                                      os.path.abspath(args.tDB),
                                      os.path.abspath(args.qDB),
                                      gig,
                                      os.path.abspath(args.cesar_binary),
                                      args.uhq_flank)
        # add some flags if required
        job = job + " --mask_stops" if args.mask_stops else job
        job = job + " --check_loss" if args.check_loss else job
        job = job + " --no_fpi" if args.no_fpi else job

        # add U12 introns data if this gene has them:
        job = job + f" --u12 {os.path.abspath(args.u12)}" if u12_this_gene else job

        all_jobs[job] = gig

    eprint(f"\nThere are {len(all_jobs.keys())} jobs in total.")
    eprint("Splitting the jobs.")
    # split jobs in buckets | compute proportions
    filled_buckets = fill_buckets(buckets, all_jobs)
    prop_sum = sum([k * len(v) for k, v in filled_buckets.items()])
    # estimate proportion of a bucket in the runtime
    buckets_prop = {k: (k * len(v)) / prop_sum for k, v in filled_buckets.items()} \
        if 0 not in filled_buckets.keys() else {0: 1.0}
    eprint("Bucket proportions are:")
    eprint("\n".join([f"{k} -> {v}" for k, v in buckets_prop.items()]))
    # get number of jobs for each bucket
    bucket_jobs_num = {k: math.ceil(args.jobs_num * v) for k, v in buckets_prop.items()}
    # save jobs, get comb lines
    to_combine = save_jobs(filled_buckets, bucket_jobs_num, args.jobs_dir)
    # save combined jobs, combined is a file containing paths to separate jobs
    os.mkdir(args.results) if not os.path.isdir(args.results) else None
    os.mkdir(args.check_loss) if args.check_loss \
        and not os.path.isdir(args.check_loss) else None

    f = open(args.combined, "w")
    for num, comb in enumerate(to_combine, 1):
        basename = os.path.basename(comb).split(".")[0]
        results_path = os.path.abspath(os.path.join(args.results, basename + ".bdb"))
        combined_command = f"{CESAR_RUNNER} {comb} {results_path}"
        if args.check_loss:
            loss_data_path = os.path.join(args.check_loss,
                                          f"{basename}.inact_mut.txt")
            combined_command += f" --check_loss {loss_data_path}"
        if args.rejected_log:
            log_path = os.path.join(args.rejected_log, f"{num}.txt")
            combined_command += f" --rejected_log {log_path}"
        f.write(combined_command + "\n")
    f.close()

    # save skipped genes if required
    if args.skipped_genes:
        skipped = skipped_1 + skipped_2 + skipped_3
        f = open(args.skipped_genes, "w")
        # usually we have gene + reason why skipped
        # we split them with tab
        f.write("\n".join(["\t".join(x) for x in skipped]) + "\n")
        f.close()

    f = open(args.paralogs_log, "w")
    # save IDs of paralogous projections
    for k, v in chain_gene_field.items():
        if v != "PARALOG":
            continue
        gene_ = f"{k[1]}.{k[0]}\n"
        f.write(gene_)
    f.close()

    eprint(f"Estimated: {dt.now() - t0}")
    sys.exit(0)


if __name__ == "__main__":
    main()
