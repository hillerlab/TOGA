#!/usr/bin/env python3
"""Split jobs for the chain classification.

Creates files containing chain to genes lines.
You can create either one file to call the pipeline on a PC
or several files to call chain runner multiple times in parallel on cluster.
"""
import argparse
import os
import sys
import subprocess
import random
from datetime import datetime as dt
from modules.chain_bed_intersect import chain_bed_intersect
from modules.common import parts
from modules.common import die
from modules.common import setup_logger
from modules.common import to_log
from version import __version__

__author__ = "Bogdan Kirilenko, 2020."
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# global
WORK_DATA = {}  # script-related data
LOCATION = os.path.abspath(os.path.dirname(__file__))
CHAIN_RUNNER = os.path.join(LOCATION, "chain_runner.py")
t0 = dt.now()


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("chain_file", type=str, help="Chain file for local alignments.")
    # app.add_argument("chain_index", type=str, help="Chain sqlite 3 index db")
    app.add_argument("bed_file", type=str, help="Bed file, gene annotations.")
    app.add_argument("bed_index", type=str, help="Indexed bed")
    app.add_argument("--log_file", type=str, help="Path to logfile")
    app.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Don't print to console"
    )
    app.add_argument("--parallel_logs_dir", type=str, help="Path to dir storing logs from each cluster job")
    app.add_argument(
        "--jobs_num",
        "--jn",
        type=int,
        default=800,
        help="Number of cluster jobs, 800 as default.",
    )
    app.add_argument(
        "--job_size",
        type=int,
        default=None,
        help="How many jobs to put into one cluster job."
        "If defined, --jobs_Num is ignored.",
    )
    app.add_argument(
        "--jobs",
        "-j",
        type=str,
        default="chain_classification_jobs",
        help="Directory to save lists with chains and "
        "intersected genes, chain_classification_jobs as default.",
    )
    app.add_argument(
        "--jobs_file",
        "-c",
        type=str,
        default="jobs_file",
        help="File containing combined jobs, jobs_file as default.",
    )
    app.add_argument(
        "--results_dir",
        "-r",
        type=str,
        default="results",
        help="Redirect stdout from cluster job to this dir, " "results_dir as default",
    )
    app.add_argument(
        "--errors_dir",
        "-e",
        type=str,
        default=None,
        help="Redirect stderr from cluster job to " "this dir, None as default",
    )
    app.add_argument(
        "--make_index",
        "-i",
        action="store_true",
        dest="make_index",
        help="Make index file.",
    )
    # second-part related stuff
    app.add_argument(
        "--index_file",
        "-b",
        type=str,
        help="BDB file containing "
        "chains. If not assigned use [chain_file].bdb as default.",
    )
    app.add_argument(
        "--ref", type=str, default="hg38", help="Reference species, hg38 as default."
    )
    app.add_argument(
        "--vv", action="store_true", dest="vv", help="Add -v flag to unit commands."
    )
    app.add_argument(
        "--rejected", default=None, help="Track rejected genes in the file given"
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def call_proc(cmd):
    """Call a subprocess and catch errors."""
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        die(f"Error! Process {cmd} died! Abort.")


def check_args(args):
    """Check if args are correct, fill global dict."""
    WORK_DATA["vv"] = True if args.vv else False

    try:  # check the directories, create if it is necessary
        os.mkdir(args.jobs) if not os.path.isdir(args.jobs) else None
        os.mkdir(args.results_dir) if not os.path.isdir(args.results_dir) else None
        os.mkdir(args.errors_dir) if args.errors_dir and not os.path.isdir(
            args.errors_dir
        ) else None
        WORK_DATA["jobs"] = args.jobs
        WORK_DATA["results_dir"] = args.results_dir
        WORK_DATA["errors_dir"] = args.errors_dir

    except FileNotFoundError as grepexc:  # a one of those tasks failed
        to_log(f"Arguments are corrupted!\n{str(grepexc)}")
        die("Cannot create one of the directories requested.")

    # define about chain and bed files
    WORK_DATA["chain_file"] = (
        args.chain_file
        if os.path.isfile(args.chain_file)
        else die(f"Error! Chain file {args.chain_file} is wrong!")
    )

    WORK_DATA["bed_file"] = (
        args.bed_file
        if os.path.isfile(args.bed_file)
        else die(f"Error! Bed file {args.bed_file} is wrong!")
    )
    to_log(f"split_chain_jobs: Use bed file {args.bed_file} and chain file {args.chain_file}")

    # look for .ID.bb file
    index_file = (
        args.index_file
        if args.index_file
        else args.chain_file.replace(".chain", ".chain_ID_position")
    )

    if os.path.isfile(index_file):  # check if bb file is here
        WORK_DATA["index_file"] = index_file
    # elif args.make_index:  # create index if not exists
    #     idbb_cmd = f"/modules/chain_bdb_index.py {args.chain_file} {index_file}"
    #     call_proc(idbb_cmd)
    #     WORK_DATA["index_file"] = index_file
    else:  # die
        die(
            f"Error! Cannot find index file at {index_file}\n"
            "Please define it manually"
        )

    # define the number of jobs
    if args.job_size:  # easy:
        WORK_DATA["job_size"] = args.job_size
        WORK_DATA["jobs_num"] = None
    else:  # we must compute how many jobs to put into one cluster job
        WORK_DATA["job_size"] = None
        WORK_DATA["jobs_num"] = args.jobs_num
    WORK_DATA["bed_index"] = args.bed_index

    # some defaults
    WORK_DATA["jobs_file"] = args.jobs_file
    WORK_DATA["ref"] = args.ref
    # check if we are on cluster
    WORK_DATA["on_cluster"] = True
    to_log("split_chain jobs: the run data overview is:\n")
    for k, v in WORK_DATA.items():
        to_log(f"* {k}: {v}")


def get_chroms():
    """Read bed file and extract chroms."""
    f = open(WORK_DATA["bed_file"], "r")
    chroms_lst = []  # put ALL chroms here
    for line in f:  # I need the first field each line
        chroms_lst.append(line.split("\t")[0])
    f.close()  # return unique chroms
    return list(set(chroms_lst))


def get_intersections():
    """Make an array of intersections between genes and alignments."""
    to_log("split_chain_jobs: searching for intersections between reference transcripts and chains")
    # this function will get all chain X bed track intersections
    chain_genes_raw, skipped = chain_bed_intersect(
        WORK_DATA["chain_file"], WORK_DATA["bed_file"]
    )
    chain_genes = {k: ",".join(v) + "," for k, v in chain_genes_raw.items()}
    del chain_genes_raw
    # skipped genes do not intersect any chain
    # please note that chains with too low score are likely filtered out
    to_log(f"split_chain_jobs: chains-to-transcripts dict contains {len(chain_genes)} records")
    to_log(f"split_chain_jobs: skipped {len(skipped)} transcripts that do not intersect any chain")
    return chain_genes, skipped


def get_template():
    """Create a template for chain classification commands."""
    template = CHAIN_RUNNER + " {0} "
    # in case of using a nodes-associated disk I cannot use original filenames
    bed_to_template = WORK_DATA["bed_index"]
    bdb_to_template = WORK_DATA["chain_file"]
    template += " {0} {1}".format(bed_to_template, bdb_to_template)
    template += " -v" if WORK_DATA["vv"] else ""
    return template


def make_commands(intersection):
    """Shuffle the data and create joblist."""
    # order - list of chain IDS
    order = list(intersection.keys())
    # need to randomize the order of chains
    # most likely their list was sorted by chain ID
    # and the higher is chainID the longer it is
    # with randomization all jobs will require more or less the same amount of time
    # otherwise we get some incredibly long jobs and some too short
    random.shuffle(order)
    # fill the command list with chain ids and genes
    commands = [f"{c}\t{intersection.get(c)}" for c in order]
    return commands


def split_commands(commands):
    """Split the commands into N cluster jobs."""
    to_log(f"split_chain_jobs: preparing {len(commands)} commands")
    if WORK_DATA["job_size"]:  # size of cluster job is pre-defined
        job_size = WORK_DATA["job_size"]
    else:  # if was not defined - compute the size of cluster job
        job_size = (len(commands) // WORK_DATA["jobs_num"]) + 1
    to_log(f"split_chain_jobs: command size of {job_size} for each cluster job")
    batch = parts(commands, n=job_size)
    to_log(f"split_chain_jobs: results in {len(batch)} cluster jobs")
    return batch


def save_rejected_genes(skipped, filename):
    """Save rejected genes."""
    f = open(filename, "w")
    for line in skipped:
        f.write(f"{line[0]}\t{line[1]}\n")
    f.close()


def save(template, batch, logs_dir=None):
    """Save the cluster jobs, create jobs_file file."""
    filenames = {}  # collect filenames of cluster jobs
    for num, jobs in enumerate(batch):
        # define the path for the job
        job_path = os.path.join(WORK_DATA["jobs"], f"part_{num}")
        filenames[num] = job_path  # need these paths for jobs_file file
        # put the \n-separated jobs into the template
        # save this finally
        with open(job_path, "w") as f:
            f.write("\n".join(jobs) + "\n")

    # save the jobs_file file
    # add > line for stdout and 2> for stderr (if required)
    f = open(WORK_DATA["jobs_file"], "w")
    for num, path in filenames.items():
        cmd = template.format(path)
        stdout_part = f"> {WORK_DATA['results_dir']}/{num}.txt"
        if logs_dir:
            logs_part = f" --log_file {logs_dir}/chain_runner_{num}.log"
        else:
            logs_part = ""
        jobs_file_line = f"{cmd} {logs_part} {stdout_part}\n"
        f.write(jobs_file_line)
    f.close()


def main():
    """Entry point."""
    args = parse_args()
    setup_logger(args.log_file, write_to_console=not args.quiet)
    check_args(args)  # check if all the files, dependencies etc are correct
    intersections, skipped = get_intersections()  # intersect chains and beds
    # # extract genes that are not intersected by any chain
    # # or chrom is not aligned at all
    # missing_genes_not_ali = [x[0] for x in skipped]
    if args.rejected:
        # skipped: genes that do not intersect any chain
        save_rejected_genes(skipped, args.rejected)
    commands = make_commands(intersections)  # shuffle and create set of commands
    batch = split_commands(commands)  # split the commands into cluster jobs
    template = get_template()
    save(template, batch, logs_dir=args.parallel_logs_dir)  # save jobs and a jobs_file file
    to_log("split_chain_jobs: estimated time: {0}".format(dt.now() - t0))
    sys.exit(0)


if __name__ == "__main__":
    main()
