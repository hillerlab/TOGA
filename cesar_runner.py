#!/usr/bin/env python3
"""Run a batch of CESAR jobs and save the output."""
import argparse
import sys
import subprocess

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(msg + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("jobs_file", help="File containing a list of CESAR wrapper jobs")
    app.add_argument("output", help="BDB file containing CESAR wrapper output")
    app.add_argument("--check_loss", default=None,
                     help="File to save gene loss data if requested")
    app.add_argument("--rejected_log", default=None,
                     help="Log gene rejection events")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def call_job(cmd):
    """Call job, continue loop if fails."""
    try:  # try to call this job
        cmd_out = subprocess.check_output(cmd, shell=True).decode("utf-8")
        return cmd_out
    except subprocess.CalledProcessError:
        eprint(f"{cmd} FAILED")
        return 1  # send failure signal


def main():
    """Entry point."""
    args = parse_args()
    # read jobs
    with open(args.jobs_file, "r") as f:
        # text file, a command per line
        jobs = [x.rstrip() for x in f.readlines()]
    jobs_num = len(jobs)

    out = open(args.output, "w")  # handle output file
    gene_loss_data = []  # list to keep gene loss detector out
    rejected = []  # keep genes that were skipped at this stage + reason

    for num, job in enumerate(jobs):
        eprint(f"Calling:\n{job}")
        # catch job stdout
        job_out = call_job(job)
        if job_out == 1:
            # a job failed with code 1 -> send the signal upstream
            # abort execution, write what job exactly failed
            sys.stderr.write(f"Error! Job {job} failed!\n")
            sys.exit(1)

        # job looks like ./CESAR_wrapper.py
        # GENE CHAINS BEDFILE CHAINFILE TWOBIT FILES
        if job_out.startswith(">>>STOP_CODON>>>"):
            # CESAR wrapper detected in-frame stop codon in REFERENCE
            # There was not --mask_stops flag -> job cannot get proceed
            rejected.append(f"{job}\tIN FRAME STOP CODON\n")
            continue

        gene = job.split()[1]  # just extract the gene from command

        if args.check_loss:
            # processing job output, there is CESAR out, meta-data + inact. mut
            job_out_lines = job_out.split("\n")
            # lines starting with # and ! -> inact mut scanner output
            job_gene_loss = "\n".join([line for line in job_out_lines
                                       if line.startswith("#")
                                       or line.startswith("!")])
            gene_loss_data.append((gene, job_gene_loss))
            # all other lines -> processed CESAR output
            job_out_cesar_lines = [line for line in job_out_lines
                                   if not line.startswith("#")
                                   and not line.startswith("!")]
            job_out = "\n".join(job_out_cesar_lines)

        # write output
        out.write(f"#{gene}\n")
        out.write(f"{job_out}\n")
        eprint(f"{num + 1} / {jobs_num} done", end="\r")

    out.close()

    if args.check_loss:
        # need to save gene loss data also
        if len(gene_loss_data) == 0:
            # this might happen
            gene_loss_data.append(("None", "No inactivating mutations detected"))
        f = open(args.check_loss, "w")
        for elem in gene_loss_data:
            # write gene: inact scanner report
            gene = elem[0]
            g_loss_report = elem[1]
            f.write(f"GENE: {gene}\n")
            f.write(g_loss_report)
            f.write("\n\n")  # separator for diff genes
        f.close()

    if args.rejected_log:
        # save list of unprocessed genes + reasons
        f = open(args.rejected_log, "w")
        f.write("".join(rejected))
        f.close()

    sys.exit(0)


if __name__ == "__main__":
    main()
