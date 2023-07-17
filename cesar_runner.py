#!/usr/bin/env python3
"""Run a batch of CESAR jobs and save the output."""
import argparse
import os.path
import sys
import subprocess
from subprocess import PIPE
from modules.common import to_log
from modules.common import setup_logger
from version import __version__

__author__ = "Bogdan M. Kirilenko"

MAX_ATTEMPTS = 2
ZERO_CODE = 0
ERR_CODE = 1
FRAGM_CHAIN_ISSUE_CODE = 2

MODULE_NAME_FOR_LOG = "cesar_runner"


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("jobs_file", help="File containing a list of CESAR wrapper jobs")
    app.add_argument("output", help="BDB file containing CESAR wrapper output")
    app.add_argument(
        "--check_loss", default=None, help="File to save gene loss data if requested"
    )
    app.add_argument("--log_file", help="Main log file")
    app.add_argument("--rejected_log", default=None, help="Log gene rejection events")
    app.add_argument("--unproc_log", "--ul", default=None, help="Log unprocessed genes")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def call_job(cmd):
    """Call job, continue loop if fails."""
    attempts = 0
    # try 3 times
    err_msg = ""
    while attempts < MAX_ATTEMPTS:
        # cmd_out = subprocess.check_output(cmd, shell=True).decode("utf-8")
        p = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        b_stdout, b_stderr = p.communicate()
        rc = p.returncode
        cmd_out = b_stdout.decode("utf-8")
        err_msg = b_stderr.decode("utf-8").replace("\n", " ")
        if rc == ZERO_CODE:
            return cmd_out, ZERO_CODE
        elif rc == FRAGM_CHAIN_ISSUE_CODE:
            err_msg = f"CESAR_wrapper.py detected that fragments overlap for {cmd}, abort"
            return err_msg, FRAGM_CHAIN_ISSUE_CODE
        else:
            to_log(f"!!{err_msg}")
            to_log(f"!!FAILED COMMAND: {cmd}")
            to_log(f"!!At attempt {attempts}")
            attempts += 1
    to_log(f"!!FAILED TO EXECUTE COMMAND {cmd} in {MAX_ATTEMPTS} attempts")
    return err_msg, ERR_CODE  # send failure signal


def __job_to_transcript(job):
    """Extract transcript ID from job."""
    fields = job.split()
    return fields[1]


def main():
    """Entry point."""
    args = parse_args()
    setup_logger(args.log_file, write_to_console=False)
    filename = os.path.basename(args.jobs_file)
    # read jobs
    with open(args.jobs_file, "r") as f:
        # text file, a command per line
        jobs = [x.rstrip() for x in f.readlines()]
    jobs_num = len(jobs)
    log_prefix = f"{MODULE_NAME_FOR_LOG}::{filename}"
    to_log(f"{log_prefix}: started executing {jobs_num} jobs")
    unprocessed_genes = []

    out = open(args.output, "w")  # handle output file
    gene_loss_data = []  # list to keep gene loss detector out
    rejected = []  # keep genes that were skipped at this stage + reason

    for num, job in enumerate(jobs, 1):
        to_log(f"{log_prefix}: calling job {job}")
        # catch job stdout
        job_out, rc = call_job(job)
        to_log(f"{log_prefix}: return code: {rc}")
        if rc == FRAGM_CHAIN_ISSUE_CODE:
            # very special case -> nothig we can do
            # mark as missnig, I guess
            to_log(f"{log_prefix}: WARNING fragment chains overlap for {job}")
            rejected.append(f"{job}\tfragment chains overlap\n")
            unprocessed_genes.append(__job_to_transcript(job))
            continue
        if rc == 1:
            # a job failed with code 1 -> send the signal upstream
            # abort execution, write what job exactly failed
            # there are rare cases where CESAR fails
            # these cases usually contain rubbish
            to_log(f"{log_prefix}: CESAR JOB FAILED: {job}")
            to_log(f"{log_prefix}: ERROR MESSAGE: {job_out}")
            rejected.append(f"{job}\tCESAR JOB FAILURE\t{job_out}\n")
            continue

        # job looks like ./CESAR_wrapper.py
        # GENE CHAINS BEDFILE CHAINFILE TWOBIT FILES
        if job_out.startswith(">>>STOP_CODON>>>"):
            # CESAR wrapper detected in-frame stop codon in REFERENCE
            # There was not --mask_stops flag -> job cannot get proceed
            to_log(f"{log_prefix}: REFERENCE IN FRAME STOP CODON FOUND: {job}")
            rejected.append(f"{job}\tIN FRAME STOP CODON\n")
            continue

        gene = job.split()[1]  # just extract the gene from command

        if args.check_loss:
            # processing job output, there is CESAR out, meta-data + inact. mut
            job_out_lines = job_out.split("\n")
            # lines starting with # and ! -> inact mut scanner output
            job_gene_loss = "\n".join(
                [line for line in job_out_lines
                 if line.startswith("#") or line.startswith("!")]
            )
            gene_loss_data.append((gene, job_gene_loss))
            # all other lines -> processed CESAR output
            job_out = "\n".join(
                [
                    line
                    for line in job_out_lines
                    if not line.startswith("#") and not line.startswith("!")
                ]
            )

        # write output
        out.write(f"#{gene}\n")
        out.write(f"{job_out}\n")

    out.close()

    to_log(f"{log_prefix}: saving output for joblist")
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

    if args.unproc_log and len(unprocessed_genes) > 0:
        f = open(args.unproc_log, "w")
        for elem in unprocessed_genes:
            f.write(f"{elem}\n")
        f.close()


if __name__ == "__main__":
    main()
