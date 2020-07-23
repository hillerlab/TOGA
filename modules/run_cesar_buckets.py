#!/usr/bin/env python3
"""Use this module in case if you have several CESAR buckets."""
import argparse
import os
import sys
import subprocess
from subprocess import PIPE
import time
from datetime import datetime as dt

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

PUSH_INTERVAL = 20  # seconds; between each push
# ITER_LIMIT = 3000  # number of para check cycles
ITER_DURATION = 120  # seconds between each para check
kill = lambda process: process.kill()


def run_cesar_buckets(cesar_buckets, project_name, cesar_combined, wd):
    """Run CESAR jobs in several buckets."""
    t0 = dt.now()
    time_now = str(t0).split()[1].split(":")
    temp_files = []
    buckets = [int(x) for x in cesar_buckets.split(",") if x != ""]
    paras = []  # collect para project names
    processes = []  # subprocess objects
    proj_names = {b: "{0}_CESAR_at_{1}_{2}_b{3}".format(project_name, time_now[0], time_now[1], b)
                  for b in buckets}
    grep_bucket_templ = "cat {0} | grep _{1}.bdb"
    for b in buckets:  # call each bucket separately
        grep_bucket_cmd = grep_bucket_templ.format(cesar_combined, b)
        try:
            bucket_tasks = subprocess.check_output(grep_bucket_cmd, shell=True).decode("utf-8")
        except subprocess.CalledProcessError:
            continue  # trash
        if len(bucket_tasks) < 10:
            continue  # also nothing
        tasks_num = len(bucket_tasks.split("\n"))
        bucket_file = "cesar_comb_{0}_queue".format(b)
        bucket_path = os.path.join(wd, bucket_file)
        temp_files.append(bucket_path)
        with open(bucket_path, "w") as f:
            f.write(bucket_tasks)
        memoryMb = b * 1000
        para_proc = "para make {0} {1} -q=day -memoryMb={2} -maxNumResubmission 3"\
                    "".format(proj_names[b], bucket_path, memoryMb)
        paras.append(proj_names[b])
        p = subprocess.Popen(para_proc, shell=True)
        processes.append(p)
        sys.stderr.write("Pushed {0} cluster jobs in {1}\n".format(tasks_num, proj_names[b]))
        time.sleep(PUSH_INTERVAL)

    iter_num = 0
    while True:
        # YES, I know it is a bad idea to use such the constructions
        # run until all jobs are done
        all_done = True
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
    return temp_files

if __name__ == "__main__":
    # for tests
    app = argparse.ArgumentParser()
    app.add_argument("cesar_combined")
    app.add_argument("project_name")
    app.add_argument("cesar_buckets")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    project_name = args.project_name[:-1] if args.project_name.endswith("/") else args.project_name
    wd = os.path.join(os.getcwd(), project_name),
    _ = run_cesar_buckets(args.cesar_buckets, args.project_name, args.cesar_combined, wd)
