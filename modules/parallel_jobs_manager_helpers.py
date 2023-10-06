"""Helper function related to parallelisation."""
import time
import os
from constants import Constants
from jinja2 import FileSystemLoader, Environment
from modules.common import to_log

__author__ = "Bogdan M. Kirilenko"

ITER_DURATION = 60  # CESAR jobs check interval


def monitor_jobs(jobs_managers, die_if_sc_1=False):
    """Monitor parallel jobs if many batches run simultaneously."""
    to_log(f"## Stated polling cluster jobs until they done")
    iter_num = 0
    while True:  # Run until all jobs are done (or crashed)
        all_done = True  # default val, re-define if something is not done
        for job_manager in jobs_managers:
            # check if each process is still running
            rc = job_manager.check_status()
            if rc is None:
                all_done = False
        if all_done:
            to_log("### CESAR jobs done ###")
            break
        else:
            to_log(f"Polling iteration {iter_num}; already waiting {ITER_DURATION * iter_num} seconds.")
            time.sleep(ITER_DURATION)
            iter_num += 1

    if any(jm.return_code != 0 for jm in jobs_managers) and die_if_sc_1 is True:
        # some para/nextflow job died: critical issue
        # if die_if_sc_1 is True: terminate the program
        err = "Error! Some para/nextflow processes died!"
        # TODO: think about the best error class
        raise AssertionError(err)


def load_template(filename):
    templateLoader = FileSystemLoader(Constants.PARA_TEMPLATES_DIR)
    templateEnv = Environment(loader=templateLoader, autoescape=True, trim_blocks=True)
    return templateEnv.get_template(filename)
