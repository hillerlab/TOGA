"""Utility class for Toga manager."""
import os
import shutil
import sys
from datetime import datetime as dt

from modules.common import to_log


class TogaUtil:
    @staticmethod
    def generate_project_name():
        """Generate project name automatically."""
        today_and_now = dt.now().strftime("%Y.%m.%d_at_%H:%M:%S")
        project_name = f"TOGA_project_on_{today_and_now}"
        return project_name

    @staticmethod
    def log_python_version():
        to_log(f"# python interpreter path: {sys.executable}")
        to_log(f"# python interpreter version: {sys.version}")

    @staticmethod
    def append_technical_err_to_predefined_class(transcripts_path, out_path):
        """Append file with predefined classifications."""
        if not os.path.isfile(transcripts_path):
            # in this case, we don't have transcripts with tech error
            # can simply quit the function
            return
        with open(transcripts_path, "r") as f:
            transcripts_list = [x.rstrip() for x in f]
        f = open(out_path, "a")
        for elem in transcripts_list:
            f.write(f"TRANSCRIPT\t{elem}\tM\n")
        f.close()

    @staticmethod
    def merge_directory_content(dir_name, output, ignore_empty=False):
        """Merge all files in a directory into one."""
        files_list = os.listdir(dir_name)
        if len(files_list) == 0 and ignore_empty is False:
            sys.exit(f"Error! {dir_name} is empty")
        elif len(files_list) == 0 and ignore_empty is True:
            # in this case we allow empty directories
            # just remove the directory and return
            shutil.rmtree(dir_name)
            return
        buffer = open(output, "w")
        for filename in files_list:
            path = os.path.join(dir_name, filename)
            with open(path, "r") as f:
                content = f.read()
            buffer.write(content)
        buffer.close()

    @staticmethod
    def terminate_parallel_processes(jobs_managers):
        to_log(f"KeyboardInterrupt: terminating {len(jobs_managers)} running parallel processes")
        for job_manager in jobs_managers:
            job_manager.terminate_process()