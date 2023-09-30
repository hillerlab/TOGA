#!/usr/bin/env python3
"""Strategy pattern implementation to handle parallel jobs.

Provides implementations for nextflow and para strategies.
Please feel free to implement your custom strategy if
neither nextflow nor para satisfy your needs.

WIP, to be enabled later.
"""
from abc import ABC, abstractmethod
import subprocess
import os
import shutil
from modules.common import to_log
from version import __version__

__author__ = "Bogdan M. Kirilenko"


class ParallelizationStrategy(ABC):
    """
    Abstract base class for a parallelization strategy.
    """
    def __init__(self):
        self._process = None

    @abstractmethod
    def execute(self, joblist_path, manager_data, label, wait=False, **kwargs):
        """
        Execute the jobs in parallel.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        :param wait: Boolean -> controls whether run blocking or not
        """
        pass

    @abstractmethod
    def check_status(self):
        """
        Check the status of the jobs.

        :return: Status of the jobs.
        """
        pass

    def terminate_process(self):
        """Terminates the associated process"""
        if self._process:
            self._process.terminate()


class NextflowStrategy(ParallelizationStrategy):
    """
    Concrete strategy for parallelization using Nextflow.
    """
    CESAR_CONFIG_TEMPLATE_FILENAME = "call_cesar_config_template.nf"
    CHAIN_CONFIG_TEMPLATE_FILENAME = "extract_chain_features_config.nf"
    CHAIN_JOBS_PREFIX = "chain_feats__"
    CESAR_JOBS_PREFIX = "cesar_jobs__"
    CESAR_CONFIG_MEM_TEMPLATE = "${_MEMORY_}"

    def __init__(self):
        super().__init__()
        self._process = None
        self.joblist_path = None
        self.manager_data = None
        self.label = None
        self.nf_project_path = None
        self.keep_logs = False
        self.use_local_executor = None
        self.nextflow_config_dir = None
        self.nextflow_logs_dir = None
        self.memory_limit = 16
        self.nf_master_script = None
        self.config_path = None
        self.return_code = None

    def execute(self, joblist_path, manager_data, label, wait=False, **kwargs):
        """Implementation for Nextflow."""
        # define parameters
        self.joblist_path = joblist_path
        self.manager_data = manager_data
        self.label = label
        self.memory_limit = int(kwargs.get("memory_limit", "16"))

        self.nf_project_path = manager_data.get("nextflow_dir", None)  # in fact, contains NF logs
        self.keep_logs = manager_data.get("keep_nf_logs", False)
        self.use_local_executor = manager_data.get("local_executor", False)
        self.nf_master_script = manager_data["NF_EXECUTE"]  # NF script that calls everything
        self.nextflow_config_dir = manager_data.get("nextflow_config_dir", None)
        self.config_path = self.__create_config_file()
        # create the nextflow process

        cmd = f"nextflow {self.nf_master_script} --joblist {joblist_path}"
        if self.config_path:
            cmd += f" -c {self.config_path}"

        log_dir = manager_data["logs_dir"]
        os.mkdir(log_dir) if not os.path.isdir(log_dir) else None
        log_file_path = os.path.join(manager_data["logs_dir"], f"{label}.log")
        with open(log_file_path, "w") as log_file:
            to_log(f"Parallel manager: pushing job {cmd}")
            self._process = subprocess.Popen(cmd,
                                             shell=True,
                                             stdout=log_file,
                                             stderr=log_file,
                                             cwd=self.nf_project_path)
        if wait:
            self._process.wait()

    def __create_config_file(self):
        """Create config file and return path to it if needed"""
        if self.use_local_executor:
            # for local executor, no config file is needed
            return None
        if self.label.startswith(self.CHAIN_JOBS_PREFIX):
            config_path = os.path.abspath(os.path.join(self.nextflow_config_dir,
                                                       self.CHAIN_CONFIG_TEMPLATE_FILENAME))
            return config_path
        if self.label.startswith(self.CESAR_JOBS_PREFIX):
            # need to craft CESAR joblist first
            config_template_path = os.path.abspath(os.path.join(self.nextflow_config_dir,
                                                                self.CESAR_CONFIG_TEMPLATE_FILENAME))
            with open(config_template_path, "r") as f:
                cesar_config_template = f.read()
            config_string = cesar_config_template.replace(self.CESAR_CONFIG_MEM_TEMPLATE,
                                                          f"{self.memory_limit}")
            config_filename = f"cesar_config_{self.memory_limit}_queue.nf"
            toga_temp_dir = self.manager_data["temp_wd"]
            config_path = os.path.abspath(os.path.join(toga_temp_dir, config_filename))
            with open(config_path, "w") as f:
                f.write(config_string)
            return config_path
        self.use_local_executor = True  # ??? should not be reachable normally
        return None  # using local executor again

    def check_status(self):
        """Check if nextflow jobs are done."""
        if self.return_code:
            return self.return_code
        running = self._process.poll() is None
        if running:
            return None
        self.return_code = self._process.returncode
        # the process just finished
        # nextflow provides a huge and complex tree of log files
        # remove them if user did not explicitly ask to keep them
        # if not self.keep_logs and self.nf_project_path:
        #     # remove nextflow intermediate files
        #     shutil.rmtree(self.nf_project_path) if os.path.isdir(self.nf_project_path) else None
        if self.config_path and self.label.startswith(self.CESAR_JOBS_PREFIX):
            # for cesar TOGA creates individual config files
            os.remove(self.config_path) if os.path.isfile(self.config_path) else None
        return self.return_code


class ParaStrategy(ParallelizationStrategy):
    """
    Concrete strategy for parallelization using Para.

    Para is rather an internal Hillerlab tool to manage slurm.

    """

    def __init__(self):
        super().__init__()
        self._process = None
        self.return_code = None

    def execute(self, joblist_path, manager_data, label, wait=False, **kwargs):
        """Implementation for Para."""
        cmd = f"para make {label} {joblist_path} "
        if "queue_name" in kwargs:
            queue_name = kwargs["queue_name"]
            cmd += f" -q={queue_name} "
        # otherwise use default medium queue
        if "memory_limit" in kwargs:
            memory_mb = kwargs["memory_limit"] * 1000  # para uses MB instead of GB
            cmd += f" --memoryMb={memory_mb}"
        # otherwise use default para's 10Gb

        log_dir = manager_data["logs_dir"]
        os.mkdir(log_dir) if not os.path.isdir(log_dir) else None
        log_file_path = os.path.join(manager_data["logs_dir"], f"{label}.log")
        with open(log_file_path, "w") as log_file:
            self._process = subprocess.Popen(cmd, shell=True, stdout=log_file, stderr=log_file)
        if wait:
            self._process.wait()

    def check_status(self):
        """Check if Para jobs are done."""
        if self.return_code:
            return self.return_code
        running = self._process.poll() is None
        if not running:
            self.return_code = self._process.returncode
            return self.return_code
        else:
            return None


class SnakeMakeStrategy(ParallelizationStrategy):
    """
    Not implemented class for Snakemake strategy.
    Might be helpful for users experiencing issues with Nextflow.
    """
    def __int__(self):
        self._process = None
        self.return_code = None
        raise NotImplementedError("Snakemake strategy is not yet implemented")

    def execute(self, joblist_path, manager_data, label, wait=False, **kwargs):
        raise NotImplementedError("Snakemake strategy is not yet implemented")

    def check_status(self):
        raise NotImplementedError("Snakemake strategy is not yet implemented")


class CustomStrategy(ParallelizationStrategy):
    """
    Custom parallel jobs execution strategy.
    """

    def __init__(self):
        super().__init__()
        self._process = None
        self.return_code = None
        raise NotImplementedError("Custom strategy is not implemented -> pls see documentation")

    def execute(self, joblist_path, manager_data, label, wait=False, **kwargs):
        """Custom implementation.

        Please provide your implementation of parallel jobs' executor.
        Jobs are stored in the joblist_path, manager_data is a dict
        containing project-wide TOGA parameters.

        The method should build a command that handles executing all the jobs
        stored in the file under joblist_path. The process object is to be
        stored in the self._process. It is recommended to create a non-blocking subprocess.

        I would recommend to store the logs in the manager_data["logs_dir"].
        Please have a look what "manager_data" dict stores -> essentially, this is a
        dump of the whole Toga class attributes.

        If your strategy works well, we can include it in the main repo.
        """
        raise NotImplementedError("Custom strategy is not implemented -> pls see documentation")

    def check_status(self):
        """Check if Para jobs are done.

        Please provide implementation of a method that checks
        whether all jobs are done.

        To work properly, the method should return None if the process is still going.
        Otherwise, return status code (int)."""
        raise NotImplementedError("Custom strategy is not implemented -> pls see documentation")


class ParallelJobsManager:
    """
    Class for managing parallel jobs using a specified parallelization strategy.
    """

    def __init__(self, strategy: ParallelizationStrategy):
        """
        Initialize the manager with a parallelization strategy.

        :param strategy: The parallelization strategy to use.
        """
        self.strategy = strategy
        self.return_code = None

    def execute_jobs(self, joblist_path, manager_data, label, **kwargs):
        """
        Execute jobs in parallel using the specified strategy.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        """
        self.strategy.execute(joblist_path, manager_data, label, **kwargs)

    def check_status(self):
        """
        Check the status of the jobs using the specified strategy.

        :return: Status of the jobs.
        """
        return self.strategy.check_status()

    def terminate_process(self):
        """Terminate associated process."""
        self.strategy.terminate_process()
