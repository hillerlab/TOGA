#!/usr/bin/env python3
"""Strategy pattern implementation to handle parallel jobs.

Provides implementations for nextflow and para strategies.
Please feel free to implement your custom strategy if
neither nextflow nor para satisfy your needs.
"""
from abc import ABC, abstractmethod
import subprocess
import os


class ParallelizationStrategy(ABC):
    """
    Abstract base class for a parallelization strategy.
    """

    @abstractmethod
    def execute(self, joblist_path, manager_data, label, **kwargs):
        """
        Execute the jobs in parallel.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        """
        pass

    @abstractmethod
    def check_status(self):
        """
        Check the status of the jobs.

        :return: Status of the jobs.
        """
        pass


class NextflowStrategy(ParallelizationStrategy):
    """
    Concrete strategy for parallelization using Nextflow.
    """

    def execute(self, joblist_path, manager_data, label, **kwargs):
        """Implementation for Nextflow."""
        pass

    def check_status(self):
        """Check if nextflow jobs are done."""
        pass


class ParaStrategy(ParallelizationStrategy):
    """
    Concrete strategy for parallelization using Para.
    """

    def __init__(self):
        self._process = None

    def execute(self, joblist_path, manager_data, label, **kwargs):
        """Implementation for Para."""
        cmd = f"para make {label} {joblist_path} "
        if "queue_name" in kwargs:
            queue_name = kwargs["queue_name"]
            cmd += f" -q={queue_name} "
        if "memory_mb" in kwargs:
            memory_mb = kwargs["memory_mb"]
            cmd += f" --memoryMb={memory_mb}"

        log_file_path = os.path.join(manager_data["logs_dir"], f"{label}.log")
        with open(log_file_path, "w") as log_file:
            self._process = subprocess.Popen(cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)

    def check_status(self):
        """Check if Para jobs are done."""
        running = self._process.poll() is None
        if not running:
            return self._process.returncode
        else:
            return None


class CustomStrategy(ParallelizationStrategy):
    """
    Custom parallel jobs execution strategy.
    """

    def execute(self, joblist_path, manager_data, label, **kwargs):
        """Custom implementation.

        Please provide your implementation of parallel jobs executor.
        Jobs are stored in the joblist_path, manager_data is a dict
        containing project-wide TOGA parameters."""
        raise NotImplementedError("Custom strategy is not implemented -> pls see documentation")

    def check_status(self):
        """Check if Para jobs are done.

        Please provide implementation of a method that checks
        whether all jobs are done."""
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

    def execute_jobs(self, joblist_path, manager_data, label, **kwargs):
        """
        Execute jobs in parallel using the specified strategy.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        """
        self.strategy.execute(joblist_path, manager_data, label, kwargs)

    def check_jobs_status(self):
        """
        Check the status of the jobs using the specified strategy.

        :return: Status of the jobs.
        """
        return self.strategy.check_status()
