#!/usr/bin/env python3
"""Strategy pattern implementation to handle parallel jobs.

Provides implementations for nextflow and para strategies.
Please feel free to implement your custom strategy if
neither nextflow nor para satisfy your needs.
"""
from abc import ABC, abstractmethod
import subprocess


class ParallelizationStrategy(ABC):
    """
    Abstract base class for a parallelization strategy.
    """

    @abstractmethod
    def execute(self, joblist_path, manager_data, label):
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

    def execute(self, joblist_path, manager_data, label):
        """Implementation for Nextflow."""
        pass

    def check_status(self):
        """Check if nextflow jobs are done."""
        pass


class ParaStrategy(ParallelizationStrategy):
    """
    Concrete strategy for parallelization using Para.
    """

    def execute(self, joblist_path, manager_data, label):
        """Implementation for Para."""
        pass

    def check_status(self):
        """Check if Para jobs are done."""
        pass


class CustomStrategy(ParallelizationStrategy):
    """
    Custom parallel jobs execution strategy.
    """

    def execute(self, joblist_path, manager_data, label):
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

    def execute_jobs(self, joblist_path, manager_data, label):
        """
        Execute jobs in parallel using the specified strategy.

        :param joblist_path: Path to the joblist file.
        :param manager_data: Data from the manager class.
        :param label: Label for the run.
        """
        self.strategy.execute(joblist_path, manager_data, label)

    def check_jobs_status(self):
        """
        Check the status of the jobs using the specified strategy.

        :return: Status of the jobs.
        """
        return self.strategy.check_status()
