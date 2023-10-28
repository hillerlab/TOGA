#!/usr/bin/env python3
"""Module to manage versioning."""

__author__ = "Bogdan M. Kirilenko"

class Version:
    def __init__(self, major, minor, patch, metadata=None):
        self.major = major
        self.minor = minor
        self.patch = patch
        self.metadata = metadata
        self.version_repr = f"{major}.{minor}.{patch}"
        if self.metadata:
            self.version_repr += f".{self.metadata}"

    def update_readme(self, filename="README.md"):
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            for line in lines:
                if "img.shields.io/badge/version-" in line:
                    line = f'![version](https://img.shields.io/badge/version-{self.version_repr}-blue)\n'
                f.write(line)

    def check_changelog(self, filename="VersionHistory.md"):
        with open(filename, "r") as f:
            header_lines = [x for x in f if x.startswith("# ") and x.endswith(" #\n")]
        header_lines_with_this_v = [x for x in header_lines if f" {self.version_repr} " in x]
        if len(header_lines_with_this_v) == 0:
            print(f"Warning! The version {self.version_repr} is absent in the {filename}")

    def __repr__(self):
        return self.version_repr

    def to_string(self):
        return self.version_repr


__version__ = Version(1, 1, 7, metadata="dev")

if __name__ == "__main__":
    print(f"TOGA version: {__version__}")
    __version__.update_readme()
    __version__.check_changelog()
