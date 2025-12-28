# -*- coding: utf-8 -*-
"""
Database Progress Printer Utility

Provides unified progress output for database fetch operations across RADIS.
Addresses GitHub Issue #868: Same printing when downloading from different databases.

Created on 2025-12-28

@author: RADIS contributors
"""

import sys
from typing import Optional

from tqdm import tqdm


class DatabaseProgressPrinter:
    """Unified progress output for database fetch operations.

    Provides consistent formatting across all database operations (HITRAN, HITEMP,
    ExoMol, GEISA) with configurable verbose levels.

    Parameters
    ----------
    database_name : str
        Name of the database (e.g., "HITEMP", "ExoMol", "HITRAN", "GEISA")
    molecule : str
        Molecule being fetched (e.g., "CO2", "H2O")
    verbose : int or bool
        Verbosity level:
        - 0 or False: Silent (no output)
        - 1 or True: Normal output (default)
        - 2+: Detailed/debug output

    Example
    -------
    ::

        printer = DatabaseProgressPrinter("HITEMP", "CO2", verbose=1)
        printer.header("Fetching 2024 database")
        printer.section("Downloading")
        printer.info("3 files needed, 2 cached")
        with printer.progress_bar(3, "Downloading chunks") as pbar:
            for i in range(3):
                # download work
                pbar.update(1)
        printer.success("Complete: Added HITEMP-CO2-2024 to radis.json")
    """

    # Header width for consistent formatting
    HEADER_WIDTH = 80

    def __init__(
        self, database_name: str, molecule: str, verbose: int = 1, version: str = ""
    ):
        self.database_name = database_name.upper()
        self.molecule = molecule
        self.verbose = int(verbose) if isinstance(verbose, (int, bool)) else 1
        self.version = version

    def _print(self, message: str, min_verbose: int = 1):
        """Internal print method that respects verbose level."""
        if self.verbose >= min_verbose and sys.stdout is not None:
            print(message)
            sys.stdout.flush()

    def header(self, action: str = "", details: str = ""):
        """Print operation header with visual separator.

        Parameters
        ----------
        action : str
            Action being performed (e.g., "Fetching", "Downloading")
        details : str, optional
            Additional details for the header
        """
        if self.verbose < 1:
            return

        separator = "-" * self.HEADER_WIDTH
        version_str = f" {self.version}" if self.version else ""
        title = f"{self.molecule} - {self.database_name}{version_str}"
        if action:
            title += f" - {action}"
        if details:
            title += f" {details}"

        self._print(separator)
        self._print(title)
        self._print(separator)

    def section(self, title: str, level: int = 1):
        """Start a new section (Download, Parse, etc.).

        Parameters
        ----------
        title : str
            Section title (e.g., "Download", "Parsing", "Processing")
        level : int
            Minimum verbose level to show this section
        """
        if self.verbose >= level:
            self._print(f"\n{title}:")

    def info(self, message: str, level: int = 1, indent: int = 0):
        """Print info message if verbose >= level.

        Parameters
        ----------
        message : str
            Message to display
        level : int
            Minimum verbose level required to show message
        indent : int
            Indentation level (0=none, 1=bullet, 2=sub-bullet)
        """
        if self.verbose >= level:
            prefix = ""
            if indent == 1:
                prefix = "- "
            elif indent >= 2:
                prefix = "  " * (indent - 1) + "- "
            self._print(f"{prefix}{message}")

    def success(self, message: str):
        """Print success message.

        Parameters
        ----------
        message : str
            Success message to display
        """
        self._print(f"\n{message}")

    def warning(self, message: str):
        """Print warning message.

        Parameters
        ----------
        message : str
            Warning message to display
        """
        self._print(f"Warning: {message}")

    def progress_bar(
        self,
        total: int,
        desc: str = "",
        unit: str = "file",
        unit_scale: bool = False,
        unit_divisor: int = 1024,
        leave: bool = True,
        position: int = 0,
        disable: Optional[bool] = None,
    ) -> tqdm:
        """Return a configured tqdm progress bar.

        Parameters
        ----------
        total : int
            Total number of iterations
        desc : str
            Description text for the progress bar
        unit : str
            Unit name (e.g., "file", "line", "B" for bytes)
        unit_scale : bool
            If True, scale units automatically (useful for bytes)
        unit_divisor : int
            Divisor for unit scaling (default 1024 for bytes)
        leave : bool
            If True, keep progress bar after completion
        position : int
            Position for nested progress bars
        disable : bool or None
            If True, disable progress bar. If None, use verbose level.

        Returns
        -------
        tqdm
            Configured progress bar instance
        """
        if disable is None:
            disable = self.verbose < 1

        return tqdm(
            total=total,
            desc=desc,
            unit=unit,
            unit_scale=unit_scale,
            unit_divisor=unit_divisor,
            leave=leave,
            position=position,
            disable=disable,
        )

    def download_progress(
        self,
        total_bytes: int,
        desc: str = "Downloading",
    ) -> tqdm:
        """Return a progress bar configured for download operations.

        Parameters
        ----------
        total_bytes : int
            Total download size in bytes
        desc : str
            Description text

        Returns
        -------
        tqdm
            Progress bar configured for byte downloads
        """
        return self.progress_bar(
            total=total_bytes,
            desc=desc,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        )

    def parsing_progress(
        self,
        total_files: int,
        desc: str = "Processing files",
    ) -> tqdm:
        """Return a progress bar configured for file parsing operations.

        Parameters
        ----------
        total_files : int
            Total number of files to process
        desc : str
            Description text

        Returns
        -------
        tqdm
            Progress bar configured for file processing
        """
        return self.progress_bar(
            total=total_files,
            desc=desc,
            unit="file",
        )

    def lines_progress(
        self,
        total_lines: Optional[int],
        desc: str = "Parsing lines",
    ) -> tqdm:
        """Return a progress bar configured for line parsing operations.

        Parameters
        ----------
        total_lines : int or None
            Total number of lines expected (None if unknown)
        desc : str
            Description text

        Returns
        -------
        tqdm
            Progress bar configured for line parsing
        """
        return self.progress_bar(
            total=total_lines if total_lines else 0,
            desc=desc,
            unit=" lines",
            unit_scale=True,
        )

    def download_summary(
        self, files_needed: int, files_total: int, files_cached: int = 0
    ):
        """Print download summary.

        Parameters
        ----------
        files_needed : int
            Number of files that need to be downloaded
        files_total : int
            Total number of files
        files_cached : int
            Number of files already cached
        """
        if files_needed == 0:
            self.info("All files already cached.", indent=1)
        else:
            msg = f"Download {files_needed} file(s) missing out of {files_total}."
            if files_cached > 0:
                msg += f" ({files_cached} cached)"
            self.info(msg, indent=1)

    def parsing_summary(self, lines_loaded: int, time_elapsed: Optional[float] = None):
        """Print parsing summary.

        Parameters
        ----------
        lines_loaded : int
            Number of lines loaded
        time_elapsed : float, optional
            Time taken in seconds
        """
        msg = f"Lines loaded: {lines_loaded:,}"
        if time_elapsed is not None:
            msg += f" ({time_elapsed:.1f}s)"
        self.info(msg, indent=1)

    def complete(self, database_entry: str = ""):
        """Print completion message.

        Parameters
        ----------
        database_entry : str
            Name of database entry added to radis.json
        """
        if database_entry:
            self.success(f"Added {database_entry} database in radis.json")
        else:
            self.success("Complete")


# Convenience function for quick usage
def get_progress_printer(
    database_name: str, molecule: str, verbose: int = 1, version: str = ""
) -> DatabaseProgressPrinter:
    """Create a DatabaseProgressPrinter instance.

    Parameters
    ----------
    database_name : str
        Name of the database
    molecule : str
        Molecule being fetched
    verbose : int
        Verbosity level
    version : str
        Database version

    Returns
    -------
    DatabaseProgressPrinter
        Configured printer instance
    """
    return DatabaseProgressPrinter(database_name, molecule, verbose, version)
