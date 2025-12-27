# -*- coding: utf-8 -*-
"""Created on Wed Oct 18 15:38:14 2017.

@author: erwan

-------------------------------------------------------------------------------
"""

import sys
from time import time


class ProgressBar:
    """Writes completion status and expended time.

    Set it up by creating the object, then calling
    :py:meth:`~radis.misc.progress_bar.ProgressBar.update` and
    :py:meth:`~radis.misc.progress_bar.ProgressBar.done`.

    Parameters
    ----------
    N: int
        (expected) number of iterations
    active: bool
        if ``False``, do not show anything (tip : feed it a ``verbose`` argument)

    Other Parameters
    ----------------
    t0: float
        initializes starting time at ``t0`` (useful for successive loops)

    Example
    -------
    add a progress bar in a loop::

        pb = ProgressBar(N)
        for i, ... in enumerate(...):
            (...)
            pb.update(i)
        pb.done()

    See test in progress_bar.py
    """

    # Todo: One day extend for multiprocss with several progress values?
    # https://stackoverflow.com/questions/7392779/is-it-possible-to-print-a-string-at-a-certain-screen-position-inside-idle

    def __init__(self, N, active=True, t0=None):

        self.t0 = time()
        if t0 is not None:
            self.t0 -= t0
        self.N = N
        self.active = active

    def set_active(self, active=True):
        """Option to activate/deactivate the ProgressBar.

        Used not to make it appear on small processes (based on a
        condition) without changing most of the code
        """
        self.active = active

    def update(self, i, modulo=1, message=""):
        """Update the completion status i/N and time spent.

        Parameters
        ----------
        i: int
            current iteration
        modulo: int
            if higher than ``1``, skip some iterations.
        message: str
            add a custom message. Tip: evaluate your own variables with
            f'{my_variable}' strings

        Example
        -------

        ::


        """
        if not self.active:
            return

        N = self.N
        t0 = self.t0
        if i % modulo == 0:
            # Show raw counts when total unknown (N=0 or None)
            if N is None or N == 0:
                if t0 is None:
                    msg = f"Parsed {i:,} lines\t{message}"
                else:
                    msg = f"({time() - t0:.0f}s)\tParsed {i:,} lines\t{message}"
            else:
                if t0 is None:
                    msg = f"{i / N * 100:.1f}%\t{message}"
                else:
                    msg = f"({time() - t0:.0f}s)\t{i / N * 100:.1f}%\t{message}"

            if sys.stdout is not None:
                sys.stdout.write("\r" + msg)
                sys.stdout.flush()

    def done(self):
        """Close the Progress bar."""
        if not self.active:
            return

        # Skip final update if total is unknown
        if self.N is not None and self.N != 0:
            self.update(self.N)
        # make new line
        if sys.stdout is not None:
            sys.stdout.write("\n")
            sys.stdout.flush()


# %% Helper functions for database fetch progress display


def print_database_header(
    database, molecule, wavenum_min=None, wavenum_max=None, isotope=None, **kwargs
):
    """Print a standardized header for database fetch operations.

    Parameters
    ----------
    database : str
        Database name (e.g., 'HITRAN', 'HITEMP 2024', 'ExoMol')
    molecule : str
        Molecule name
    wavenum_min, wavenum_max : float, optional
        Wavenumber range in cm-1
    isotope : str, optional
        Isotope identifier
    **kwargs : dict
        Additional info to display (e.g., local_folder, database_version)
    """
    print("-" * 80)
    header = f"{molecule} - {database}"
    if wavenum_min is not None and wavenum_max is not None:
        header += f" - {wavenum_min:.1f}-{wavenum_max:.1f} cm-1"
    print(header)
    print("-" * 80)

    if isotope is not None:
        print(f"Isotope: {isotope}")
    for key, value in kwargs.items():
        label = key.replace("_", " ").capitalize()
        print(f"{label}: {value}")


def print_database_section(section_name):
    """Print a section header with underline formatting.

    Parameters
    ----------
    section_name : str
        Section name (e.g., 'Download', 'Processing', 'Parsing')
    """
    print(f"\n\x1b[4m{section_name}:\x1b[0m")


# %% Now tested in radis/radis/test/test_misc.py
# def test_progress_bar(*args, **kwargs):
#    ''' Minimal example of a progress bar '''
#
#    from radis.misc.progress_bar import ProgressBar
#    from time import sleep
#    from numpy.random import rand
#
#    print('Testing progress bar')
#
#    a = 0
#    r = list(range(1000))
#    N = len(r)
#    pb = ProgressBar(N)
#    for i in r:
#        pb.update(i,modulo=10)
#        a += i
#        sleep(rand()*3e-3)
#    pb.done()
#
#    return True  # nothing implemented
#
#
# if __name__ == '__main__':
#    test_progress_bar()
