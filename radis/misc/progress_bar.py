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
            if t0 is None:
                msg = "{0:.1f}%\t{1}".format(i / N * 100, message)
            else:
                msg = "({0:.0f}s)\t{1:.1f}%\t{2}".format(
                    time() - t0, i / N * 100, message
                )

            if sys.stdout is not None:
                sys.stdout.write("\r" + msg)
                sys.stdout.flush()

    def done(self):
        """Close the Progress bar."""
        if not self.active:
            return

        self.update(self.N)
        # make new line
        if sys.stdout is not None:
            sys.stdout.write("\n")
            sys.stdout.flush()


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
