"""
`Profiler` class, for printing time spent on various steps during Spectrum
calculation under :py:class:`~radis.lbl.factory.SpectrumFactory` based on verbose level.

Also stores Spectrum calculation time dependent parameters, under the attribute
    :py:attr:`~radis.misc.profiler.Profiler.dict_time`

Routine Listing
---------------

- :meth:`~radis.misc.profiler.Profiler.start`
- :meth:`~radis.misc.profiler.Profiler.stop`
- :meth:`~radis.misc.profiler.Profiler._print`

-------------------------------------------------------------------------------
"""

from radis.misc.printer import printg
from time import time

class Profiler(object):
    """A class to store Spectrum calculation time dependent parameters, under the attribute
    :py:attr:`~radis.misc.profiler.Profiler.dict_time` of
    :py:class:`~radis.lbl.factory.SpectrumFactory`.

    It also hold functions to print all the entities based on verbose value.

    See Also
    --------

    :py:attr:`~radis.lbl.loader.DatabankLoader.input`,
    :py:attr:`~radis.misc.profiler.Profiler,

    """

    def __init__(self, verbose):
        super(Profiler, self).__init__()

        # Dev: Init here to be found by autocomplete
        self.initial = {}
        self.dict_time = {}
        self.verbose = verbose

    def start(self, key, verbose, details=""):
        if __debug__:
            self.initial[key] = {
                "start_time": time(),
                "verbose": verbose,
                "details": details,
            }

    def stop(self, key):
        if __debug__:
            items = self.initial.pop(key)
            self.dict_time[key] = time() - items["start_time"]
            if self.verbose >= items["verbose"]:
                self._print(key, items["verbose"], items["details"])

    def _print(self, key, verbose, details):
        if verbose == 1:
            if key == None:
                print(details)
            else:
                print(details, "in {0:.1f}s".format(self.dict_time[key]))
        elif verbose >= 2:
            if key == None:
                printg("..." * (verbose - 1), details)
            else:
                printg(
                    "..." * (verbose - 1),
                    details,
                    "in {0:.2f}s".format(self.dict_time[key]),
                )