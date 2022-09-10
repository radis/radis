"""
`Profiler` class, for printing time spent on various steps during Spectrum
calculation under :py:class:`~radis.lbl.factory.SpectrumFactory` based on verbose level.

Also stores Spectrum calculation time dependent parameters, under the attribute
    :py:attr:`~radis.misc.profiler.Profiler.final`

Routine Listing
---------------

- :meth:`~radis.misc.profiler.Profiler.add_entry`
- :meth:`~radis.misc.profiler.Profiler.start`
- :meth:`~radis.misc.profiler.Profiler.add_time`
- :meth:`~radis.misc.profiler.Profiler.stop`
- :meth:`~radis.misc.profiler.Profiler._print`

-------------------------------------------------------------------------------
"""

from collections import OrderedDict
from time import perf_counter

from radis.misc.printer import printg


class Profiler(object):
    """A class to store Spectrum calculation time dependent parameters, under the attribute
    :py:attr:`~radis.misc.profiler.Profiler.final` of
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
        self.verbose = verbose
        self.final = OrderedDict()
        self.relative_time_percentage = {}

    # Creates profiler dictionary structure
    def add_entry(self, dictionary, key, verbose, count):
        if count == verbose:
            if key in dictionary:
                return
            else:
                dictionary[key] = {}

            return

        self.add_entry(dictionary[list(dictionary)[-1]], key, verbose, count + 1)

    def start(self, key, verbose_level, optional=""):
        if __debug__:
            self.initial[key] = {
                "start_time": perf_counter(),
                "verbose_level": verbose_level,
            }

            self.add_entry(self.final, key, verbose_level, 1)

        if len(optional) != 0 and self.verbose >= verbose_level:
            print(optional)

    # Adds time calculated for each key in profiler
    def add_time(self, dictionary, key, verbose, count, time_calculated):
        if count == verbose:
            executed_more_than_once = 0  # Checks if input key was executed before

            # Checking if value of key is empty or float
            if isinstance(dictionary[key], float):
                time_calculated = dictionary[key] + time_calculated
                executed_more_than_once = 1
            else:
                try:
                    if "value" in dictionary[key] and isinstance(
                        dictionary[key]["value"], float
                    ):
                        time_calculated = dictionary[key]["value"] + time_calculated
                        executed_more_than_once = 2
                except KeyError:
                    pass

            # Adding values
            if executed_more_than_once == 2:
                dictionary[key]["value"] = time_calculated
            elif executed_more_than_once == 1:
                dictionary[key] = time_calculated
            else:
                if len(dictionary[key]) != 0:
                    dictionary[key].update({"value": time_calculated})
                else:
                    dictionary[key] = time_calculated
            return

        self.add_time(
            dictionary[list(dictionary)[-1]], key, verbose, count + 1, time_calculated
        )

    def stop(self, key, details):
        if __debug__:
            items = self.initial.pop(key)
            time_calculated = perf_counter() - items["start_time"]

            if items["verbose_level"] == 1:
                # Profiler ends; Deserializing to Dictionary format
                self.final = dict(self.final)

            self.add_time(self.final, key, items["verbose_level"], 1, time_calculated)

            if self.verbose >= items["verbose_level"]:
                self._print(
                    items["verbose_level"],
                    details,
                    time_calculated=time_calculated,
                )

    def _print(self, verbose_level, details, time_calculated):

        if verbose_level == 1:
            print("{0:.2f}s -".format(time_calculated), details)
        elif verbose_level >= 2:
            printg(
                "..." * (verbose_level - 1),
                "{0:.2f}s -".format(time_calculated),
                details,
            )
