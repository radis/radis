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

from time import time

from radis.misc.printer import printg


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
        self.relative_time_percentage = {}

    def start(self, key, verbose):
        if __debug__:
            self.initial[key] = {
                "start_time": time(),
                "verbose": verbose,
            }

    def stop(self, key, details=""):
        if __debug__:
            items = self.initial.pop(key)
            # Storing time with verbose level
            self.dict_time[key] = [time() - items["start_time"], items["verbose"]]
            if self.verbose >= items["verbose"]:
                self._print(items["verbose"], details, key)

    def _print(self, verbose, details, key=None):
        if verbose == 1:
            if key == None:
                print(details)
            else:
                print("{0:.2f}s -".format(self.dict_time[key][0]), details)
        elif verbose >= 2:
            if key == None:
                printg("..." * (verbose - 1), details)
            else:
                printg(
                    "..." * (verbose - 1),
                    "{0:.2f}s -".format(self.dict_time[key][0]),
                    details,
                )

    def percentage_distribution(self):
        """Computes and stores the percentage of time spent by each step of a particular verbose level.
        Stored in sorted order:
        {
            verbose_level: {
                "key_1": [percentage_time, time],
                "key_2": [percentage_time, time],
            }
        }

        """
        # {verbose_level: sum of time}
        total_sum_verbose_wise = {}
        # {verbose_level: {key: value}}
        verbose_distribution = {}

        if self.verbose:
            upper_limit = self.verbose + 1

            for i in range(1, upper_limit):
                temp_sum = 0
                first = 1
                for j in self.dict_time:
                    if i in self.dict_time[j]:
                        if first == 1:
                            verbose_distribution[i] = {}
                            first = 0
                        verbose_distribution[i][j] = [self.dict_time[j][0]]
                        temp_sum += self.dict_time[j][0]
                total_sum_verbose_wise[i] = temp_sum

            # Adding percentage of time taken by each key
            for i in range(1, upper_limit):
                for j in verbose_distribution[i]:
                    verbose_distribution[i][j].insert(
                        0,
                        verbose_distribution[i][j][0] / total_sum_verbose_wise[i] * 100,
                    )

            # Sorting on the basis of value (Descending)
            for i in range(1, upper_limit):
                verbose_distribution[i] = sorted(
                    verbose_distribution[i].items(),
                    key=lambda kv: (kv[1], kv[0]),
                    reverse=True,
                )

            # Storing it in a parameter
            self.relative_time_percentage = verbose_distribution
