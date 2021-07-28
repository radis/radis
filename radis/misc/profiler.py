"""
`Profiler` class, for printing time spent on various steps during Spectrum
calculation under :py:class:`~radis.lbl.factory.SpectrumFactory` based on verbose level.

Also stores Spectrum calculation time dependent parameters, under the attribute
    :py:attr:`~radis.misc.profiler.Profiler.final`

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
        self.current_verbose2_counter = ""  # Holds current Verbose level 2 key
        self.final = {}
        self.flag = False  # For checking transit between verbose level 3 to 2
        self.relative_time_percentage = {}

    def start(self, key, verbose_level, optional=""):
        if __debug__:
            self.initial[key] = {
                "start_time": time(),
                "verbose_level": verbose_level,
            }

        if verbose_level == 1:
            self.final[key] = {"value": None}
        else:
            if verbose_level == 2:
                self.current_verbose2_counter = key
                self.final[list(self.final)[-1]].update({key: {}})
            else:
                self.final[list(self.final)[-1]][self.current_verbose2_counter].update(
                    {key: {}}
                )

        if len(optional) != 0 and self.verbose >= verbose_level:
            print(optional)

    def stop(self, key, details):
        if __debug__:
            items = self.initial.pop(key)
            time_calculated = time() - items["start_time"]
            if items["verbose_level"] == 1:
                self.final[key]["value"] = time_calculated
            else:

                if items["verbose_level"] == 2:
                    if self.flag:
                        self.final[list(self.final)[-1]][key].update(
                            {"value": time_calculated}
                        )
                    else:
                        self.final[list(self.final)[-1]][key] = time_calculated
                    self.flag = False
                else:
                    self.final[list(self.final)[-1]][self.current_verbose2_counter][
                        key
                    ] = time_calculated
                    self.flag = True

            if self.verbose >= items["verbose_level"]:
                self._print(
                    items["verbose_level"],
                    details,
                    key,
                    time_calculated=time_calculated,
                )

    def _print(self, verbose_level, details, key, time_calculated):

        if verbose_level == 1:
            print("{0:.2f}s -".format(time_calculated), details)
        elif verbose_level >= 2:
            printg(
                "..." * (verbose_level - 1),
                "{0:.2f}s -".format(time_calculated),
                details,
            )

    def percentage_distribution(self):
        """Computes and stores the percentage of time spent by each step of a particular verbose level.
        Stored in sorted order:
        {
            verbose_level: [
                ("key_1", [percentage_time, time]),
                ("key_2", [percentage_time, time]),
            ]
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
                for j in self.final:
                    if i in self.final[j]:
                        if first == 1:
                            verbose_distribution[i] = {}
                            first = 0
                        verbose_distribution[i][j] = [self.final[j][0]]
                        temp_sum += self.final[j][0]
                total_sum_verbose_wise[i] = temp_sum

            # Adding percentage of time taken by each key
            for i in range(1, upper_limit):
                for j in verbose_distribution[i]:
                    try:
                        verbose_distribution[i][j].insert(
                            0,
                            verbose_distribution[i][j][0]
                            / total_sum_verbose_wise[i]
                            * 100,
                        )
                    except ZeroDivisionError:
                        # Adding percentage of time = 0
                        verbose_distribution[i][j].insert(0, 0)

            # Sorting on the basis of value (Descending)
            for i in range(1, upper_limit):
                verbose_distribution[i] = sorted(
                    verbose_distribution[i].items(),
                    key=lambda kv: (kv[1], kv[0]),
                    reverse=True,
                )

            # Storing it in a parameter
            self.relative_time_percentage = verbose_distribution
