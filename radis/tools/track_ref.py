# -*- coding: utf-8 -*-
"""
Experiment module to track references used in the code; and return
the appropriate citations

To be made an independant package in the future
"""


class RefTracker(dict):
    """Track bibliographic references by storing the DOIs in a dictionary.

    Prints them nicely using the :py:mod:`habanero` CrossRef API package.

    Examples
    --------
    ::

        rt = RefTracker()

        # in your code:
        rt.add("10.1016/a.test.2021.012345", "algorithm")
        rt.add("10.1016/a.test.2021.012345", "post-processing")
        rt.add("10.1016/a.test.2021.678910", "data retrieval")
        rt.add("10.1016/a.test.2021.111213", "data retrieval")

        # user part:
        rt.cite()

    Returns::

        -----------------------------------
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.012345,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.012345 in bibentry format }

        Used for data retrieval
        -----------------------
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.678910,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.678910 in bibentry format }
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.111213,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.111213 in bibentry format }

    Other Example
    -------------

    Can be initialized from an existing dict (useful when serializing with JSON) ::

        ref_dict = {'10.1016/a.test.2021.012345': ['algorithm', 'post-processing'],
             '10.1016/a.test.2021.678910': ['data retrieval'],
             '10.1016/a.test.2021.111213': ['data retrieval']}

        rt = RefTracker(**ref_dict)
        rt.cite()

    Returns::

        -----------------------------------
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.012345,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.012345 in bibentry format }

        Used for data retrieval
        -----------------------
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.678910,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.678910 in bibentry format }
        @article{DUMMY_KEY,
                doi = 10.1016/a.test.2021.111213,
                title = PLACEHOLDER FOR 10.1016/a.test.2021.111213 in bibentry format }


    .. minigallery:: radis.tools.track_ref.RefTracker

    See Also
    --------
    :py:meth:`~radis.spectrum.spectrum.Spectrum.cite`


    """

    def __init__(self, **kwargs):
        super(RefTracker).__init__()

        if kwargs:
            for k, v in kwargs.items():
                if not isinstance(v, list):
                    v = [v]
                for vi in v:
                    self.add(k, vi)

    def add(self, doi, why=""):
        if not doi in self:
            self[doi] = [why]
        else:
            if why not in self[doi]:  # no duplicates
                self[doi].append(why)

    def cite(self, format="bibentry"):
        """See more in :py:func:`habanero.content_negotiation`"""

        # Grouping dictionary by value :
        by_use = {}
        for doi, why in self.items():
            why_str = ", ".join(why)
            by_use[why_str] = (
                [doi] if why_str not in by_use.keys() else by_use[why_str] + [doi]
            )

        for why, dois in sorted(by_use.items()):
            msg = f"Used for {why}"
            print(msg)
            print("-" * len(msg))
            for doi in dois:
                print(self._content(doi, format=format))
            print("")

    def _content(self, doi, format):
        if doi.startswith("10.1016/a.test"):
            return f"@article{{DUMMY_KEY,\n\tdoi = {doi}, \n\ttitle = PLACEHOLDER FOR {doi} in {format} format }}"

        try:
            from habanero import cn
        except ImportError:
            raise ImportError(
                "Printing citations require the `habanero` CrossRef API package. Install with `pip install habanero`"
            )
        return cn.content_negotiation(ids=doi, format=format)


if __name__ == "__main__":

    import pytest

    print("Testing track_ref:", pytest.main(["../test/tools/test_track_ref.py"]))
