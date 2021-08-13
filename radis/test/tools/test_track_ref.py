# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:27:22 2021

@author: erwan
"""

from radis.db.references import doi
from radis.tools.track_ref import RefTracker


def test_reftracker(verbose=True, *args, **kwargs):

    rt = RefTracker()

    # in your code:
    rt.add("10.1016/a.test.2021.012345", "algorithm")
    rt.add("10.1016/a.test.2021.012345", "post-processing")
    rt.add("10.1016/a.test.2021.678910", "data retrieval")
    rt.add("10.1016/a.test.2021.111213", "data retrieval")

    # user part:
    if verbose:
        rt.cite()

    # Init from a dict directly
    ref_dict = dict(rt)
    rt2 = RefTracker(**ref_dict)

    assert dict(rt) == dict(rt2)


def test_citations_in_eq_spectrum(verbose=True, *args, **kwargs):

    from radis import calc_spectrum

    s = calc_spectrum(
        1900,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        mole_fraction=0.1,
        path_length=1,  # cm
        databank="hitran",  # or use 'hitemp'
    )

    if verbose:
        s.cite()

    assert doi["RADIS-2018"] in s.references
    assert doi["DIT-2020"] in s.references
    assert doi["HITRAN-2016"] in s.references
    assert doi["TIPS-2020"] in s.references


if __name__ == "__main__":
    verbose = True

    # test_reftracker()
    # test_citations_in_eq_spectrum()

    from radis import calc_spectrum

    s = calc_spectrum(
        1900,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tvib=2000,  #
        Trot=300,
        mole_fraction=0.1,
        path_length=1,  # cm
        databank="hitran",  # or use 'hitemp'
    )

    # if verbose:
    #     s.cite()

    assert doi["RADIS-2018"] in s.references
    assert doi["DIT-2020"] in s.references
    assert doi["HITRAN-2016"] in s.references
    assert doi["TIPS-2020"] in s.references

    assert doi["Guelachvili-1983"] in s.references
