# -*- coding: utf-8 -*-
"""
Summary
-------

Runs tests for neq/misc so that they can be accessed by pytest (and hopefully
the CI test suite)

Examples
--------

Run all tests::

    pytest      # (in command line, in project folder)

Run only fast tests (i.e: tests that a  'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------

"""

import pytest


@pytest.mark.fast
def test_utils(*args, **kwargs):
    from radis.misc.utils import _test

    assert _test(*args, **kwargs)


@pytest.mark.fast
def test_progress_bar(*args, **kwargs):
    """Minimal example of a progress bar."""

    from time import sleep

    from numpy.random import rand

    from radis.misc.progress_bar import ProgressBar

    print("Testing progress bar")

    a = 0
    r = list(range(200))
    N = len(r)
    pb = ProgressBar(N)
    for i in r:
        pb.update(i, modulo=10)
        a += i
        sleep(rand() * 3e-3)
    pb.done()


@pytest.mark.fast
def test_norm():
    import numpy as np

    from radis.misc import norm

    a = np.random.rand(100)
    b = np.random.rand(100)

    c = norm(a)
    d = norm(c, normby=b)

    result = norm(d, how="mean")

    assert isinstance(result, np.ndarray)


@pytest.mark.fast
def test_normOn():

    import numpy as np

    from radis.misc import norm_on

    I = np.array([10, 9, 8, 7, 6, 5, 0, 9])
    w = np.arange(8)

    expected = np.array([1.25, 1.125, 1.0, 0.875, 0.75, 0.625, 0.0, 1.125])
    result = norm_on(I, w, wmin=2, wmax=5)

    assert np.allclose(result, expected)


@pytest.mark.fast
def test_allclose():

    import numpy as np

    from radis.misc import array_allclose

    a = np.random.rand(20)
    b = a + 0.1
    c = b
    c[0] = float("NaN")
    assert not array_allclose(a, b)
    assert not array_allclose(a, c)


@pytest.mark.fast
def test_nantrapz():

    import numpy as np

    from radis.misc import nantrapz

    I = np.array([10, 9, 8, 7, float("NaN"), 5, 0, 9])
    w = np.arange(8)

    result = nantrapz(I, w, dx=1.0, axis=-1)
    assert np.isclose(result, 44.5, rtol=1e-5)


def _run_testcases(verbose=True, *args, **kwargs):

    test_utils(verbose=verbose, *args, **kwargs)
    test_progress_bar()
    test_nantrapz()
    test_allclose()
    test_normOn()
    test_norm()


if __name__ == "__main__":
    print(("Testing misc.py: ", _run_testcases(verbose=True)))
