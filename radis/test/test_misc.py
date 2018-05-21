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

from __future__ import print_function, absolute_import, division, unicode_literals
import pytest


@pytest.mark.fast
def test_config(*args, **kwargs):
    from radis.misc.config import _test
    return _test(*args, **kwargs)


@pytest.mark.fast
def test_utils(*args, **kwargs):
    from radis.misc.utils import _test
    return _test(*args, **kwargs)


def _run_testcases(verbose=True, *args, **kwargs):

    test_config(verbose=verbose, *args, **kwargs)
    test_utils(verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    print(('Testing misc.py: ', _run_testcases(verbose=True)))
