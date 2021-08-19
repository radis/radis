"""
Run all examples in `radis/examples <https://github.com/radis/radis/tree/develop/examples>`__
"""

import os
import runpy
from os.path import abspath, dirname, join

import pytest

EXAMPLE_FOLDER = abspath(join(dirname(__file__), "..", "..", "examples"))

scripts = [
    join(EXAMPLE_FOLDER, k) for k in os.listdir(EXAMPLE_FOLDER) if k.endswith(".py")
]


@pytest.mark.parametrize("script", scripts)
def test_script_execution(script):
    """Run all examples in `radis/examples <https://github.com/radis/radis/tree/develop/examples>`__"""
    # ensure Matplotlib is interactive (figures won't block the rest of the execution)
    import matplotlib.pyplot as plt

    plt.ion()
    # run:
    runpy.run_path(script, init_globals=locals())


if __name__ == "__main__":
    print(scripts)
    pytest.main(["test_examples.py", "--pdb"])
