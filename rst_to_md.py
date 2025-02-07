"""RST to MD conversion.

- auto-convert README.rst to long_description, removing some sphinx-only syntax
so it can be rendered by PyPi
- read version number from __version__.txt

Notes
-----

For developers:

when creating a new version, just update the __version__.txt file

to register it on Pypi see register.py::

    python register.py
"""
import io
import re
import sys
from os.path import abspath, dirname, exists, join

# from setuptools import Extension, find_packages, setup #for cython - unused after v0.15
# from setuptools import setup

# Build description from README (PyPi compatible)
# -----------------------------------------------

# Utils to format RST
def yield_sphinx_only_markup(lines):
    """Cleans-up Sphinx-only constructs (ie from README.rst), so that *PyPi*
    can format it properly.

    To check for remaining errors, install ``sphinx`` and run::

        python setup.py --long-description | sed -file 'this_file.sed' | rst2html.py  --halt=warning

    :param file_inp:     a `filename` or ``sys.stdin``?
    :param file_out:     a `filename` or ``sys.stdout`?`

    References
    ----------

    https://stackoverflow.com/questions/16367770/my-rst-readme-is-not-formatted-on-pypi-python-org

    Notes
    -----

    Check output with::

        python setup.py --long-description | rst2html.py > output.html
    """
    substs = [
        ## Selected Sphinx-only Roles.
        #
        (r":abbr:`([^`]+)`", r"\1"),
        (r":ref:`([^`]+)`", r"`\1`_"),
        (r":term:`([^`]+)`", r"**\1**"),
        (r":dfn:`([^`]+)`", r"**\1**"),
        (r":(samp|guilabel|menuselection):`([^`]+)`", r"``\2``"),
        ## Sphinx-only roles:
        #        :foo:`bar`   --> foo(``bar``)
        #        :a:foo:`bar` XXX afoo(``bar``)
        #
        # (r'(:(\w+))?:(\w+):`([^`]*)`', r'\2\3(``\4``)'),
        (r":(\w+):`([^`]*)`", r"\1(``\2``)"),
        ## Sphinx-only Directives.
        #
        (r"\.\. doctest", r"code-block"),
        (r"\.\. plot::", r".. "),
        (r"\.\. seealso", r"info"),
        (r"\.\. glossary", r"rubric"),
        (r"\.\. figure::", r".. "),
        ## Other
        #
        (r"\|version\|", r"x.x.x"),
        ## added to make RADIS docs Pypi compatible
        #        (r'\.\. image::',          r'.. '),
        #        (r'\.\. |CO2| replace:: CO\ :sub:`2`',          r'.. '),
    ]

    regex_subs = [(re.compile(regex, re.IGNORECASE), sub) for (regex, sub) in substs]

    def clean_line(line):
        try:
            for (regex, sub) in regex_subs:
                line = regex.sub(sub, line)
        except Exception as ex:
            print(("ERROR: %s, (line(%s)" % (regex, sub)))
            raise ex

        return line

    for line in lines:
        yield clean_line(line)


# (note: README.rst has been converted to README.md by register.py, and cleaned afterwards )
description = "A fast line-by-line code for high-resolution infrared molecular spectra"
readme_path = join(abspath(dirname(__file__)), "README.md")
if not exists(readme_path):
    long_description = description
else:
    with io.open(readme_path, encoding="utf-8") as f:
        long_description = f.read()

# Read version number from file
with open(join(dirname(__file__), "radis", "__version__.txt")) as version_file:
    __version__ = version_file.read().strip()


# %% Cython extensions:
# Ensures RADIS still builds if without Cython


class BuildFailed(Exception):
    pass


# from distutils.errors import CCompilerError, DistutilsExecError, DistutilsPlatformError #for cython - unused after v0.15


def show_message(*lines):
    """Note : will only happen if user installs with `pip install -v`"""
    print("=" * 74, file=sys.stderr)
    for line in lines:
        print(line, file=sys.stderr)
    print("=" * 74, file=sys.stderr)
