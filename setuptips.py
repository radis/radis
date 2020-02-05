# -*- coding: utf-8 -*-
"""
Utils for setup.py
"""

from __future__ import absolute_import
from __future__ import print_function
import re

# Utils to format RST
def yield_sphinx_only_markup(lines):
    """
    Cleans-up Sphinx-only constructs (ie from README.rst),
    so that *PyPi* can format it properly.

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


if __name__ == "__main__":
    import codecs

    readme_lines = codecs.open("README.rst", encoding="utf-8").readlines()
    long_description = "".join(yield_sphinx_only_markup(readme_lines))
    print(long_description)
