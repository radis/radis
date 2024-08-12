"""Install file for RADIS.

Typical install procedure, plus:

- auto-convert README.rst to long_description, removing some sphinx-only syntax
so it can be rendered by PyPi
- read version number from __version__.txt
- some dependencies should be installed manually (typically: all packages with
compiled components such as numpy, pandas, etc.)


Examples
--------

Install (normal, use-only)::

    python setup.py install

Or (create an alias, so you can still edit)::

    python setup.py develop

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
from setuptools import setup

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


# Deprecated - See Main install routine for radis<0.15
# def get_ext_modules(with_binaries):
#     """
#     Parameters
#     ----------
#     with_binaries: bool
#         if False, do not try to build Cython extensions
#     """
#     # Based on code from https://github.com/pallets/markupsafe/blob/main/setup.py

#     print(sys.version)
#     ext_modules = []
#     cmdclass = {}

#     if not with_binaries:
#         # skip building extensions
#         return {"cmdclass": cmdclass, "ext_modules": ext_modules}

#     # TO-DO: set language level

#     #### Cython is now completely removed - August 2024 - Radis 0.15 ###
#     # try:
#     #     import cython
#     # except (ModuleNotFoundError) as err:
#     #     raise BuildFailed(
#     #         "Cython not found : Skipping all Cython extensions...!"
#     #     ) from err

#     # print("Cython " + cython.__version__)

#     # from Cython.Distutils import build_ext

#     # class build_ext_subclass(build_ext):
#     #     def build_extensions(self):
#     #         c = self.compiler.compiler_type
#     #         copt = {
#     #             "msvc": ["/openmp", "/Ox", "/fp:fast", "/favor:INTEL64"],
#     #             "mingw32": ["-fopenmp", "-O3", "-ffast-math", "-march=native"],
#     #         }

#     #         lopt = {"mingw32": ["-fopenmp"]}

#     #         print("Compiling with " + c + "...")
#     #         try:
#     #             for e in self.extensions:
#     #                 e.extra_compile_args = copt[c]
#     #         except (KeyError):
#     #             pass
#     #         try:
#     #             for e in self.extensions:
#     #                 e.extra_link_args = lopt[c]
#     #         except (KeyError):
#     #             pass
#     #         try:
#     #             build_ext.build_extensions(self)
#     #         except (CCompilerError, DistutilsExecError, DistutilsPlatformError) as err:
#     #             raise BuildFailed() from err

#     # from numpy import get_include

#     # ext_modules.append(
#     #     Extension(
#     #         "radis_cython_extensions",
#     #         sources=[
#     #             "./radis/cython/radis_cython_extensions.pyx",
#     #         ],
#     #         include_dirs=[get_include()],
#     #         language="c++",
#     #         extra_link_args=[],
#     #         define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
#     #     )
#     # )

#     # cmdclass["build_ext"] = build_ext_subclass

#     return {"cmdclass": cmdclass, "ext_modules": ext_modules}


#%% Main install routine
def run_setup():
    setup(
        name="radis",
        version=__version__,
        description=description,
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/radis/radis",
        author="Erwan Pannier",
        author_email="erwan.pannier@gmail.com",
        license="GNU Lesser General Public License v3 (LGPLv3)",
        keywords=[
            "spectrum",
            "infrared",
            "spectra",
            "radiation",
            "nonequilibrium",
            "spectroscopy",
            "molecules",
            "HITRAN",
            "hitemp",
            "exomol",
            "line-by-line",
        ],
        # packages=find_packages(), #this misses the gpu folder, see https://github.com/radis/radis/issues/681
        install_requires=[
            # Although it seems like a duplicate of requirements.txt, it is NECESSARY to tell pip what are the dependencies, issue #680.
            # (requirements.txt is only for the user, not for pip, see https://packaging.python.org/discussions/install-requires-vs-requirements/ and https://stackoverflow.com/questions/14399534/reference-requirements-txt-for-the-install-requires-kwarg-in-setuptools-setup-py)
            ### From environment.yml
            "astropy>=4.3.1",  # Unit aware calculations
            "astroquery>=0.4.6",  # to fetch HITRAN databases
            "beautifulsoup4>=4.10.0",  # parse ExoMol website
            "cantera>=2.5.1",  # for chemical equilibrium computations
            "configparser",
            "habanero>=1.2.0",  # CrossRef API to retrieve data from doi
            "h5py>=3.2.1",  # load HDF5
            "joblib",  # for parallel loading of SpecDatabase
            "lmfit",  # for new fitting modules
            "matplotlib",
            "numpy",
            "numba",  # just-in-time compiler
            "pandas",
            "plotly>=2.5.1",  # for line survey HTML output
            "psutil",  # to get user RAM
            "tables",  # for pandas to HDF5 export - WARNING named "pytables" in conda
            "scipy>=1.4.0",
            "seaborn",  # other matplotlib themes
            "termcolor",  # terminal colors
            "specutils",
            ### From requirements.txt
            "lxml",  # parser used for ExoMol website
            "hjson",  # Json with comments (for default_radis.json)
            "publib",  # Plotting styles for Matplotlib. Version, update to 0.4.0 needed for matplotlib==3.8. However, only 0.3.2 is compatible with python==3.8 (see #647)
            "hitran-api",  # HAPI, used to access TIPS partition functions
            "peakutils",
            "ruamel.yaml",
            "json-tricks>=3.15.0",  # to deal with non jsonable formats
            "mpldatacursor",
            "nvidia-cufft-cu11",
            "periodictable",
            'vaex-core ; python_version < "3.11"',
            'vaex-hdf5 ; python_version < "3.11"',
            'vaex-viz ; python_version < "3.11"',
        ],
        extras_require={
            "dev": [
                "numpydoc",  # for Jedi (autocompletion) to recognize
                "black>=20.8b1",  # for code-linting in accordance to PEP8
                "isort",  # for sorting imports
                "pre-commit",  # to enforce Black before each commit
                "pytest",  # to run test suite
                "ipython>=7.0.0",  # useful for fast debugging
            ],
            "docs": [
                "sphinx-autodoc-annotation",
                "sphinx_autodoc_defaultargs>=0.1.2",
                "sphinx>=1.7.0",
                "astroquery>=0.3.9",
                "sphinxcontrib-apidoc",
                "sphinx-gallery",
                "lmfit",
                "pytest",  # Sphinx autodoc also parses test suite
            ],
        },
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
            "Topic :: Scientific/Engineering",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.8",  # end of life, Novembre 2024, https://devguide.python.org/versions/
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Operating System :: OS Independent",
        ],
        # **get_ext_modules(with_binary), #see Main install routine for radis<0.15
        include_package_data=True,  # add non .py data files in MANIFEST.in
        # package_data={'radis': ['radis/phys/units.txt']},
        zip_safe=False,  # impossible as long as we have external files read with __file__ syntax
        platforms="any",
    )


# %% Run Main install routine
run_setup()
# %% Run Main install routine for radis<0.15 - Cython is not implemented anymore
# try:
#     run_setup(with_binary=False)
# except BuildFailed:
#     import traceback

#     traceback.print_exc()
#     show_message(
#         "WARNING: RADIS C-extension could not be compiled, speedups",
#         " are not enabled.",
#         "Failure information, if any, is above.",
#         "Retrying the build without the C extension now.",
#     )
#     run_setup(with_binary=False)
#     show_message(
#         "WARNING: RADIS C-extension could not be compiled, speedups"
#         " are not enabled.",
#         "Plain-Python build succeeded.",
#     )
