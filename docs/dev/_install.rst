.. _label_install:
Install
-------

This is the developper documentation, if you want to access the latest features, modify the code, or contribute.
For the user documentation, see the :ref:`use-only install <label_install_useonly>`.

1. Fork the `RADIS repository <https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_ and use that url for the clone step below. This will make submitting changes easier in the long term for you.

2. In a terminal, run::

    git clone https://github.com/<YOUR GITHUB USERNAME>/radis.git

3. For a best experience, it is recommended to install Radis on top of an Anaconda Python distribution, in an
isolated environment. We recommend to use the faster, local-optimization ``libmamba`` solver for Anaconda.
You can create a radis environment with all dependencies with::

    cd radis
    conda env create --file environment.yml --solver=libmamba

4. Then install Radis in the `radis-env` environment::

    conda activate radis-env
    pip install -e . -v

- The ``-e`` (editable) argument creates a link from the local folder ``./`` into Python
  site-packages.

- For development, you can install additional packages with ``pip install -e .[dev]``
  These include tools for testing and linting your code.

- For documentation development, use ``pip install -e .[docs]`` to install
  additional packages needed for building the documentation.

Dependencies are managed in two files:
- ``pyproject.toml``: for pip installation
- ``environment.yml``: for conda installation

These files are kept in sync through automated tests. When adding new dependencies,
make sure to add them to both files.

To make sure the install worked, run the :ref:`first example <label_first_example>`
from the Quick Start page. Then, you're all set.

Update
------

Whenever you want to update your local copy of Radis, run::

    git pull
    pip install -e .

For development or documentation updates, use ``pip install -e .[dev]`` or ``pip install -e .[docs]`` respectively.

.. note::
    The command `pip install -e .` is different from `pip install radis --upgrade`. The former will install the latest version from the local folder, which may be more recent than the latest release on PyPI.


