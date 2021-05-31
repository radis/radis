.. _label_install:

Install
-------

You can either install RADIS from `pip`, the Python package manager. But if
you want to access the latest features, or modify the code and contribute,
we suggest that you Fork the source code from `GitHub <https://github.com/radis/radis>`_.

**Use-only version** : cant modify the code

In a terminal, run::

    pip install --user radis -v

The 'pip' module has to be installed (by default if you've installed Python
with Anaconda).

If you later want to update the package::

    pip install radis --upgrade

**Developer version**: to modify the code and contribute to the
project.

We suggest that you fork the `RADIS repository <https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_ and use that url for the clone step below. This will make submitting changes easier in the long term for you:

In a terminal, run::

    git clone https://github.com/<GITHUB USERNAME>/radis.git
    cd radis
    pip install -e .[dev] -v

- The ``-e`` (editable) argument creates a link from the local folder ``./`` into Python
  site-packages.

- ``[dev]`` is optional, but installs additional developer packages that may be useful for testing and
  linting your code.

To make sure the install worked, run the :ref:`first example <label_first_example>`
from the Quick Start page. Then, you're all set.



For a best experience, it is recommended to install Radis on an Anaconda Python distribution, in an
isolated environment. You can create a radis environment with all dependencies with::

    cd radis
    conda env create --file environment.yml

Then install radis in the `radis` environment::

    conda activate radis
    pip install -e .



Test your changes
-----------------

If you want to modify the source code, you need to ensure that you don't break
any of the existing tests.
Refer to the :ref:`Test Section <label_dev_test>` to learn how to run the
tests locally.




Update your changes
-------------------

Submit a `Pull Request <https://github.com/radis/radis/pulls>`__ from GitHub.

Online tests will be run automatically. They will check for:

- Physical test cases
- Code format (see Code linting below)


Code linting
------------

Radis follows `Black <https://black.readthedocs.io/en/stable/>`__ style for code linting to
maintain consistent coding style across modules. Code style is checked using CI services
which run automatically on each pull request. **Black** is automatically installed when radis
is set-up in developer mode.

To format any file/files::

    black /path/to/file/or/directory/

You can include Black coding style `directly in some text editors <https://github.com/psf/black#editor-integration>`__

Alternatively, Black coding style can be checked automatically before each commit. For that all you need to do is to run the following command once::

    cd radis
    pre-commit install

On each commit, format will be fixed if it was incorrect. All you need to do is to commit a second time. Exemple::

    $ git commit -am "test"
    black....................................................................Failed
    - hook id: black
    - files were modified by this hook

    reformatted [ALL FAILING FILES]
    All done!
    1 file reformatted.

    $ git commit -am "test"
    black....................................................................Passed
    [develop XXX] test
     1 file changed, 1 insertion(+)

Note that pre-commit will always require you to commit again after a test was failed, because `it's safer <https://github.com/pre-commit/pre-commit/issues/532>`__. If for any reason you want to skip formatting you can commit with the ``--no-verify`` `argument <https://git-scm.com/docs/git-commit>`__.




Update
------

With Pip you can keep the package up-to-date with::

    pip install radis --upgrade

If using the latest developer version (cloned from `GitHub <https://github.com/radis/radis>`_ and installed with `pip install -e .[dev]`), use git to `pull` the latest changes.

Help
----

If you encounter any problem, please `report an Issue <https://github.com/radis/radis/issues?utf8=%E2%9C%93&q=is%3Aissue>`_ on GitHub.

You can also ask for advice on the `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__
or the community chat:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

.. image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
    :target: https://radis.github.io/slack-invite/
    :alt: Slack
