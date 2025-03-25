.. _label_developer_guide:

===============
Developer Guide
===============

.. _label_developer_contribute:

RADIS is an open-source project, and therefore anyone can contribute, whether you already
know spectroscopy or not. The project is organized around a
`GitHub repository <https://github.com/radis/radis/>`__.

Getting Started
=============

New to open source? Welcome! 👋 Here's how to get started:

1. Look for `Good First Issues <https://github.com/radis/radis/contribute>`__ - these are specifically tagged to help new contributors.
2. The `Documentation TODO List <https://github.com/radis/radis/issues/77>`__ is also a great place to start.

Making Changes
============

We use the `GitHub Flow <https://guides.github.com/introduction/flow/index.html>`__ for all code changes:

1. Fork the repository
2. Create a branch for your feature
3. Make your changes
4. Open a Pull Request (PR)
5. Address review comments
6. Get merged!

For small changes (like fixing typos), you can edit directly on GitHub and open a PR.

Code Style and Linting
--------------------

We use the `Black <https://black.readthedocs.io/en/stable/>`__ and `Flake8 <https://flake8-nb.readthedocs.io/en/latest/>`__ coding style to ensure consistent code formatting.
Code style is checked using CI services
which run automatically on each pull request. **Black** is automatically installed when radis
is set-up in developer mode.

To format any file/files::

    black /path/to/file/or/directory/

You can include Black coding style `directly in some text editors <https://github.com/psf/black#editor-integration>`__

Black coding style can be checked automatically before each commit. For that all you need to do is to run the following command once::

    cd radis
    pre-commit install

On each commit, format will be fixed if it was incorrect. All you need to do is to commit a second time. Example::

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



Keeping Your Fork Updated
----------------------

If you're a regular contributor, keep your fork in sync with the main repository:

1. Add the upstream remote::

    git remote add upstream git://github.com/radis/radis.git
    git fetch upstream

2. Create new branches from the latest upstream version::

    git branch -b [NEW_BRANCH] upstream/develop

3. Update your branch with upstream changes::

    git pull --rebase upstream [NEW_BRANCH]

4. Push to your fork::

    git push -u origin [NEW-BRANCH]

Sources
=======
.. include:: _install.rst

Dependency Management
===================

RADIS uses modern Python packaging with ``pyproject.toml`` for pip installations and ``environment.yml`` for conda installations.
These files are kept in sync through automated tests.

Adding Dependencies
-----------------

When adding new dependencies:

1. Add them to both ``pyproject.toml`` and ``environment.yml``
2. Run the sync test locally: ``python radis/test/test_sync_dependencies.py``
3. Choose the appropriate dependency group:
   - Core dependencies: Add to both files' main section
   - Development tools: Add to ``[project.optional-dependencies].dev`` in ``pyproject.toml``
   - Documentation tools: Add to ``[project.optional-dependencies].docs`` in ``pyproject.toml``

Package Structure
---------------

The package is structured to exclude tests and documentation from the main installation:
- Tests are available when installing with ``pip install -e .[dev]``
- Documentation tools are available with ``pip install -e .[docs]``
- The main installation (``pip install radis`` or ``pip install -e .``) includes only the core package

Update your changes online (push)
---------------------------------

Submit a `Pull Request <https://github.com/radis/radis/pulls>`__ from GitHub.

Online tests will be run automatically. They will check for:

- Physics test cases, to ensure that the code is still working as expected (see :ref:`Test Section <label_dev_test>` to run them locally).
- Code format (see :ref:`Code Style and Linting <label_linting>`).
- Dependency synchronization between ``pyproject.toml`` and ``environment.yml``

.. include:: _test.rst

.. include:: _architecture.rst

Help
----

If you encounter any problem:

1. Check if it's a known issue in our `GitHub Issues <https://github.com/radis/radis/issues>`__
2. If not, please `open a new issue <https://github.com/radis/radis/issues/new/choose>`__ with:
   - A quick summary
   - Steps to reproduce
   - Expected vs actual behavior
   - Any relevant notes or attempted solutions

You can also ask for advice on the community chat:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

.. image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
    :target: https://radis.github.io/slack-invite/
    :alt: Slack
