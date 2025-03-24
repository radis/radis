.. _label_developer_guide:

===============
Developer Guide
===============

.. _label_developer_contribute:

RADIS is an open-source project, and therefore anyone can contribute, whether you already
know spectroscopy or not. The project is organized around a
`GitHub repository <https://github.com/radis/radis/>`__.

You want to become a `contributor <https://github.com/radis/radis/graphs/contributors>`__? Welcome!
First follow the install procedure for developpers below. Then, you can have a look at all the `GitHub opened issues <https://github.com/radis/radis/issues>`__.
and tackle specifically one of the `Good First Issues <https://github.com/radis/radis/contribute>`__.

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
- Code format (see :ref:`Code linting below <label_linting>`).
- Dependency synchronization between ``pyproject.toml`` and ``environment.yml``

.. include:: _test.rst

.. include:: _linting.rst

.. include:: _architecture.rst

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
