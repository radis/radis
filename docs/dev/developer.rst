.. _label_developer_guide:

===============
Developer Guide
===============

.. _label_developer_contribute:

RADIS is an open-source project, and therefore anyone can contribute, whether you already
know spectroscopy or not. The project is organized around a
`GitHub repository <https://github.com/radis/radis/>`__.

Installation for developers
=======
.. include:: _install.rst


GitHub repository
=======
.. include:: _github.rst


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

Testing RADIS
=============
.. include:: _test.rst

Architecture
============
.. include:: _architecture.rst

Help
===================

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
