.. _label_developer_guide:

===============
Developer Guide
===============

.. _label_developer_contribute:

RADIS is an open-source project, and therefore anyone can contribute, whether you already
know spectroscopy or not. The project is organized around a
`GitHub repository <https://github.com/radis/radis/>`__.

You want to become a `contributor <https://github.com/radis/radis/graphs/contributors>`__? Welcome!
First follow the install procedure for developpers below. Then, you can either:

- have a look at all the `GitHub opened issues <https://github.com/radis/radis/issues>`__.
- tackle specifically one of the `Good First Issues <https://github.com/radis/radis/contribute>`__.

Sources
=======
.. include:: _install.rst

Update your changes online (push)
---------------------------------

Submit a `Pull Request <https://github.com/radis/radis/pulls>`__ from GitHub.

Online tests will be run automatically. They will check for:

- Physics test cases, to ensure that the code is still working as expected (see :ref:`Test Section <label_dev_test>` to run them locally).
- Code format (see :ref:`Code linting below <label_linting>`).

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
