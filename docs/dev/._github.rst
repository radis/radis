.. _contributor-github:


Getting Started
----------------------

New to open source? Welcome! 👋 Here's how to get started:

1. Look for `Good First Issues <https://github.com/radis/radis/contribute>`__ - these are specifically tagged to help new contributors.
2. The `Documentation TODO List <https://github.com/radis/radis/issues/77>`__ is also a great place to start.

Making Changes
----------------------

We use the `GitHub Flow <https://guides.github.com/introduction/flow/index.html>`__ for all code changes:

1. Fork the repository
2. Create a branch for your feature
3. Make your changes
4. Open a Pull Request (PR)
5. Address review comments
6. Get merged!

For small changes (like fixing typos), you can edit directly on GitHub and open a PR.

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


Contributing Securely: HITRAN Credentials
----------------------

Some RADIS features (e.g., fetching HITEMP databases from HITRAN) require a valid
HITRAN account. These credentials are provided to GitHub Actions as secrets:

- ``HITRAN_EMAIL``
- ``HITRAN_PASSWORD``

Due to GitHub’s security model, repository secrets are **never exposed to forks**,
see also the `GitHub documentation <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`__ .
If you push changes to your own fork and run workflows there, those workflows will
**not** have access to ``HITRAN_EMAIL`` or ``HITRAN_PASSWORD``.

Option 1: Continue working in your fork (with your own secrets)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer to work in your fork:

1. Navigate to **Settings → Secrets and variables → Actions** on your fork’s GitHub page.
2. Add the following repository secrets:

   - ``HITRAN_EMAIL`` → your HITRAN account email
   - ``HITRAN_PASSWORD`` → your HITRAN account password

3. Workflows in your fork will then use *your* credentials when fetching HITEMP databases.

.. note::

   You must have a valid HITRAN registration before adding these.


Option 2: Work on a branch inside the main ``radis/radis`` repository. Applies to regular contributors.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Request **Write**, **Maintain**, or **Admin** access from the RADIS maintainers.
2. Push your work to a branch inside the main repository, for example
   ``feature/my-new-test``.
3. When workflows run there, they will have access to the official secrets.

RADIS will use these environment variables for authentication.


Code Style and Linting
--------------------

We use the `Black <https://black.readthedocs.io/en/stable/>`__ and `Flake8 <https://flake8-nb.readthedocs.io/en/latest/>`__ coding style to ensure consistent code formatting.
Code style is checked using CI services which run automatically on each pull request.
The linting is executed by the file linting.yml in the .github/workflows folder.
 **Black** is automatically installed when radis is set-up in developer mode.

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
