.. _contributor-hitran-secrets:

Some RADIS features (e.g., fetching HITEMP databases from HITRAN) require a valid
HITRAN account. These credentials are provided to GitHub Actions as secrets:

- ``HITRAN_EMAIL``
- ``HITRAN_PASSWORD``

Why secrets are not available in forks
--------------------------------------

Due to GitHub’s security model, repository secrets are **never exposed to forks**.
If you push changes to your own fork and run workflows there, those workflows will
**not** have access to ``HITRAN_EMAIL`` or ``HITRAN_PASSWORD``.

Options for contributors
------------------------

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
