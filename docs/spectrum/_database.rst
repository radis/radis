Spectrum Database
=================

:py:class:`~radis.spectrum.spectrum.Spectrum` objects can be stored/loaded to/from
.spec JSON files using the :py:meth:`~radis.spectrum.spectrum.Spectrum.store` method
and the :py:func:`~radis.tools.database.load_spec` function.

.. minigallery:: radis.SpecDatabase

It is also possible to set up a :class:`~radis.tools.database.SpecDatabase`
which reads all .spec files in a folder. The :class:`~radis.tools.database.SpecDatabase`
can then be connected to a :class:`~radis.lbl.factory.SpectrumFactory` so that
spectra already in the database are not recomputed, and that new calculated spectra
are stored in the folder

Example::

    db = SpecDatabase(r"path/to/database")     # create or loads database

    db.update()  # in case something changed
    db.see(['Tvib', 'Trot'])   # nice print in console

    s = db.get('Tvib==3000')[0]  # get a Spectrum back
    db.add(s)  # update database (and raise error because duplicate!)

A :py:class:`~radis.tools.database.SpecDatabase` can also be used to
compare the physical and computation parameters of all spectra in a folder.
Indeed, whenever the database is loaded, a summary ``.csv`` file
is generated that contains all conditions and can be read, for instance,
with Excel.

Example::

    from radis import SpecDatabase
    SpecDatabase(r".")     # this generates a .csv file in the current folder


The examples below show some actions that can be performed on a database:


Iterate over all Spectra in a database
--------------------------------------

Both methods below are equivalent. Directly iterating over the database::

    db = SpecDatabase('.')
    for s in db:
        print(s.name)

Or using the :meth:`~radis.tools.database.SpecList.get` method with
no filtering condition::

    db = SpecDatabase('.')
    for s in db.get():
        print(s.name)

You can also use dictionary-like methods: :meth:`~radis.tools.database.SpecList.keys`,
:meth:`~radis.tools.database.SpecList.values` and :meth:`~radis.tools.database.SpecList.items`
where Spectrum are returned under a ``{path:Spectrum}`` dictionary.


Filter spectra that match certain conditions
--------------------------------------------

If you want to get Spectra in your database that match certain conditions
(e.g: a particular temperature), you may want to have a look at the
:py:meth:`~radis.tools.database.SpecList.get`,
:py:meth:`~radis.tools.database.SpecList.get_unique` and
:py:meth:`~radis.tools.database.SpecList.get_closest` methods


Fit an experimental spectrum against precomputed spectra
--------------------------------------------------------

The :py:meth:`~radis.tools.database.SpecDatabase.fit_spectrum` method
of :py:class:`~radis.tools.database.SpecDatabase` can be used to
return the spectrum of the database that matches the best an experimental
spectrum::

    s_exp = experimental_spectrum(...)
    db = SpecDatabase('...')
    db.fit_spectrum(s_exp)

By default :py:meth:`~radis.tools.database.SpecDatabase.fit_spectrum` uses
the :py:func:`~radis.spectrum.compare.get_residual` function. You can use
an customized function too (below: to get the transmittance)::

    from radis import get_residual
    db.fit_spectrum(s_exp, get_residual=lambda s_exp, s: get_residual(s_exp, s, var='transmittance'))

You don't necessarily need to precompute spectra to fit an experimental spectrum.
You can find an example of :ref:`multi temperature fitting script <label_examples_multitemperature_fit>`
in the Example pages, which shows the evolution of the spectra in real-time. You can get inspiration from there!

Updating a database
-------------------

Update all spectra in current folder with a new condition ('author'), making
use of the :meth:`~radis.tools.database.SpecList.items` method::

    from radis import SpecDatabase
    db = SpecDatabase('.')
    for path, s in db.items():
        s.conditions['author'] = 'me'
        s.store(path, if_exists_then='replace')

You may also be interested in the :py:meth:`~radis.tools.database.SpecList.map`
method.


When not to use a Database
--------------------------

If you simply want to store and reload one :class:`~radis.spectrum.spectrum.Spectrum`
object, no need to use a database: you better use the :meth:`~radis.spectrum.spectrum.Spectrum.store`
method and :func:`~radis.tools.database.load_spec` function.

Databases prove useful only when you want to filter precomputed Spectra based on
certain conditions.
