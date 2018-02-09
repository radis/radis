*****************
Spectrum Database
*****************

Spectrum objects can be stored/loaded to/from .spec JSON files using the 
:meth:`~radis.spectrum.spectrum.Spectrum.store` method and the 
:func:`~radis.tools.database.load_spec` function. 

It is also possible to set up a :class:`~radis.tools.database.SpecDatabase` 
which reads all .spec files in a folder. The :class:`~radis.tools.database.SpecDatabase`  
can then be connected to a :class:`~radis.los.factory.SpectrumFactory` so that 
spectra already in the database are not recomputed, and that new calculated spectra 
are stored in the folder

Example::

    db = SpecDatabase(r"path/to/database")     # create or loads database

    db.update()  # in case something changed
    db.see(['Tvib', 'Trot'])   # nice print in console

    s = db.get('Tvib==3000')[0]  # get a Spectrum back
    db.add(s)  # update database (and raise error because duplicate!)



