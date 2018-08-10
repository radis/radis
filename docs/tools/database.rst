*****************
Spectrum Database
*****************

:class:`~radis.spectrum.spectrum.Spectrum` objects can be stored/loaded to/from 
.spec JSON files using the :meth:`~radis.spectrum.spectrum.Spectrum.store` method 
and the :func:`~radis.tools.database.load_spec` function. 

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

A :class:`~radis.tools.database.SpecDatabase` can also be used to simply
compare the physical and computation parameters of all spectra in folder. 
Indeed, whenever the database is generated, a summary ``.csv`` file 
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

Or using the :meth:`~radis.tools.database.SpecDatabase.get` method with 
no filtering condition::

    db = SpecDatabase('.')
    for s in db.get():
        print(s.name)
        
You can also use dictionary-like methods: :meth:`~radis.tools.database.SpecDatabase.keys`,
:meth:`~radis.tools.database.SpecDatabase.values` and :meth:`~radis.tools.database.SpecDatabase.items`
where Spectrum are returned under a ``{path:Spectrum}`` dictionary.
        
        
Filter spectra that match certain conditions 
--------------------------------------------

If you want to get Spectra in your database that match certain conditions 
(e.g: a particular temperature), you may want to have a look at the 
:meth:`~radis.tools.database.SpecDatabase.get`, 
:meth:`~radis.tools.database.SpecDatabase.get_unique` and 
:meth:`~radis.tools.database.SpecDatabase.get_closest` methods


Updating a database
-------------------

Update all spectra in current folder with a new condition ('author'), making 
use of the :meth:`~radis.tools.database.SpecDatabase.items` method::
                   
    from radis import SpecDatabase
    db = SpecDatabase('.')
    for path, s in db.items():
        s.conditions['author'] = 'me'
        s.store(path, if_exists_then='replace')
                
You may also be interested in the :meth:`~radis.tools.database.SpecDatabase.map` 
method. 


When not to use a Database 
--------------------------

If you simply want to store and reload one :class:`~radis.spectrum.spectrum.Spectrum` 
object, no need to use a database: you better use the :meth:`~radis.spectrum.spectrum.Spectrum.store` 
method and :func:`~radis.tools.database.load_spec` function.

Databases prove useful only when you want to filter precomputed Spectra based on 
certain conditions.  
