# -*- coding: utf-8 -*-
"""
Created on Mon May  7 17:34:52 2018

@author: erwan
"""

from os.path import exists
from shutil import rmtree

import matplotlib.pyplot as plt
import pytest

from radis.lbl import SpectrumFactory
from radis.misc.printer import printm
from radis.test.utils import getTestFile, setup_test_line_databases


def test_retrieve_from_database(
    plot=False, verbose=True, warnings=True, *args, **kwargs
):
    """Test autoretrieve from a database:

    first generate an empty :py:class:`~radis.tools.database.SpecDatabase`
    associated to a :py:class`~radis.lbl.factory.SpectrumFactory`,
    then calculate a first spectrum, then calculate it again and make sure
    it is retrieved from the database
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    temp_database_name = "temp_spec_database"

    try:
        if verbose:
            printm(">>> test_retrieve_from_database")

        assert not exists(temp_database_name)

        setup_test_line_databases()  # add HITEMP-CO2-TEST in ~/radis.json if not there

        sf = SpectrumFactory(
            2284.2,
            2284.6,
            wstep=0.001,  # cm-1
            pressure=20 * 1e-3,  # bar
            cutoff=0,
            path_length=0.1,
            mole_fraction=400e-6,
            molecule="CO2",
            isotope=[1],
            medium="vacuum",
            truncation=5,
            export_populations="rovib",
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.init_databank(
            "HITEMP-CO2-TEST",
            load_columns="noneq",
        )  # unlike load_databank, will automatically be built when needed
        db = sf.init_database(temp_database_name, autoretrieve=True)

        assert len(db) == 0

        # Calculate a first spectrum
        s1 = sf.non_eq_spectrum(2000, 1000)

        assert len(db) == 1
        # Note that the Spectrum in database is not s1 because populations
        # are not stored by default

        # hack: now change the `autoretrieve` status to `force` to make sure
        # the spectrum is retrieved. Else, an error will be raised
        sf.autoretrievedatabase = "force"

        # Calculate spectrum under the same conditions
        s2 = sf.non_eq_spectrum(2000, 1000)

        # Note that s1 == s2 won't work because populations are not stored
        # by default in the database
        assert s1.compare_with(s2, spectra_only=True, plot=plot)

        return True

    finally:
        rmtree(temp_database_name)


@pytest.mark.fast
def test_ignore_cached_files():
    """
    Previous implementation of RADIS saved the cached h5 files generated while reading the
    dataset in the same directory from where the data was being read. Using a wildcard input
    such as `path = "cdsd_hitemp_09_frag*"` in such case led to the cached files present
    in directory to also being loaded and treated as the dataset files. This resulted in
    an error due to the differences in the way data is stored in h5 files versus in dataset
    files such as par, txt, etc.

    Reference: `https://github.com/radis/radis/issues/121`
    """

    sf = SpectrumFactory(wavenum_min=2000, wavenum_max=3000, pressure=1)

    file_dir = getTestFile("cdsd_hitemp_09_fragment.txt")
    test_file = file_dir[:-8] + "*"
    sf.load_databank(path=test_file, format="cdsd-hitemp", parfuncfmt="hapi")

    try:
        sf.load_databank(path=test_file, format="cdsd-hitemp", parfuncfmt="hapi")
    except UnicodeDecodeError as err:
        raise UnicodeDecodeError(
            "Couldn't load database the 2nd time. This may be due to cache files trying to be read as normal files"
        ) from err


@pytest.mark.fast
def test_ignore_irrelevant_files(*args, **kwargs):
    """
    Implemented in https://github.com/radis/radis/pull/185

    - Test that calculating with cache files outside the spectral range they
    are properly ignored
    """

    from radis.misc.warning import EmptyDatabaseError, IrrelevantFileWarning

    # Regenerate .h5 cache file
    sf = SpectrumFactory(wavenum_min=2280, wavenum_max=2290)
    test_file = getTestFile("cdsd_hitemp_09_fragment.txt")
    sf.load_databank(
        path=test_file, format="cdsd-hitemp", parfuncfmt="hapi", db_use_cached="regen"
    )
    assert exists(test_file.replace(".txt", ".h5"))
    # Also note that there was no EmptyDatabaseError : file was properly loaded!

    # Load same .h5 cache file in another spectral range
    # ... Expect an error in use_cached = 'force' (file properly detected as
    # irrelevant)
    sf2 = SpectrumFactory(wavenum_min=100000, wavenum_max=100002)
    with pytest.raises(IrrelevantFileWarning):
        sf2.load_databank(
            path=test_file,
            format="cdsd-hitemp",
            parfuncfmt="hapi",
            db_use_cached="force",
        )

    # Again without 'force'
    # ... Expect no IrrelevantFile error, however range should be empty.
    with pytest.raises(EmptyDatabaseError):
        sf3 = SpectrumFactory(wavenum_min=100000, wavenum_max=100002)
        sf3.load_databank(path=test_file, format="cdsd-hitemp", parfuncfmt="hapi")


@pytest.mark.fast
def test_custom_abundance(verbose=True, plot=False, *args, **kwargs):
    """Test calculatinos with custom abundances

    set the terrestrial abundance to False to
    force computation of linestrength from Einstein coefficients instead
    of rescaling parameters

    :py:meth:`~radis.lbl.base.DatabankLoader.get_abundance`,
    :py:meth:`~radis.lbl.base.DatabankLoader.set_abundance`
    """

    sf = SpectrumFactory(
        2284.2,
        2284.6,
        wstep=0.001,  # cm-1
        pressure=20 * 1e-3,  # bar
        cutoff=0,
        path_length=0.1,
        mole_fraction=400e-6,
        molecule="CO2",
        isotope="1,2",
        medium="vacuum",
        truncation=5,
        verbose=0,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITEMP-CO2-TEST", load_columns="all")

    #%% Compute a spectrum ; set the terrestrial abundance to False to
    # force computation of linestrength from Einstein coefficients instead
    # of rescaling parameters
    s = sf.eq_spectrum(2000)
    sf.molparam.terrestrial_abundances = False
    s2 = sf.eq_spectrum(2000)
    assert (
        s != s2
    )  # compare spectra; we should see small differences due to different methods
    assert s.compare_with(s2, "abscoeff", rtol=0.5e-2, plot=False)
    # assert np.allclose(2 * s2.get("abscoeff")[1], s.get("abscoeff")[1])

    if plot:
        from radis import plot_diff

        plot_diff(s, s2)

    #%% Now restart by changing the abundance itself ; as a user would do

    # ... check abundance get/set functions work well
    abundance12_0 = sf.get_abundance("CO2", [1, 2])

    if verbose:
        print("Abundance: CO2[1]", sf.get_abundance("CO2", 1))
        print("Abundance: CO2[1,2]", sf.get_abundance("CO2", [1, 2]))

    sf.set_abundance("CO2", [1, 2], abundance12_0 / 3)
    assert all(sf.get_abundance("CO2", [1, 2]) == abundance12_0 / 3)
    sf.set_abundance("CO2", 1, abundance12_0[0] / 3)
    sf.set_abundance("CO2", 2, abundance12_0[1] / 3)
    assert all(sf.get_abundance("CO2", [1, 2]) == abundance12_0 / 3)

    # Below we compute a spectrum with the abundance divided by 3; and then
    # check that the resulting absorption coefficient is 3x smaller

    s3 = sf.eq_spectrum(2000)
    assert s.compare_with((3 * s3.take("abscoeff")), "abscoeff", rtol=0.5e-2, plot=True)


def _run_testcases(verbose=True, plot=False):

    test_retrieve_from_database(plot=plot, verbose=verbose)
    test_ignore_cached_files()
    test_ignore_irrelevant_files(verbose=verbose)
    test_custom_abundance()


if __name__ == "__main__":
    _run_testcases()
