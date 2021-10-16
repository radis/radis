# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests::

    pytest      # (in command line, in project folder)

Run only fast tests (i.e: tests that a  'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------


"""

import os
from os.path import basename, exists

import numpy as np
import pytest
from numpy import allclose, linspace

from radis.phys.convert import nm2cm
from radis.spectrum import Spectrum, calculated_spectrum

fig_prefix = basename(__file__) + ": "

# %% Test routines


def test_spectrum_creation_method(*args, **kwargs):
    import pytest

    w = np.linspace(300, 700, 1000)
    k = (np.random.rand(1000) + 0.3) ** 10
    T = np.exp(-k * 1)

    # Good inputs:
    # ... format of quantities :
    Spectrum({"wavelength": w, "abscoeff": k, "transmittance_noslit": T}, wunit="nm")
    Spectrum({"abscoeff": (w, k), "transmittance_noslit": (w, T)}, wunit="nm")

    Spectrum({"wavelength": w, "abscoeff": k, "transmittance_noslit": T}, wunit="nm")

    # ... units:
    Spectrum({"wavespace": w, "abscoeff": k}, wunit="nm")
    Spectrum({"wavelength": w, "abscoeff": k}, wunit="nm")

    # Bad inputs:
    with pytest.raises(AssertionError) as err:
        Spectrum(
            {"wavelength": w, "abscoeff": k, "transmittance_noslit": T}
        )  # wunit not defined
    with pytest.raises(AssertionError) as err:
        Spectrum(
            {"wavenumber": w, "abscoeff": k, "transmittance_noslit": T}
        )  # wunit not defined

    # ... wavespace defined multiple times
    with pytest.raises(AssertionError) as err:
        Spectrum({"wavelength": w, "wavenumber": 1e7 / w, "abscoeff": k}, wunit="nm")
    with pytest.raises(AssertionError) as err:
        Spectrum({"wavelength": w, "wavespace": 1e7 / w, "abscoeff": k}, wunit="nm")

    # ... format of quantities :
    with pytest.raises(AssertionError) as err:
        Spectrum({"wavenumber": w, "abscoeff": np.hstack((k, k))}, wunit="cm-1")
    assert "Input arrays should have the same length" in str(err.value)

    # ... units badly defeined :
    with pytest.raises(AssertionError) as err:
        Spectrum({"wavespace": w, "abscoeff": k})
    assert "waveunit ('nm', 'cm-1'?) has to be defined" in str(err.value)

    with pytest.raises(AssertionError) as err:
        Spectrum({"abscoeff": (w, k)})
    assert "waveunit ('nm', 'cm-1'?) has to be defined" in str(err.value)

    with pytest.raises(AssertionError):
        Spectrum({"wavenumber": w, "abscoeff": k}, wunit="nm")

    with pytest.raises(AssertionError):
        Spectrum({"wavelength": w, "abscoeff": k}, wunit="cm-1")


def test_spectrum_get_methods(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Test all spectrum methods on a Spectrum generated in Specair"""

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    from radis.tools.slit import get_FWHM

    if plot and close_plots:
        import matplotlib.pyplot as plt

        plt.close("all")

    s = load_spec(getTestFile("N2C_specair_380nm.spec"), binary=False)

    # general methods
    if verbose:
        print(s)
    dir(s)

    # access properties
    assert s.get_name() == "N2C_specair_380nm"
    assert all(
        s.get_radiance_noslit(Iunit="W/m2/sr/nm")
        == s.get("radiance_noslit", Iunit="W/m2/sr/nm")[1]
    )
    assert all(nm2cm(s.get_wavelength(medium="vacuum")) == s.get_wavenumber())
    assert np.isclose(s.get_power(unit="W/cm2/sr"), 2631.6288408588148)
    assert s.get_waveunit() == "nm"
    assert np.isclose(
        s.get_power(unit="W/cm2/sr"),
        s.get_integral("radiance_noslit", Iunit="W/cm2/sr/nm"),
    )
    assert s.get_conditions()["Tgas"] == 1500
    assert len(s.get_vars()) == 2
    assert s.is_at_equilibrium() == False
    assert s.is_optically_thin() == False

    # Check applied slit has the correct width
    s.apply_slit(0.5)
    wslit, Islit = s.get_slit()
    wstep = np.diff(wslit)[0]
    assert np.isclose(get_FWHM(*s.get_slit()), 0.5, atol=1.1 * wstep)

    if plot:
        s.plot_slit()

    if verbose:
        print("Tested Spectrum methods:")
        print("...print(Spectrum)")
        print(".. get_name()")
        print(".. get_radiance_noslit() vs get()")
        print(".. get_wavelength() vs get_wavenumber")
        print(".. get_power()")
        print(".. get_waveunit()")
        print(".. get_power() vs get_integral()")
        print(".. get_conditions()")
        print(".. get_vars()")
        print(".. is_at_equilibrium()")
        print(".. is_optically_thin()")
        print(".. get_slit()")


@pytest.mark.fast
def test_trimming(verbose=True, *args, **kwargs):
    """Test :py:meth:`radis.spectrum.spectrum.Spectrum.trim`"""

    from radis.misc.arrays import count_nans

    # trim both sides
    w = np.linspace(300, 600)
    T = np.ones_like(w)
    T[:15] = np.nan
    T[-10:] = np.nan
    s = Spectrum.from_array(w, T, "transmittance_noslit", wunit="nm", unit="")

    Nnans = count_nans(T)
    assert len(s) == 50  # before trimming
    s.trim()
    assert len(s) == 50 - Nnans

    # trim one side
    w = np.linspace(300, 600)
    T = np.ones_like(w)
    T[:15] = np.nan
    s = Spectrum.from_array(w, T, "transmittance_noslit", wunit="nm", unit="")

    Nnans = count_nans(T)
    assert len(s) == 50  # before trimming
    s.trim()
    assert len(s) == 50 - Nnans


@pytest.mark.fast
def test_copy(verbose=True, *args, **kwargs):
    """Test that a Spectrum is correctly copied

    We compare a Spectrum that has:
    - all available spectral quantities
    - a slit
    - many calculation conditions
    - no populations
    - no lines
    """

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))

    s.update()
    s.apply_slit(1.5)
    s2 = s.copy()

    # Test spectrum in general
    assert s == s2
    assert s is not s2

    # Test all quantities in detail
    for var in s._q.keys():
        assert np.allclose(s._q[var], s2._q[var], equal_nan=True)
        assert not (s._q[var] is s2._q[var])

    if verbose:
        print("Tested that s2 == s (but s2 is not s) after Spectrum copy")


def test_populations(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    """Test that populations in a Spectrum are correctly read"""

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()  # dont get stuck with Matplotlib if executing through pytest

    if plot and close_plots:
        import matplotlib.pyplot as plt

        plt.close("all")

    import pytest

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    # get a spectrum
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))

    # Get all populations
    pops = s.get_populations(molecule="CO", isotope=1)
    assert np.isclose(pops["Ia"], 0.986544)

    # Get vib levels
    df_vib = s.get_vib_levels()
    assert len(df_vib) == 49

    # Get rovib levels (dont exist: should fail!)
    with pytest.raises(KeyError):  # expected behavior
        s.get_rovib_levels()

    if plot:
        s.plot_populations()

    if verbose:
        print("test_populations: OK")


def test_store_functions(verbose=True, *args, **kwargs):
    """Test some store / retrieve functions"""

    from radis.spectrum.models import transmittance_spectrum
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    temp_file = "test_radis_tempfile_transmittance.txt"
    assert not exists(temp_file)

    # Test that the transmittance stored on .txt and loaded again match
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update()
    try:
        s.savetxt(temp_file, "transmittance_noslit", wunit="nm_vac")
        w, T = np.loadtxt(temp_file).T
    finally:
        os.remove(temp_file)

    s2 = transmittance_spectrum(w, T, wunit="nm_vac")
    assert s.compare_with(s2, spectra_only="transmittance_noslit", plot=False)

    # TODO: add test that ensures we can load a binary file without binary=True
    # (and the warning should be captured)

    return True


@pytest.mark.fast
def test_intensity_conversion(verbose=True, *args, **kwargs):
    """Test conversion of intensity cm-1 works:

    - conversion of mW/sr/cm2/nm -> mW/sr/cm2/cm-1

    """

    from radis import planck, planck_wn

    w_nm = linspace(300, 3000)
    w_cm = nm2cm(w_nm)
    I_nm = planck(w_nm, T=6000, unit="mW/sr/cm2/nm")

    s = calculated_spectrum(
        w_nm,
        I_nm,
        wunit="nm_vac",
        Iunit="mW/sr/cm2/nm",
    )
    # mW/sr/cm2/nm -> mW/sr/cm2/cm-1
    w, I = s.get("radiance_noslit", Iunit="mW/sr/cm2/cm-1")
    I_cm = planck_wn(w_cm, T=6000, unit="mW/sr/cm2/cm-1")
    assert allclose(I_cm, I, rtol=1e-3)


# def test_emissivity_conversion(verbose=True, *args, **kwargs):
#    ''' Test conversion of intensity cm-1 works:
#
#    - conversion of mW/sr/cm2/nm -> mW/sr/cm2/cm-1
#
#    '''
#
#    from radis.phys.units import convert_emi2nm
#    from radis.test.utils import getTestFile
#    from radis.tools.database import load_spec
#
#    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
#    s.update()
#
#    # mW/sr/cm3/nm -> mW/sr/cm3/cm-1
#
#    w, I = s.get('emissivity_noslit', Iunit='mW/sr/cm3/cm-1')
#    I_cm = convert_emi2nm(I, 'mW/sr/cm3/cm-1', 'mW/sr/m3/µm')
#

# TODO: finish implementing emissivity_conversino above


def test_rescaling_function(verbose=True, *args, **kwargs):
    """Test rescaling functions"""

    from radis.test.utils import getTestFile

    s = Spectrum.from_txt(
        getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt"),
        quantity="radiance_noslit",
        wunit="nm",
        unit="mW/cm2/sr/µm",  # Specair units: mW/cm2/sr/µm
        conditions={
            "Tvib": 3000,
            "Trot": 1200,
            "path_length": 1,  # arbitrary
            "medium": "air",
        },
        populations={"molecules": {"N2C": 1e13}},  # arbitrary
        # (just an example)
    )
    s.update(optically_thin=True)
    s.rescale_path_length(10)

    assert np.isclose(s.get_radiance_noslit(Iunit="mW/cm2/sr/nm")[0], 352.57305783248)


def test_resampling_function(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Test resampling functions

    Get a Spectrum calculated in cm-1, then resample on a smaller range in cm-1,
    and in approximately the same range (but in nm). Check that all 3 overlap
    """
    # %%
    from radis.spectrum import get_residual_integral
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    if plot and close_plots:
        import matplotlib.pyplot as plt

        plt.close("all")

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.name = "original"

    # Test resampling without changing anything
    s_cm = s.resample(s.get_wavenumber(), "cm-1", inplace=False)
    s_nm = s.resample(s.get_wavelength(), "nm", inplace=False)
    assert s == s_cm
    assert np.allclose(s.get("abscoeff", wunit="nm"), s_nm.get("abscoeff", wunit="nm"))

    # Test resampling on new ranges
    s2 = s.copy()
    s2b = s.copy()
    s3 = s.copy()
    s2.resample(np.linspace(4500, 4700, 10000), unit="nm_vac")
    s2b.resample(np.linspace(4500, 4700, 10000), unit="nm")
    s3.resample(np.linspace(2127.2, 2227.7, 10000), unit="cm-1")
    s2.name = "resampled in nm (vacuum)"
    s2b.name = "resampled in nm (air)"
    s3.name = "resampled in cm-1"

    s.compare_with(
        s2,
        plot=plot,
        title="Residual: {0:.2g}".format(
            get_residual_integral(s, s2, "abscoeff", ignore_nan=True)
        ),
    )
    s.compare_with(
        s2b,
        plot=plot,
        title="Residual: {0:.2g}".format(
            get_residual_integral(s, s2b, "abscoeff", ignore_nan=True)
        ),
    )
    s.compare_with(
        s3,
        plot=plot,
        title="Residual: {0:.2g}".format(
            get_residual_integral(s, s3, "abscoeff", ignore_nan=True)
        ),
    )

    assert get_residual_integral(s, s2, "abscoeff", ignore_nan=True) < 1e-3
    assert get_residual_integral(s, s2b, "abscoeff", ignore_nan=True) < 1e-3
    assert get_residual_integral(s, s3, "abscoeff", ignore_nan=True) < 1e-5


def test_resampling_nan_function(verbose=True, *args, **kwargs):
    """Test resampling functions, for Spectra with nans"""
    from radis import get_residual
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    plot = True

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True).crop(
        2170, 2180, "cm-1"
    )
    s.rescale_path_length(10)  # cm
    s.name = "original"

    # We want to resample on a small range:
    w_exp = s.get_wavenumber()[700:-700][::5]
    # ... add a shift:
    w_exp += np.diff(s.get_wavenumber())[0] * 2 / 3

    s.update("transmittance_noslit")
    s.apply_slit(3, "nm")

    assert s.has_nan()

    # s2.resample(w_exp)
    s2 = s.resample(w_exp, inplace=False)
    s2.name = "resampled"

    assert get_residual(s2, s, "transmittance") < 1e-6

    if plot:
        from radis import plot_diff

        plot_diff(s, s2, show_points=True)


@pytest.mark.fast
def test_noplot_different_quantities(*args, **kwargs):
    """Prevents User Errors: Ensures an error is raised if plotting different
    quantities on the same graph"""

    import matplotlib.pyplot as plt

    plt.ion()

    from radis import load_spec
    from radis.test.utils import getTestFile

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update()
    s.plot("abscoeff", nfig="test_noplot_different_quantities")
    with pytest.raises(ValueError):  # expect an error
        s.plot("emisscoeff", nfig="same")

    plt.close("test_noplot_different_quantities")


@pytest.mark.fast
def test_plot_by_parts(plot=True, *args, **kwargs):
    """Test :py:func:`~radis.misc.plot.split_and_plot_by_parts`
    and plot_by_parts=True in :py:meth:`~radis.spectrum.spectrum.Spectrum.plot`
    """

    import matplotlib.pyplot as plt

    plt.ion()

    from radis import load_spec
    from radis.test.utils import getTestFile

    load_spec(getTestFile("CO2_measured_spectrum_4-5um.spec"), binary=True).plot(
        plot_by_parts=True, nfig="plot by parts (non continuous spectrum)"
    )

    if not plot:
        plt.close("plot by parts (non continuous spectrum)")


@pytest.mark.fast
def test_normalization(*args, **kwargs):

    from radis import Radiance, load_spec
    from radis.test.utils import getTestFile

    # Generate the equivalent of an experimental spectrum
    s = load_spec(getTestFile(r"CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update()  # add radiance, etc.
    s.apply_slit(0.5)  # nm
    s = Radiance(s)

    # Test normalization
    assert s.units["radiance"] != ""
    s.normalize()
    assert s.max() != 1
    s.normalize(inplace=True)
    assert s.max() == 1
    assert s.normalize().units["radiance"] == ""

    s2 = s.normalize(normalize_how="area")
    assert np.isclose(s2.get_integral("radiance", wunit=s2.get_waveunit()), 1)

    #
    s3 = s.normalize(wrange=((2125, 2150)), normalize_how="area")
    assert np.isclose(
        s3.crop(2125, 2150).get_integral("radiance", wunit=s3.get_waveunit()), 1
    )


@pytest.mark.fast
def test_sort(*args, **kwargs):
    """Test :py:meth:`~radis.spectrum.spectrum.Spectrum.sort`"""

    from radis import load_spec
    from radis.test.utils import getTestFile

    s_exp = load_spec(getTestFile("CO2_measured_spectrum_4-5um.spec"), binary=True)

    from radis.misc.arrays import is_sorted

    assert not is_sorted(s_exp.get("radiance")[0])
    assert is_sorted(s_exp.sort().get("radiance")[0])


# %%


def _run_testcases(
    plot=True,
    close_plots=False,
    verbose=True,
    debug=False,
    warnings=True,
    *args,
    **kwargs
):
    """Test procedures

    Parameters
    ----------

    debug: boolean
        swamps the console namespace with local variables. Default ``False``

    """

    # Test all Spectrum methods
    # -------------------------
    test_spectrum_creation_method(*args, **kwargs)
    test_spectrum_get_methods(
        debug=debug,
        verbose=verbose,
        plot=plot,
        close_plots=close_plots,
        *args,
        **kwargs
    )
    test_copy(verbose=verbose, *args, **kwargs)
    test_trimming(*args, **kwargs)
    test_populations(
        verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs
    )
    test_store_functions(verbose=verbose, *args, **kwargs)

    # Test populations
    # ----------
    test_populations(verbose=verbose, *args, **kwargs)

    # Test conversion of intensity cm-1 works
    # -------------
    test_intensity_conversion(debug=debug, verbose=verbose, *args, **kwargs)

    # Test updating / rescaling functions (no self absorption)
    # ---------
    test_rescaling_function(debug=debug, *args, **kwargs)

    test_resampling_function(
        debug=debug, plot=plot, close_plots=close_plots, *args, **kwargs
    )
    test_resampling_nan_function(verbose=verbose, *args, **kwargs)

    test_normalization(*args, **kwargs)

    # Test plot firewalls:
    test_noplot_different_quantities(*args, **kwargs)

    # Test plot by parts
    test_plot_by_parts(plot=plot, *args, **kwargs)
    test_sort()

    return True


if __name__ == "__main__":
    print(("Test spectrum: ", _run_testcases(debug=False, close_plots=False)))
