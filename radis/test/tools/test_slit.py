# -*- coding: utf-8 -*-
"""
Test slit methods

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

    pytest -m fast


-------------------------------------------------------------------------------

"""


from os.path import basename
from warnings import catch_warnings, filterwarnings

import matplotlib.pyplot as plt
import numpy as np
import pytest
from numpy import abs, cos, linspace, pi, sqrt, tan

from radis.lbl.factory import SpectrumFactory
from radis.misc.arrays import nantrapz
from radis.misc.plot import fix_style, set_style
from radis.misc.printer import printm
from radis.phys.convert import dcm2dnm, dnm2dcm
from radis.phys.units import is_homogeneous

# from radis.lbl import SpectrumFactory
from radis.spectrum.models import calculated_spectrum, transmittance_spectrum
from radis.spectrum.spectrum import _cut_slices
from radis.test.utils import setup_test_line_databases
from radis.tools.database import load_spec
from radis.tools.slit import (
    get_effective_FWHM,
    get_FWHM,
    import_experimental_slit,
    offset_dilate_slit_function,
)

fig_prefix = basename(__file__) + ": "

# %% --------------------------------------------------------------------------
#  Test cases
# -----------------------------------------------------------------------------


def _clean(plot, close_plots):
    if plot:
        plt.ion()  # dont get stuck with Matplotlib if executing through pytest
        set_style()
        if close_plots:
            plt.close("all")


@pytest.mark.fast
def test_all_slit_shapes(
    FWHM=0.4, verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Test all slit generation functions and make sure we get the expected FWHM"""

    _clean(plot, close_plots)

    # get spectrum
    from radis.spectrum.spectrum import Spectrum
    from radis.test.utils import getTestFile

    s = Spectrum.from_txt(
        getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt"),
        quantity="radiance_noslit",
        wunit="nm",
        unit="mW/cm2/sr/µm",
    )
    wstep = np.diff(s.get_wavelength())[0]

    # Plot all slits
    # ... gaussian
    s.apply_slit(FWHM, unit="nm", shape="gaussian", plot_slit=plot)
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2 * wstep)

    # ... triangular
    s.apply_slit(FWHM, unit="nm", shape="triangular", plot_slit=plot)
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2 * wstep)

    # ... trapezoidal
    s.apply_slit(
        (FWHM * 0.9, FWHM * 1.1), unit="nm", shape="trapezoidal", plot_slit=plot
    )
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2 * wstep)

    #    # ... trapezoidal
    #    s.apply_slit(FWHM, unit='nm', shape='trapezoidal', plot_slit=plot, norm='max')
    #    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=1.1*wstep)

    # ... experimental
    s.apply_slit(getTestFile("slitfunction.txt"), unit="nm", plot_slit=plot)
    assert np.isclose(get_effective_FWHM(*s.get_slit()), FWHM, atol=0.01)
    # note that we're applying a slit function measured at 632.5 nm to a Spectrum
    # at 4.7 µm. It's just good for testing the functions

    #    # ... experimental, convolve with max
    #    s.apply_slit(getTestFile('slitfunction.txt'), unit='nm', norm_by='max', plot_slit=plot)
    #    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=1.1*wstep)

    if verbose:
        print("\n>>> _test_all_slits yield correct FWHM (+- wstep) : OK\n")

    return True  # nothing defined yet


@pytest.mark.fast
def test_slit_unit_conversions_spectrum_in_cm(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Test that slit is consistently applied for different units

    Assert that:

    - calculated FWHM is the one that was applied

    """

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    _clean(plot, close_plots)

    # %% Get a Spectrum (stored in cm-1)
    s_cm = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)

    s_cm.rescale_mole_fraction(1)  # just because it makes better units
    s_cm.update()
    wstep = s_cm.conditions["wstep"]

    assert s_cm.get_waveunit() == "cm-1"  # ensures it's stored in cm-1

    for shape in ["gaussian", "triangular"]:

        # Apply slit in cm-1
        slit_cm = 2
        s_cm.name = "Spec in cm-1, slit {0:.2f} cm-1".format(slit_cm)
        s_cm.apply_slit(slit_cm, unit="cm-1", shape=shape, mode="same")
        # ... mode=same to keep same output length. It helps compare both Spectra afterwards
        # in cm-1 as that's s.get_waveunit()
        fwhm = get_FWHM(*s_cm.get_slit())
        assert np.isclose(slit_cm, fwhm, atol=2 * wstep)

        # Apply slit in nm this time
        s_nm = s_cm.copy()
        w_cm = s_nm.get_wavenumber()
        slit_nm = dcm2dnm(slit_cm, w_cm[len(w_cm) // 2])
        s_nm.name = "Spec in cm-1, slit {0:.2f} nm".format(slit_nm)
        s_nm.apply_slit(slit_nm, unit="nm", shape=shape, mode="same")

        plotargs = {}
        if plot:
            plotargs["title"] = "test_slit_unit_conversions: {0} ({1} cm-1)".format(
                shape, slit_cm
            )
        s_cm.compare_with(
            s_nm,
            spectra_only="radiance",
            rtol=1e-3,
            verbose=verbose,
            plot=plot,
            **plotargs
        )


@pytest.mark.fast
def test_slit_unit_conversions_spectrum_in_nm(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Test that slit is consistently applied for different units

    Assert that:

    - calculated FWHM is the one that was applied

    """

    from radis.spectrum.spectrum import Spectrum
    from radis.test.utils import getTestFile

    _clean(plot, close_plots)

    # %% Get a Spectrum (stored in nm)

    s_nm = Spectrum.from_txt(
        getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt"),
        quantity="radiance_noslit",
        wunit="nm",
        unit="mW/cm2/sr/µm",
        conditions={"self_absorption": False},
    )

    with catch_warnings():
        filterwarnings(
            "ignore", "Condition missing to know if spectrum is at equilibrium:"
        )
        # just because it makes better units
        s_nm.rescale_path_length(1, 0.001)

    wstep = np.diff(s_nm.get_wavelength())[0]

    assert s_nm.get_waveunit() == "nm"  # ensures it's stored in cm-1

    for shape in ["gaussian", "triangular"]:

        # Apply slit in nm
        slit_nm = 0.5
        s_nm.name = "Spec in nm, slit {0:.2f} nm".format(slit_nm)
        s_nm.apply_slit(slit_nm, unit="nm", shape=shape, mode="same")
        # ... mode=same to keep same output length. It helps compare both Spectra afterwards
        # in cm-1 as that's s.get_waveunit()
        fwhm = get_FWHM(*s_nm.get_slit())
        assert np.isclose(slit_nm, fwhm, atol=2 * wstep)

        # Apply slit in nm this time
        s_cm = s_nm.copy()
        w_nm = s_nm.get_wavelength()
        slit_cm = dnm2dcm(slit_nm, w_nm[len(w_nm) // 2])
        s_cm.name = "Spec in nm, slit {0:.2f} cm-1".format(slit_cm)
        s_cm.apply_slit(slit_cm, unit="cm-1", shape=shape, mode="same")

        plotargs = {}
        if plot:
            plotargs["title"] = "test_slit_unit_conversions: {0} ({1} nm)".format(
                shape, slit_nm
            )
        s_nm.compare_with(
            s_cm,
            spectra_only="radiance",
            rtol=1e-3,
            verbose=verbose,
            plot=plot,
            **plotargs
        )

    # %%


def test_convoluted_quantities_units(*args, **kwargs):
    """Test that units are correctly convoluted after convolution"""

    from radis.test.utils import getTestFile

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"), binary=True)

    s.update(verbose=False)

    assert s.units["radiance_noslit"] == "mW/cm2/sr/nm"
    assert s.units["transmittance_noslit"] == ""

    s.apply_slit(0.5, norm_by="area", verbose=False)

    assert s.units["radiance"] == "mW/cm2/sr/nm"
    assert s.units["transmittance"] == ""

    s.apply_slit(0.5, norm_by="max", verbose=False)

    assert is_homogeneous(s.units["radiance"], "mW/cm2/sr")
    assert s.units["transmittance"] == "nm"  # whatever that means


@pytest.mark.fast
def test_against_specair_convolution(
    plot=True, close_plots=True, verbose=True, debug=False, *args, **kwargs
):

    _clean(plot, close_plots)

    # Test
    from radis.test.utils import getTestFile

    # Plot calculated vs convolved with slit
    # Specair units: mW/cm2/sr/µm
    w, I = np.loadtxt(getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt")).T
    s = calculated_spectrum(
        w, I, conditions={"Tvib": 3000, "Trot": 1200}, Iunit="mW/cm2/sr/µm"
    )

    if plot:
        fig = plt.figure(fig_prefix + "SPECAIR convoluted vs SPECAIR non convoluted")
        s.plot("radiance_noslit", nfig=fig.number, label="calc")
    slit_nm = 0.1
    s.apply_slit(slit_nm, shape="triangular", norm_by="area")
    if plot:
        s.plot(
            "radiance",
            nfig=fig.number,
            label="slit {0}nm".format(slit_nm),
            color="r",
            lw=2,
        )
        plt.legend()

    # Compare with Specair slit function
    s.apply_slit(slit_nm, norm_by="max")
    # Note unit conversion from mW/cm2/sr/um*nm is properly done!
    if plot:
        fig = plt.figure(fig_prefix + "convoluted RADIS vs convoluted SPECAIR")
        s.plot(
            "radiance",
            nfig=fig.number,
            label="slit {0}nm, norm with max".format(slit_nm),
            color="r",
            lw=2,
            Iunit="mW/cm2/sr",
        )
    # ... Get initial spectrum convoluted with Specair
    ws, Is = np.loadtxt(
        getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000_slit0.1.txt")
    ).T  # Specair units: mW/cm2/sr
    if plot:
        plt.plot(ws, Is, "k", label="normalized in Specair (sides not cropped)")

    # ... Test output is the same
    wu, Iu = s.get("radiance", Iunit="mW/cm2/sr")

    As = nantrapz(Is, ws)  # Todo one day: replace with get_power() function
    Au = nantrapz(Iu, wu)
    if verbose:
        print(
            (
                "Integrals should match: {0:.2f} vs {1:.2f} ({2:.2f}% error)".format(
                    As, Au, 100 * abs(As - Au) / As
                )
            )
        )
    assert np.isclose(As, Au, rtol=1e-2)

    # Test resampling
    s.resample(linspace(376, 380.6, 3000), unit="nm")
    s.apply_slit(slit_nm, norm_by="max")
    if plot:
        s.plot(
            "radiance",
            Iunit="mW/cm2/sr",
            lw=3,
            nfig=fig.number,
            color="b",
            zorder=-3,
            label="Resampled",
        )
        plt.legend()

    if verbose:
        print("\n>>>Testing spectrum slit matches Specair: OK")

    return True


@pytest.mark.fast
def test_normalisation_mode(plot=True, close_plots=True, verbose=True, *args, **kwargs):
    """Test norm_by = 'area' vs norm_by = 'max'"""

    from radis.test.utils import getTestFile

    _clean(plot, close_plots)

    # %% Compare spectra convolved with area=1 and max=1
    # Slit in nm
    # Spectrum in nm
    # Specair units: mW/cm2/sr/µm
    w, I = np.loadtxt(getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt")).T
    s = calculated_spectrum(
        w, I, conditions={"Tvib": 3000, "Trot": 1200}, Iunit="mW/cm2/sr/µm"
    )

    FWHM = 2

    s.apply_slit(FWHM, norm_by="area")  # spectrum convolved with area=1
    w_area, I_area = s.get("radiance")
    if plot:
        fig = plt.figure(fig_prefix + "Spectrum in nm + slit in nm")
        fig.clear()
        ax = fig.gca()
        s.plot(nfig=fig.number, wunit="nm", label="norm_by: area", lw=3)
    s.apply_slit(FWHM, norm_by="max")  # spectrum convolved with max=1
    w_max, I_max = s.get("radiance", wunit="nm")
    if plot:
        ax.plot(w_max, I_max / FWHM, "r", label="(norm_by:max)/FWHM")
        ax.legend(loc="best")
    assert np.allclose(I_area, I_max / FWHM, equal_nan=True)
    if verbose:
        print("equivalence of normalisation mode for spectrum in 'nm': OK")

    # %% Compare spectra convolved with area=1 and max=1
    # Slit in nm
    # Spectrum in cm-1

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update()
    # spectrum convolved with area=1
    s.apply_slit(FWHM, norm_by="area", plot_slit=plot)
    w_area, I_area = s.get("radiance")
    if plot:
        fig = plt.figure(fig_prefix + "Spectrum in cm-1 + slit in nm")
        fig.clear()
        ax = fig.gca()
        s.plot(nfig=fig.number, wunit="nm", label="norm_by: area", lw=3)
    # spectrum convolved with max=1
    s.apply_slit(FWHM, norm_by="max", plot_slit=plot)
    w_max, I_max = s.get("radiance", wunit="nm")
    if plot:
        ax.plot(w_max, I_max / FWHM, "r", label="(norm_by:max)/FWHM")
        ax.legend(loc="best")
    assert np.allclose(I_area, I_max / FWHM, equal_nan=True)
    if verbose:
        print("equivalence of normalisation mode for spectrum in 'cm-1': {0}: OK")
    assert is_homogeneous(s.units["radiance"], "mW/cm2/sr")
    if verbose:
        print(
            (
                "radiance unit ({0}) is homogeneous to 'mW/cm2/sr': OK".format(
                    s.units["radiance"]
                )
            )
        )

    return True


@pytest.mark.fast
def test_slit_energy_conservation(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """Convoluted and non convoluted quantities should have the same area
    (difference arises from side effects if the initial spectrum is not 0 on
    the sides"""

    from radis.test.utils import getTestFile

    _clean(plot, close_plots)

    if verbose:
        print("\n>>> _test_slit_energy_conservation\n")

    s = calculated_spectrum(
        *np.loadtxt(getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000_slit0.1.txt")).T,
        wunit="nm",
        Iunit="mW/cm2/sr/nm"
    )  # arbitrary)

    P = s.get_power(unit="mW/cm2/sr")
    s.apply_slit(0.5, norm_by="area")
    w, I = s.get("radiance", wunit="nm", Iunit="mW/cm2/sr/nm")
    Pc = abs(np.trapz(I[~np.isnan(I)], x=w[~np.isnan(I)]))  # mW/cm2/sr

    b = np.isclose(P, Pc, 3e-2)

    if plot:
        fig = plt.figure(fig_prefix + "energy conservation during resampling")
        s.plot(nfig=fig.number, label="{0:.1f} mW/cm2/sr".format(P))
        s.plot("radiance_noslit", nfig=fig.number, label="{0:.1f} mW/cm2/sr".format(Pc))
        plt.title("Energy conservation: {0}".format(b))
        plt.legend()
        plt.tight_layout()

    assert np.isclose(P, Pc, 3e-2)

    return True


# Function used to test Slit dispersion


def dirac(w0, width=20, wstep=0.009):
    """Return a Dirac in w0. Plus some space on the side"""

    w_side_p = np.arange(w0, w0 + width, wstep)
    w_side_n = np.arange(w0, w0 - width, -wstep)
    w = np.hstack((w_side_n[:0:-1], w_side_p))
    I = np.zeros_like(w)
    I[len(I) // 2] = 1 / wstep

    return w, I


def linear_dispersion(w, f=750, phi=-6, m=1, gr=300):
    """dlambda / dx
    Default values correspond to Acton 750i

    Parameters
    ----------

    f: focal length (mm)
         default 750 (SpectraPro 2750i)

    phi: angle in degrees (°)
        default 9

    m: order of dispersion
        default 1

    gr: grooves spacing (gr/mm)
        default 300
    """
    # correct units:
    phi *= 2 * pi / 360
    d = 1e-3 / gr
    disp = w / (2 * f) * (tan(phi) + sqrt((2 * d / m / (w * 1e-9) * cos(phi)) ** 2 - 1))
    return disp  # to nm/mm


@pytest.mark.fast
def test_linear_dispersion_effect(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """A test case to show the effect of wavelength dispersion (cf spectrometer
    reciprocal function) on the slit function

    Test succeeds if a :py:data:`~radis.misc.warning.SlitDispersionWarning`
    is correctly triggered and if the FWHM is properly dilated.
    """

    from radis.test.utils import getTestFile

    _clean(plot, close_plots)

    w_slit, I_slit = import_experimental_slit(getTestFile("slitfunction.txt"))

    if plot:
        plt.figure(fig_prefix + "Linear dispersion effect")
        plt.plot(
            w_slit,
            I_slit,
            "--k",
            label="Exp: FWHM @{0}nm: {1:.3f} nm".format(
                632.6, get_effective_FWHM(w_slit, I_slit)
            ),
        )

    from radis.misc.warning import SlitDispersionWarning

    with pytest.warns(
        SlitDispersionWarning
    ):  # expect a "large slit dispersion" warning

        # Test how slit function FWHM scales with linear_dispersion
        for w0, FWHM in zip([380, 1000, 4200, 5500], [0.393, 0.385, 0.280, 0.187]):
            w, I = dirac(w0)

            wc, Ic = offset_dilate_slit_function(
                w_slit,
                I_slit,
                w,
                slit_dispersion=linear_dispersion,
                threshold=0.01,
                verbose=False,
            )
            assert np.isclose(FWHM, get_effective_FWHM(wc, Ic), atol=0.001)
            if plot:
                plt.plot(
                    wc,
                    Ic,
                    label="FWHM @{0:.2f} nm: {1:.3f} nm".format(
                        w0, get_effective_FWHM(wc, Ic)
                    ),
                )

    if plot:
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Dirac $x$ slit function")
        plt.legend(loc="best", prop={"size": 15})
        fix_style()

    return True


@pytest.mark.fast
def test_cut_slices(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    """A test case to verify that _cut_slices does cut the spectrum into slices

    Test fails if a :py:data:`~radis.misc.warning.SlitDispersionWarning`
    is  triggered and succeeds if the spectrum is sliced into 8 segments.
    """
    from radis.misc.warning import SlitDispersionWarning

    _clean(plot, close_plots)

    threshold = 0.01
    w = np.arange(4000, 4400, 0.01)
    w_slit = np.arange(4198, 4202, 0.1)

    slices = _cut_slices(w, w_slit, linear_dispersion, threshold)

    for sl in slices:
        try:
            offset_dilate_slit_function(
                w_slit, np.ones_like(w_slit), w[sl], linear_dispersion, threshold, True
            )
        except SlitDispersionWarning:
            return False
        if plot:
            plt.plot(
                w,
                sl,
                label="Slice {0:.2f}-{1:.2f} nm , slit dispersion ratio: {2:.3f}".format(
                    w[sl][0],
                    w[sl][-1],
                    linear_dispersion(w[sl][-1]) / linear_dispersion(w[sl][0]),
                ),
            )
    if plot:
        plt.title("Cut slices, threshold (boundaries removed) = {0}".format(threshold))
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Slices (Boolean)")
        plt.legend(loc="best", prop={"size": 15})

    assert len(slices) == 8
    slices = _cut_slices(w[::-1], w_slit, linear_dispersion, threshold)
    assert len(slices) == 8

    return True


@pytest.mark.fast
def test_auto_correct_dispersion(
    verbose=True, plot=True, close_plots=True, *args, **kwargs
):
    """A test case to show the effect of wavelength dispersion (cf spectrometer
    reciprocal function) on the slit function

    Test succeeds if there is a difference between convoluted spectra (with and without slit dispersion)
    """

    from radis.test.utils import getTestFile

    _clean(plot, close_plots)

    slit_measured_632nm = getTestFile("slitfunction.txt")

    w, I = np.loadtxt(getTestFile("calc_N2C_spectrum_Trot1200_Tvib3000.txt")).T
    s = calculated_spectrum(
        w, I, conditions={"Tvib": 3000, "Trot": 1200}, Iunit="mW/cm2/sr/µm"
    )
    s2 = s.copy()

    def slit_dispersion(w):
        return linear_dispersion(w, f=750, phi=-6, m=1, gr=2400)

    s.apply_slit(slit_measured_632nm)

    if plot:
        w_slit_632, I_slit_632 = import_experimental_slit(
            getTestFile("slitfunction.txt")
        )
        w_full_range = np.linspace(w.min(), w_slit_632.max())
        plt.figure(
            "Spectrometer Dispersion (f={0}mm, phi={1}°, gr={2}".format(750, -0.6, 2400)
        )
        plt.plot(w_full_range, slit_dispersion(w_full_range))
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Reciprocal Linear Dispersion")

        # Compare 2 spectra
        s.plot(nfig="Linear dispersion effect", color="r", label="not corrected")

    s2 = s2.apply_slit(
        slit_measured_632nm, slit_dispersion=slit_dispersion, inplace=False
    )

    if plot:
        s2.plot(nfig="same", color="k", label="corrected")
        plt.legend()
        # Plot different slits:
        s2.plot_slit()
    assert np.isclose(
        s.take("radiance").max() / (s2.take("radiance").max()), 1.183, atol=0.001
    )
    return True


@pytest.mark.fast
def test_resampling(rtol=1e-2, verbose=True, plot=True, warnings=True, *args, **kwargs):
    """Test what happens when a spectrum in nm or cm-1, is convolved
    with a slit function in nm. In particular, slit function is generated
    in the spectrum unit, and spectrum is resampled if not evenly spaced"""

    if verbose:
        printm("Test auto resampling")

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    plCO = SpectrumFactory(
        wavenum_min=2230,
        wavenum_max=2260,
        mole_fraction=0.02,
        path_length=100,  # cm
        truncation=10,  # cm^-1
        wstep=0.02,
        isotope=[1, 2, 3],
        verbose=verbose,
    )
    plCO.warnings["MissingSelfBroadeningWarning"] = "ignore"
    plCO.load_databank("HITRAN-CO-TEST")
    sCO = plCO.eq_spectrum(Tgas=300)

    w_nm, T_nm = sCO.get("transmittance_noslit", wunit="nm")
    w_nm, I_nm = sCO.get("radiance_noslit", wunit="nm", Iunit="mW/cm2/sr/nm")
    sCO_nm = transmittance_spectrum(
        w_nm, T_nm, wunit="nm"
    )  # a new spectrum stored in nm
    # sCO_nm = theoretical_spectrum(w_nm, I_nm, wunit='nm', Iunit='mW/cm2/sr/nm') #  a new spectrum stored in nm

    if plot:
        fig = plt.figure(fig_prefix + "auto-resampling")
        sCO.plot(
            "transmittance_noslit",
            wunit="cm-1",
            nfig=fig.number,
            marker="o",
            color="k",
            lw=3,
            ms=10,
            label="(stored in cm-1)",
        )
        plt.title("No slit function")
        sCO_nm.plot(
            "transmittance_noslit",
            wunit="cm-1",
            nfig=fig.number,
            marker="o",
            color="r",
            label="(stored in nm)",
        )
        #            plt.xlim((2246.58, 2247.52))
        #            plt.ylim((0.87, 1.01))
        plt.legend()

    slit_function = 0.8
    slit_unit = "cm-1"
    sCO.apply_slit(slit_function, unit=slit_unit)
    sCO_nm.apply_slit(slit_function, unit=slit_unit)

    if plot:
        fig = plt.figure(fig_prefix + "auto-resampling (after convolution)")
        sCO.plot(
            "transmittance",
            wunit="cm-1",
            nfig=fig.number,
            marker="o",
            color="k",
            lw=3,
            ms=10,
            label="(stored in cm-1)",
        )
        plt.title("Slit function: {0} {1}".format(slit_function, slit_unit))
        sCO_nm.plot(
            "transmittance",
            wunit="cm-1",
            nfig=fig.number,
            marker="o",
            color="r",
            label="(stored in nm)",
        )

        #            plt.xlim((2246.58, 2247.52))
        #            plt.ylim((0.87, 1.01))
        plt.legend()

    w_conv, T_conv = sCO.get("transmittance", wunit="cm-1")
    w_nm_conv, T_nm_conv = sCO_nm.get("transmittance", wunit="cm-1")

    error = abs(
        (nantrapz(1 - T_conv, w_conv) - nantrapz(1 - T_nm_conv, w_nm_conv))
        / nantrapz(1 - T_nm_conv, w_nm_conv)
    )

    if verbose:
        printm("\n>>> _test_resampling\n")
    if verbose:
        printm(
            "Error between 2 spectra ({0:.2f}%) < {1:.2f}%: {2}".format(
                error * 100, rtol * 100, bool(error < rtol)
            )
        )
    assert bool(error < rtol)


def _run_testcases(plot=True, close_plots=False, verbose=True, *args, **kwargs):

    # Validation
    test_against_specair_convolution(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )

    # Different modes
    test_normalisation_mode(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )

    # Resampling
    test_slit_energy_conservation(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )

    # Slit dispersion
    test_linear_dispersion_effect(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )
    test_cut_slices(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )
    test_auto_correct_dispersion(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )

    # Different shapes
    test_all_slit_shapes(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )

    # Units
    test_slit_unit_conversions_spectrum_in_cm(
        verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs
    )
    test_slit_unit_conversions_spectrum_in_nm(
        verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs
    )
    test_convoluted_quantities_units(*args, **kwargs)

    test_resampling(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    print(("Testing slit.py: ", _run_testcases(plot=True)))
