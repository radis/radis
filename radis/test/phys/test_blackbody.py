# -*- coding: utf-8 -*-
"""
Test Planck functions

-------------------------------------------------------------------------------


"""

import numpy as np
from numpy import pi

from radis.phys.blackbody import planck, planck_wn, sPlanck


def test_exceptions(verbose=True, *args, **kwargs):
    """Here we test that the code actually crashes properly!"""

    import pytest

    # with both wavelength and wavenum
    with pytest.raises(ValueError):  # expected behavior
        sPlanck(
            wavelength_min=300,
            wavelength_max=2000,
            wavenum_min=300,
            wavenum_max=600,
            T=300,
            eps=1,
        )

    # with no temperature
    with pytest.raises(ValueError):  # expected behavior
        sPlanck(wavelength_min=300, wavelength_max=2000, T=None, eps=1)

    # with a wrong epsilon
    with pytest.raises(ValueError):  # expected behavior
        sPlanck(wavelength_min=300, wavelength_max=2000, T=300, eps=10)


def test_planck_nm(verbose=True, plot=True, *args, **kwargs):
    """Test blackbody with Wien's law, Stefan's law and tabulated data
    of maximum

    Reference
    ---------

    Tabulated intensity (75.987...) was checked to match SpectraPlot and
    manually calculated values

    """

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()  # dont get stuck with Matplotlib if executing through pytest

    T = 2200
    eps = 0.36

    s = sPlanck(wavelength_min=300, wavelength_max=8000, T=T, eps=eps, wstep=0.1)

    w_nm = s.get_wavelength(medium="vacuum")
    w_cm = s.get_wavenumber()

    Iunit_per_nm = "W/sr/nm/m2"
    Iunit_per_cm = "W/sr/cm-1/m2"

    if plot:
        s.plot("radiance_noslit", Iunit=Iunit_per_nm)

    I_nm = s.get_radiance_noslit(Iunit=Iunit_per_nm)
    I_cm = s.get_radiance_noslit(Iunit=Iunit_per_cm)

    # Check Wien's law
    w_nm_max = w_nm[I_nm.argmax()]
    assert np.isclose(2898 / T * 1e3, w_nm_max, rtol=1e-4)
    if verbose:
        print("Maximum emission at T={0}K: {1:.2f}nm".format(T, w_nm_max))
        print(".. test_planck_nm.py: Wien's law: OK")

    # Check Stefan's law
    sigma = 5.67e-8  # Stefan-Boltzmann constant   (W/m2/K-4)
    P_theory = eps * sigma * T**4  # Stefan's law (W/m2)
    P_calc = pi * s.get_power("W/m2/sr")  # Lambert's cosine law
    assert np.isclose(P_theory, P_calc, rtol=1e-1)
    if verbose:
        print("Blackbody radiant intensity (Stefan law): {0:.2f} W/m2".format(P_theory))
        print("Blackbody radiant intensity (calculated): {0:.2f} W/m2".format(P_calc))
        print(".. test_planck_nm.py: Stefan's law: OK")

    # Check that max is correct
    # hardcoded for 2200 K, epsilon = 0.36 (~1.3 µm)
    assert I_nm.max() == 75.98736024707178

    # Test planck and planck_wn
    assert np.allclose(planck(w_nm, T, eps, unit=Iunit_per_nm), I_nm)
    assert np.allclose(
        planck_wn(w_cm, T, eps, unit=Iunit_per_cm), I_cm, rtol=1e-2
    )  # higher tolerance because of numerical error
    # during conversion Iunit_per_nm to Iunit_per_cm

    return True


def test_planck_cm(verbose=True, plot=True, *args, **kwargs):
    """Validate Planck calculation with wavenumber

    Notes
    -----

    Earth blackbody radiation [1]_:

    - 294 K
    - eps = 0.65
    - integrated over 2pi steradian


    References
    ----------

    .. [1] Simple Earth Climate Model `NASA Smidth 2010 <https://www.giss.nasa.gov/research/briefs/schmidt_05/>`__

    """

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()  # dont get stuck with Matplotlib if executing through pytest

    T = 287.2
    eps = 0.78

    s = sPlanck(wavenum_min=10, wavenum_max=3000, T=T, eps=eps, wstep=0.1)

    w_nm = s.get_wavelength()
    w_cm = s.get_wavenumber()

    I_nm = s.get_radiance_noslit(Iunit="mW/sr/m2/nm")
    I_cm = s.get_radiance_noslit(Iunit="mW/sr/m2/cm-1")
    I_cm *= 2 * pi  # mW/m2/cm-1

    # Check Wien's law
    w_nm_max = w_nm[I_nm.argmax()]
    assert np.isclose(2898 / T * 1e3, w_nm_max, rtol=1e-3)
    if verbose:
        print("Maximum emission at T={0}K: {1:.2f}nm".format(T, w_nm_max))
        print(".. test_planck_cm.py: Wien's law: OK")

    # Check Stefan's law
    sigma = 5.67e-8  # Stefan-Boltzmann constant   (W/m2/K-4)
    P_theory = eps * sigma * T**4  # Stefan's law (W/m2)
    P_calc = pi * s.get_power("W/m2/sr")  # Lambert's cosine law
    assert np.isclose(P_theory, P_calc, rtol=1e-3)
    if verbose:
        print(
            "Earth blackbody radiant intensity (Stefan law): {0:.2f} W/m2".format(
                P_theory
            )
        )
        print(
            "Earth blackbody radiant intensity (calculated): {0:.2f} W/m2".format(
                P_calc
            )
        )
        print(".. test_planck_cm.py: Stefan's law: OK")

    if plot:
        plt.figure("test_blackbody.py: test_planck_cm")
        plt.plot(w_cm, I_cm)
        plt.title("Earth blackbody ({0} K, eps={1})".format(T, eps))
        plt.ylabel("Radiance (mW/m2/cm-1)")
        plt.xlabel("Wavenumber (cm-1)")

    return True


def _run_testcases(plot=True, verbose=True, warnings=True, *args, **kwargs):
    """Test procedures

    Parameters
    ----------

    debug: boolean
        swamps the console namespace with local variables. Default ``False``

    """

    # Test all Spectrum methods
    # -------------------------
    test_exceptions()
    test_planck_nm(verbose=verbose, plot=plot, *args, **kwargs)
    test_planck_cm(verbose=verbose, plot=plot, *args, **kwargs)

    return True


if __name__ == "__main__":
    print("test_blackbody:", _run_testcases(verbose=True, plot=True))
