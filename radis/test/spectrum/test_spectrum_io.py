# -*- coding: utf-8 -*-
"""
Test spectrum creation / export
"""

import pytest

from radis import Spectrum, get_residual, spectrum_test


@pytest.mark.fast
def test_specutils_io(verbose=True, plot=False, *args, **kwargs):
    """Test Radis/specutils  interface in particular when dealing with wavelengths

    see https://github.com/radis/radis/pull/499
    """

    s = spectrum_test().take("transmittance_noslit")

    processed_spec = Spectrum.from_array(
        *s.get("transmittance_noslit", wunit="nm"),
        "transmittance_noslit",
        wunit="nm",
        Iunit="",
    )

    with pytest.raises(ValueError) as err:

        spec_u = processed_spec.to_specutils(
            "transmittance_noslit", wunit="nm_air", Iunit="default"
        )
        assert (
            "specutil's Spectrum1D will only handle wavelengths as seen in vacuum"
            in str(err)
        )

    spec_u = processed_spec.to_specutils(
        "transmittance_noslit", wunit="nm_vac", Iunit="default"
    )

    if plot:

        import matplotlib.pyplot as plt

        s.plot("transmittance_noslit", wunit="nm_vac", lw=5)
        plt.plot(spec_u.wavelength.to("nm"), spec_u.flux, lw=3)
        plt.xlim((4591.91, 4597.57))

    assert spec_u.meta["waveunit"] == "nm_vac"  # was corrected properly
    # Reverse

    s2 = Spectrum.from_specutils(spec_u, "transmittance_noslit")

    if plot:

        s2.plot("transmittance_noslit", wunit="nm_vac", lw=1, nfig="same")

    assert get_residual(s, s2, "transmittance_noslit", ignore_nan=True) < 1e-9


@pytest.mark.fast
def test_reading_from_Matlab(verbose=True, plot=False, *args, **kwargs):
    """Test creating a Spectrum from a Matlab file with
    :py:meth:`~radis.spectrum.spectrum.Spectrum.from_mat`
    """

    from radis import Spectrum
    from radis.test.utils import getTestFile

    s = Spectrum.from_mat(
        getTestFile("trimmed_1857_VoigtCO_Minesi.mat"),
        "absorbance",
        wunit="cm-1",
        unit="",
        index=10,
    )

    import numpy as np

    assert np.isclose(s.get("absorbance")[0][100], 2011.5356487633755)
    assert np.isclose(s.get("absorbance")[1][100], 0.005560075444710532)


def test_xsc_io(plot=False, *args, **kwargs):
    """Test Spectrum created from manually downloaded hitran cross section"""

    from radis.phys.units import Unit as u
    from radis.test.utils import getTestFile

    datafile = "CH3COCH3_233.4_375.2_700.0-1780.0_13.xsc"

    s = Spectrum.from_xsc(getTestFile(datafile))

    assert s.c["Tgas"] == 233.4 * u("K")
    assert s.c["pressure"] == 375.2 * u("Torr")

    if plot:
        s.plot()


if __name__ == "__main__":
    test_specutils_io(plot=True)
    test_reading_from_Matlab(plot=True)
    test_xsc_io(plot=True)
