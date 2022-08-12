# -*- coding: utf-8 -*-
"""
Test spectrum creation / export
"""

import pytest

from radis import Spectrum, get_residual, test_spectrum


@pytest.mark.fast
def test_specutils_io(verbose=True, plot=False, *args, **kwargs):
    """Test Radis/specutils  interface in particular when dealing with wavelengths

    see https://github.com/radis/radis/pull/499
    """

    s = test_spectrum().take("transmittance_noslit")

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


if __name__ == "__main__":
    test_specutils_io(plot=True)
