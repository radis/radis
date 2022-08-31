# -*- coding: utf-8 -*-
"""
Test RTE and rescale equations from :py:meth:`radis.spectrum.spectrum.Spectrum.update`
"""


# test_equations
from numpy import allclose, exp


def test_equations(*args, **kwargs):

    import radis
    from radis.phys.blackbody import planck, planck_wn

    s = radis.test_spectrum(verbose=False)
    s.apply_slit(0.5, "nm")

    for k in list(s._q.keys()):  # reset all quantities
        if k in ["wavespace", "wavelength", "wavenumber"]:
            pass
        elif k == "abscoeff":
            new_abscoeff = s._q["abscoeff"].copy()
        else:
            del s._q[k]

    # Recompute quantities ourselves
    new_absorbance = new_abscoeff * s.conditions["path_length"]
    new_transmittance_noslit = exp(-new_absorbance)
    new_emissivity_noslit = 1 - new_transmittance_noslit
    new_radiance_noslit = new_emissivity_noslit * planck_wn(
        s.get_wavenumber(),
        s.conditions["Tgas"],
        unit="W/sr/cm2/cm-1",
    )
    new_radiance_noslit_wavelength = new_emissivity_noslit * planck(
        1e7 / s.get_wavenumber(),
        s.conditions["Tgas"],
        unit="W/sr/cm2/nm",
    )

    # Compare with arrays recomputed automatically
    assert (s.get("absorbance")[1] == new_absorbance).all()
    assert (s.get("transmittance_noslit")[1] == new_transmittance_noslit).all()
    assert (s.get("emissivity_noslit")[1] == new_emissivity_noslit).all()
    assert allclose(
        s.get("radiance_noslit", Iunit="W/sr/cm2/cm-1")[1], new_radiance_noslit
    )
    assert allclose(
        s.get("radiance_noslit", Iunit="W/sr/cm2/nm")[1], new_radiance_noslit_wavelength
    )


if __name__ == "__main__":
    test_equations()
