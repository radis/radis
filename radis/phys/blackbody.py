# -*- coding: utf-8 -*-
"""

Notes
-----

Planck functions:
- planck: planck radiation with wavelength input
- planck_wn: planck radiation with wavenumber input
- sPlanck: a RADIS :class:`~radis.spectrum.spectrum.Spectrum` blackbody object

Example
-------

Generate Earth blackbody::

    s = sPlanck(wavelength_min=3000, wavelength_max=50000,
                T=288, eps=1)
    s.plot()

-------------------------------------------------------------------------------


"""

from numpy import arange, exp, inf, ones_like, zeros_like

from radis.phys.air import air2vacuum
from radis.phys.constants import c, c_CGS, h, h_CGS, k_b, k_b_CGS
from radis.phys.units import Unit as Q_
from radis.phys.units import conv2


def planck(lmbda, T, eps=1, unit="mW/sr/cm2/nm"):
    r"""Planck function for blackbody radiation.

    .. math::
        \epsilon \frac{2h c^2}{{\lambda}^5} \frac{1}{\operatorname{exp}\left(\frac{h c}{\lambda k T}\right)-1}

    Parameters
    ----------
    λ: np.array   (nm)
       wavelength
    T: float    (K)
        equilibrium temperature
    eps: grey-body emissivity
        default 1
    unit: output unit
        default 'mW/sr/cm2/nm'

    Returns
    -------
    np.array :  (mW.sr-1.cm-2/nm)
        equilibrium radiance

    See Also
    --------
    :py:func:`~radis.tools.blackbody.sPlanck`, :py:func:`~radis.phys.blackbody.planck_wn`
    """
    k = k_b
    lbd = lmbda * 1e-9
    iplanck = (
        eps * (2 * h * c ** 2 / lbd ** 5) * (1 / (exp(h * c / (lbd * k * T)) - 1))
    )  # S.I  (W.sr-1.m-3)
    iplanck *= 1e-10  # W.sr-1.m-3 >>> mW.sr-1.cm-2.nm-1

    if Q_(unit) != Q_("mW/sr/cm2/nm"):
        iplanck = conv2(iplanck, "mW/sr/cm2/nm", unit)

    return iplanck


def planck_wn(wavenum, T, eps=1, unit="mW/sr/cm2/cm-1"):
    r"""Planck function for blackbody radiation, wavenumber version.

    .. math::
        \epsilon 2h c^2 {\nu}^3 \frac{1}{\operatorname{exp}\left(\frac{h c \nu}{k T}\right)-1}

    Parameters
    ----------
    wavenum: np.array   (cm-1)
       wavenumber
    T: float    (K)
        equilibrium temperature
    eps: grey-body emissivity
        default 1
    unit: str
        output unit. Default 'mW/sr/cm2/cm-1'


    Returns
    -------
    np.array :  default (mW/sr/cm2/cm-1)
        equilibrium radiance

    See Also
    --------
    :py:func:`~radis.tools.blackbody.sPlanck`, :py:func:`~radis.phys.blackbody.planck`
    """
    k = k_b_CGS
    h = h_CGS
    c = c_CGS

    iplanck = (
        eps
        * (2 * h * c ** 2 * wavenum ** 3)
        * (1 / (exp(h * c * wavenum / (k * T)) - 1))
    )
    # iplanck in erg/s/sr/cm2/cm-1
    iplanck *= 1e-4  # erg/s/sr/cm2/cm-1 > mW/sr/cm^2/cm-1

    if Q_(unit) != Q_("mW/sr/cm2/cm-1"):
        iplanck = conv2(iplanck, "mW/sr/cm2/cm-1", unit)

    return iplanck


# %% Predefined Spectra objects


def sPlanck(
    wavenum_min=None,
    wavenum_max=None,
    wavelength_min=None,
    wavelength_max=None,
    T=None,
    eps=1,
    wstep=0.01,
    medium="air",
    **kwargs
):
    r"""Return a RADIS :py:class:`~radis.spectrum.spectrum.Spectrum` object with blackbody radiation.

    It's easier to plug in a :py:func:`~radis.los.slabs.SerialSlabs` line-of-sight than the Planck
    radiance calculated by :py:func:`~radis.phys.blackbody.planck`.
    And you don't need to worry about units as they are handled internally.

    See :py:class:`~radis.spectrum.spectrum.Spectrum` documentation for more information

    Parameters
    ----------
    wavenum_min / wavenum_max: ():math:`cm^{-1}`)
        minimum / maximum wavenumber to be processed in :math:`cm^{-1}`.
    wavelength_min / wavelength_max: (:math:`nm`)
        minimum / maximum wavelength to be processed in :math:`nm`.
    T: float (K)
        blackbody temperature
    eps: float [0-1]
        blackbody emissivity. Default ``1``

    Other Parameters
    ----------------
    wstep: float (cm-1 or nm)
        wavespace step for calculation
    **kwargs: other keyword inputs
        all are forwarded to spectrum conditions. For instance you can add
        a 'path_length=1' after all the other arguments

    Examples
    --------
    Generate Earth blackbody::

        s = sPlanck(wavelength_min=3000, wavelength_max=50000,
                    T=288, eps=1)
        s.plot()

    Examples using sPlanck :

    .. minigallery:: radis.sPlanck

    References
    ----------
    In wavelength:

    .. math::
        \epsilon \frac{2h c^2}{{\lambda}^5} \frac{1}{\operatorname{exp}\left(\frac{h c}{\lambda k T}\right)-1}

    In wavenumber:

    .. math::
        \epsilon 2h c^2 {\nu}^3 \frac{1}{\operatorname{exp}\left(\frac{h c \nu}{k T}\right)-1}

    See Also
    --------
    :py:func:`~radis.phys.blackbody.planck`, :py:func:`~radis.phys.blackbody.planck_wn`

    """
    from radis.spectrum.spectrum import Spectrum

    # Check inputs
    if (wavelength_min is not None or wavelength_max is not None) and (
        wavenum_min is not None or wavenum_max is not None
    ):
        raise ValueError("You cannot give both wavelengths and wavenumbers")

    if wavenum_min is not None and wavenum_max is not None:
        assert wavenum_min < wavenum_max
        waveunit = "cm-1"
    else:
        assert wavelength_min < wavelength_max
        if medium == "air":
            waveunit = "nm"
        elif medium == "vacuum":
            waveunit = "nm_vac"
        else:
            raise ValueError(medium)

    if T is None:
        raise ValueError("T must be defined")

    if not (eps >= 0 and eps <= 1):
        raise ValueError("Emissivity must be in [0-1]")

    # Test range is correct:
    if waveunit == "cm-1":
        # generate the vector of wavenumbers (shape M)
        w = arange(wavenum_min, wavenum_max + wstep, wstep)
        Iunit = "mW/sr/cm2/cm-1"
        I = planck_wn(w, T, eps=eps, unit=Iunit)
    else:
        # generate the vector of lengths (shape M)
        w = arange(wavelength_min, wavelength_max + wstep, wstep)
        Iunit = "mW/sr/cm2/nm"
        if waveunit == "nm_vac":
            w_vac = w
        elif waveunit == "nm":
            w_vac = air2vacuum(w)
        # calculate planck with wavelengths in vacuum
        I = planck(w_vac, T, eps=eps, unit=Iunit)

    conditions = {"wstep": wstep}
    # add all extra parameters in conditions (ex: path_length)
    conditions.update(**kwargs)

    return Spectrum(
        quantities={
            "radiance_noslit": (w, I),
            "transmittance_noslit": (w, zeros_like(w)),
            "absorbance": (w, ones_like(w) * inf),
        },
        conditions=conditions,
        units={"radiance_noslit": Iunit, "transmittance_noslit": "", "absorbance": ""},
        cond_units={"wstep": waveunit},
        wunit=waveunit,
        name="Planck {0}K, eps={1:.2g}".format(T, eps),
    )


if __name__ == "__main__":
    from radis.test.phys.test_blackbody import _run_testcases

    _run_testcases(plot=True)
