#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CHIANTI Atomic Database Reader for RADIS
=========================================

Minimal reader for CHIANTI ``.wgfa`` (weighted oscillator strengths,
gf-values and radiative rates) format files.  Converts CHIANTI atomic
line data into a RADIS-compatible :class:`~pandas.DataFrame`.

The reader targets the ``.wgfa`` fixed-width format described in the
`CHIANTI database documentation <https://www.chiantidatabase.org/>`__.

References
----------
Dere, K.P. et al., 1997, A&AS, 125, 149 – CHIANTI atomic database.

See Also
--------
:py:func:`radis.io.hitran.fetch_hitran`,
:py:func:`radis.io.query.fetch_astroquery`
"""

import logging
import os

import pandas as pd
log = logging.getLogger(__name__)


def load_chianti(
    datapath=None,
    ion="fe_12",
    wmin=0.0,    # nm
    wmax=1e5,    # nm
):
    """Read CHIANTI atomic line data and return a RADIS-compatible DataFrame.

    Parameters
    ----------
    datapath : str, optional
        Root directory of the CHIANTI database tree.  The expected layout is
        ``<datapath>/<element>/<ion>/<ion>.wgfa``, e.g.
        ``chianti/fe/fe_12/fe_12.wgfa``.
        If *None*, defaults to the sample database shipped with RADIS at
        ``radis/db/chianti``.
    ion : str
        Ion identifier in CHIANTI notation, e.g. ``'fe_12'`` for Fe XII.
    wmin, wmax : float
        Wavelength filtering window **in nanometres (nm)**, consistent
        with the RADIS-wide convention (``fetch_nist``, ``fetch_astroquery``
        etc. all use nm).  Only transitions whose wavelength falls in
        ``[wmin, wmax]`` nm are returned.
        Defaults select all lines (0 – 100 000 nm).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:

        - ``wav`` : vacuum wavelength in **nm** (converted from Å)
        - ``gf``  : weighted oscillator strength *gf*
        - ``A``   : Einstein A coefficient (s⁻¹)
        - ``upper`` / ``lower`` : upper and lower level indices

    Raises
    ------
    FileNotFoundError
        If the ``.wgfa`` file cannot be located under *datapath*.

    Examples
    --------
    >>> from radis.io.chianti import load_chianti
    >>> df = load_chianti(ion='fe_12', wmin=19.0, wmax=20.0)
    >>> print(df[['wav', 'gf', 'A']].head())

    Notes
    -----
    CHIANTI ``.wgfa`` files use a fixed-width format.  The file is terminated
    by a sentinel line whose level indices are ``-1``.

    See https://www.chiantidatabase.org/ for the full database.
    """

    # --- resolve datapath ------------------------------------------------- #
    if datapath is None:
        import radis

        datapath = os.path.join(os.path.dirname(radis.__file__), "db", "chianti")

    # --- construct file path ---------------------------------------------- #
    ion_lower = ion.lower()
    element = ion_lower.split("_")[0]  # e.g. 'fe' from 'fe_12'
    wgfa_path = os.path.join(datapath, element, ion_lower, f"{ion_lower}.wgfa")

    if not os.path.exists(wgfa_path):
        raise FileNotFoundError(
            f"CHIANTI .wgfa file not found: {wgfa_path}\n"
            f"Expected layout: <datapath>/{element}/{ion_lower}/{ion_lower}.wgfa\n"
            f"Download the full database from https://www.chiantidatabase.org/"
        )

    # --- read fixed-width file -------------------------------------------- #
    # CHIANTI .wgfa columns (0-indexed character positions):
    #   col 0-5  : lower level index
    #   col 6-11 : upper level index
    #   col 12-25: wavelength (Å)
    #   col 26-38: gf value
    #   col 39-51: A value (s⁻¹)
    # We use pandas read_fwf with automatic width detection which works
    # well for the scientific-notation columns in .wgfa files.
    df = pd.read_fwf(
        wgfa_path,
        header=None,
        names=["lower", "upper", "wavelength_AA", "gf", "A", "_flag1", "_flag2"],
        comment="%",
    )

    # --- clean up --------------------------------------------------------- #
    # Drop the sentinel line (level indices == -1) and any non-numeric rows
    df = df[df["lower"] != -1].copy()

    for col in ("wavelength_AA", "gf", "A"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df.dropna(subset=["wavelength_AA", "gf", "A"], inplace=True)

    # Remove negative / zero wavelengths (CHIANTI uses negative λ for
    # theoretical-only lines in some ions; skip them here).
    df = df[df["wavelength_AA"] > 0]

    # --- filter by wavelength (nm → Å internally) ------------------------ #
    # wmin/wmax are in nm — convert to Å for internal filtering
    wmin_AA = wmin * 10.0
    wmax_AA = wmax * 10.0
    mask = (df["wavelength_AA"] >= wmin_AA) & (df["wavelength_AA"] <= wmax_AA)
    df = df.loc[mask].reset_index(drop=True)

    # --- build RADIS-compatible output ------------------------------------ #
    result = pd.DataFrame(
        {
            "wav": df["wavelength_AA"] * 0.1,  # Å → nm
            "gf": df["gf"],
            "A": df["A"],
            "upper": df["upper"].astype(int),
            "lower": df["lower"].astype(int),
        }
    )

    n = len(result)
    log.info("Loaded %d transitions from CHIANTI %s", n, ion)
    if n > 0:
        log.info(
            "  Wavelength range: %.2f – %.2f nm",
            result["wav"].min(),
            result["wav"].max(),
        )

    return result
