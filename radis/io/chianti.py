import os
import pandas as pd
def read_chianti(
    datapath: str,
    ion: str = "fe_12",
    wmin: float = 100.0,
    wmax: float = 200.0,
) -> pd.DataFrame:
    """
    Read CHIANTI atomic data for a given ion.

    Parameters
    ----------
    datapath : str
        Path to the root of the CHIANTI database (e.g. ``~/data/chianti``).
    ion : str
        Ion name in CHIANTI folder notation, e.g. ``"fe_12"``.
    wmin, wmax : float
        Wavelength range (Angstrom) to keep.

    Returns
    -------
    pandas.DataFrame
        Atomic line list in RADIS-compatible format.
    """

    
    ion_dir = os.path.join(datapath, "fe", ion)

    elvlc_path = os.path.join(ion_dir, f"{ion}.elvlc")
    wgfa_path = os.path.join(ion_dir, f"{ion}.wgfa")

    # TODO:
    raise NotImplementedError(
        f"CHIANTI reader not implemented yet. ELVLC={elvlc_path}, WGFA={wgfa_path}"
    )
    Notes
    -----
    Initial implementation will focus on Fe XII (``fe_12``) in the EUV range.
    """
    # TODO:
    # 1) Build ion folder path: os.path.join(datapath, "fe", ion)
    # 2) Read <ion>.elvlc (levels) into a DataFrame
    # 3) Read <ion>.wgfa (transitions) into a DataFrame
    # 4) Merge levels + transitions into one line-list DataFrame
    # 5) Map CHIANTI columns -> RADIS atomic columns
    # 6) Filter on wavelength between wmin and wmax
    __all__ = ["read_chianti"]

