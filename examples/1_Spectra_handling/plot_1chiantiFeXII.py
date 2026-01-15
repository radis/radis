"""
Compare Fe XII atomic spectrum from CHIANTI with NIST database.
Follows radis line list import, SpectrumFactory, plot_diff workflow.
Uses the new CHIANTI reader.
"""

from radis import SpectrumFactory, plot_diff
import pandas as pd
from radis import config  # for databank path

# === STEP 1: CHIANTI Reader (tumhara code) ===
def read_chianti(datapath, ion, wmin, wmax):
    """
    Read CHIANTI data for ion and wavelength range.
    Parameters:
        datapath: str - CHIANTI root directory
        ion: str - e.g. "fe_12"
        wmin, wmax: float - wavelength range in Angstrom
    Returns DataFrame with at least 'wave', 'upper', 'lower', 'A21' columns.
    """
    import os
    import pandas as pd
    
    ion_path = os.path.join(datapath, ion)
    el, stage = ion.split("_")
    
    # Find transition file (e.g., fe_12_xxx.elc)
    files = os.listdir(ion_path)
    elc_file = None
    for f in files:
        if f.endswith(".elc") and el in f:
            elc_file = os.path.join(ion_path, f)
            break
    
    if elc_file is None:
        raise FileNotFoundError(f"No .elc file found in {ion_path}")
    
    # Parse ELC file (simplified)
    df = pd.read_csv(elc_file, sep=r'\s+', comment='!', 
                     names=['wave', 'upper', 'lower', 'A21', 'loggf', 'intensity'])
    df = df[(df['wave'] >= wmin) & (df['wave'] <= wmax)]
    
    print(f"CHIANTI lines loaded: {len(df)} transitions")
    print(df.head())
    
    return df

# === MAIN EXAMPLE ===
if __name__ == "__main__":
    
    # Fe XII EUV range (Angstrom)
    wavelength_min = 193.0  # 19.3 nm
    wavelength_max = 197.0  # 19.7 nm
    
    # CHIANTI path (tumhare Desktop pe)
    chianti_root = "/Users/adityapandey/Desktop/CHIANTI_v10.0"
    
    print("=== Loading CHIANTI Fe XII lines ===")
    df_chianti = read_chianti(
        datapath=chianti_root,
        ion="fe_12",
        wmin=wavelength_min,
        wmax=wavelength_max,
    )
    
    print(f"CHIANTI head: {df_chianti.head()}")
    
    # === RADIS SpectrumFactory CHIANTI ===
    sf_chianti = SpectrumFactory(
        193.0, 197.0,           # wavelength range (Å)
        code='python',
        verbose=2,
        wavenum=False,
        isotopologues='',
    )
    
    # Load CHIANTI lines into RADIS (Å -> nm conversion)
    chianti_lines = {
        'wave': df_chianti['wave'].values / 10,     # Å -> nm
        'A': df_chianti['A21'].values,
        'upper': df_chianti['upper'].values,
        'lower': df_chianti['lower'].values,
        'label': ['CHIANTI'] * len(df_chianti)
    }
    
    sf_chianti.load_databank(chianti_lines, 'chianti')
    
    # Compute CHIANTI spectrum
    Tgas = 10000
    s_chianti = sf_chianti.eq_spectrum(Tgas=Tgas)
    print("CHIANTI spectrum computed!")
    
    # === NIST Spectrum (same range) ===
    print("=== Loading NIST Fe XII ===")
    sf_nist = SpectrumFactory(
        193.0, 197.0,           # same wavelength range
        code='python',
        verbose=2,
        wavenum=False,
    )
    sf_nist.fetch_databank('nist', species='Fe')
    
    s_nist = sf_nist.eq_spectrum(Tgas=Tgas)
    print("NIST spectrum computed!")
    
    # === COMPARISON PLOT ===
    print("Plotting CHIANTI vs NIST comparison...")
    plot_diff(s_chianti, s_nist, 'radiance_noslit', wunit='nm', 
              title='Fe XII: CHIANTI vs NIST (19.3-19.7 nm)')
    
    print("✅ CHIANTI Fe XII example complete!")
