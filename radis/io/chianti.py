#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CHIANTI Atomic Database Reader for RADIS

Minimal reader for CHIANTI .wgfa format files.
Converts CHIANTI data to RADIS-compatible format.
"""

import pandas as pd
import numpy as np
import os


class DatabankNotFound(Exception):
    """Exception raised when database file is not found"""
    pass


def load_chianti(
    datapath: str = None,
    ion: str = "fe_12",
    wmin: float = 193.0,
    wmax: float = 197.0,
) -> pd.DataFrame:
    """
    Read CHIANTI atomic line data and return in RADIS format.
    
    Parameters
    ----------
    datapath : str, optional
        Root CHIANTI database directory (containing 'fe' folder with 'fe_12', etc.)
        If None, uses: '/Users/jaiswaljaishankar05/Downloads/radis/radis/db/chianti'
    ion : str
        Ion identifier (e.g., 'fe_12' for Fe XII)
    wmin, wmax : float
        Wavelength range in Angstrom
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: wav (nm), int, A, gpp
        
    Raises
    ------
    DatabankNotFound
        If CHIANTI files not found
    
    Examples
    --------
    >>> df = load_chianti(ion='fe_12', wmin=193.0, wmax=197.0)
    >>> print(df.head())
    
    Notes
    -----
    CHIANTI format: https://www.chiantidatabase.org/
    """
    
    # Default CHIANTI path (update this to your actual path!)
    if datapath is None:
         datapath = '/Users/adityapandey/Downloads/CHIANTI_11.0.2_database'
    
    # Construct file paths
    ion_lower = ion.lower()
    element = ion_lower.split('_')[0]  # 'fe' from 'fe_12'
    ion_dir = os.path.join(datapath, element, ion_lower)
    wgfa_path = os.path.join(ion_dir, f"{ion_lower}.wgfa")
    
    # Check file exists
    if not os.path.exists(wgfa_path):
        raise DatabankNotFound(
            f"CHIANTI file not found: {wgfa_path}\n"
            f"Expected structure: {datapath}/fe/fe_12/fe_12.wgfa\n"
            f"Download from: https://www.chiantidatabase.org/"
        )
    
    # Read WGFA file (fixed-width format)
    try:
        df = pd.read_fwf(
            wgfa_path,
            skiprows=12,  # Skip header lines
            header=None,
            colspecs=[
                (0, 6),      # upper level index
                (6, 12),     # lower level index
                (12, 20),    # wavelength (Angstrom)
                (20, 30),    # gf (oscillator strength)
                (30, 40),    # A (Einstein A coefficient, s^-1)
            ],
            names=['upper', 'lower', 'wavelength', 'gf', 'A']
        )
    except Exception as e:
        raise DatabankNotFound(f"Error reading {wgfa_path}: {e}")
    
    # Clean data - convert to numeric
    for col in ['wavelength', 'gf', 'A']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Remove NaN rows
    df = df.dropna(subset=['wavelength', 'gf', 'A'])
    
    # Create RADIS-compatible DataFrame
    radis_df = pd.DataFrame({
        'wav': df['wavelength'] * 0.1,  # Angstrom to nm
        'int': df['gf'],                 # oscillator strength
        'A': df['A'],                    # Einstein A coefficient
        'gpp': 1.0,                      # placeholder (lower level degeneracy)
        'upper': df['upper'],
        'lower': df['lower']
    })
    
    # Filter by wavelength range (in Angstrom)
    mask = (df['wavelength'] >= wmin) & (df['wavelength'] <= wmax)
    result = radis_df[mask].reset_index(drop=True)
    
    print(f"✅ Loaded {len(result)} transitions from CHIANTI {ion}")
    if len(result) > 0:
        print(f"   Wavelength range: {result['wav'].min():.2f} — {result['wav'].max():.2f} nm")
        print(f"   First 5 lines:")
        print(result[['wav', 'int', 'A']].head())
    
    return result