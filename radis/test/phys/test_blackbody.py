# -*- coding: utf-8 -*-
"""

"""

import numpy as np
from radis.phys.blackbody import sPlanck, planck, planck_wn

def test_planck(verbose=True, plot=False, *args, **kwargs):
    
    T = 2200
    eps = 0.36
    
    s = sPlanck(wavelength_min=300, wavelength_max=2000, T=T, eps=eps,
                wstep=0.1)
    
    w_nm = s.get_wavelength()
    w_cm = s.get_wavenumber()
    
    Iunit_per_nm = 'W/sr/nm/m2'
    Iunit_per_cm = 'W/sr/cm_1/m2'
    
    if plot:
        s.plot('radiance_noslit', Iunit=Iunit_per_nm)
    
    
    I_nm = s.get_radiance_noslit(Iunit=Iunit_per_nm)
    I_cm = s.get_radiance_noslit(Iunit=Iunit_per_cm)

    # Check Wien's law
    w_nm_max = w_nm[I_nm.argmax()]
    assert np.isclose(2898/T*1e3, w_nm_max, rtol=1e-4)
    if verbose: print("test_blackbody.py: Wien's law: OK")
    
    # Check that max is correct
    assert I_nm.max() == 75.987331723070341  # hardcoded for 2200 K, epsilon = 0.36 (~1.3 µm)
    
    # Test planck and planck_wn
    assert np.allclose(planck(w_nm, T, eps, unit=Iunit_per_nm), I_nm)
    assert np.allclose(planck_wn(w_cm, T, eps, unit=Iunit_per_cm), I_cm,
                       rtol=1e-2)  # higher tolerance because of numerical error 
                                   # during conversion Iunit_per_nm to Iunit_per_cm
    
    return True

if __name__ == '__main__':
    print('Test planck function:', test_planck(verbose=True, plot=True))
