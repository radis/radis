# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 18:47:39 2017

@author: erwan

State-dependant and state-independant degeneracies for molecules 


"""
#TODO: Make it a JSON file


def gi(M, I):
    ''' State independant degeneracy 

    Parameters
    ----------
    
    M: int
        molecule id
        
    I: int
        isotope number

    References
    ----------

    Šimečková 2006: "Einstein A-coefficients and statistical weights for molecular
    absorption transitions in the HITRAN database". DOI 10.1016/j.jqsrt.2005.07.003

    '''

    _gi = {
        2:  # CO2
            {1: 1,   # 626
             2: 2,   # 636
             3: 1,   # 628
             4: 6,   # 627
         },
        5:   # CO
            {1: 1,   # 26
             2: 2,   # 36
             3: 1,
             4: 6,
         },
    }

    try:
        return _gi[M][I]
    except KeyError:
        raise NotImplementedError('undefined state-independant degeneracy for ' +
                                  'molecule[isotope]: {0}[{1}]'.format(M, I))


def gs(M, I):
    ''' State dependant degeneracy

    Parameters
    ----------
    
    M: int
        molecule id
        
    I: int
        isotope number

    References
    ----------

    Šimečková 2006: "Einstein A-coefficients and statistical weights for molecular
    absorption transitions in the HITRAN database". DOI 10.1016/j.jqsrt.2005.07.003

    CO2::

        For the 12C16O2, 13C16O2, 12C18O2 isotopologues, the statistical weights gs
        of the symmetric s and antisymmetric a rotational levels are 1 and 0, respectively
        (see Fig. 99 of Ref. [24], Fig. 17-6 of Ref. [26]). For the 16O12C18O, 16O12C17O,
        16O13C18O, 16O13C17O, 17O12C18O isotopologues, the statistical weights gs
        of all the rotational levels are 1"

    '''

    _gs = {
        2:  # CO2
            {1: (1,0),   # 626       (normally 1:0, but negative levels dont exist)
             2: (1,0),   # 636       (normally 1:0, but negative levels dont exist)
             3: 1,       # 628 
             4: 1,       # 627
             },
        5:   # CO
            {1: 1,   # 26
             2: 1,   # 36
             3: 1,
             4: 1,
             },
    }

    try:
        return _gs[M][I]
    except KeyError:
        raise NotImplementedError('undefined state-dependant degeneracy for ' +
                                  'molecule[isotope]: {0}[{1}]'.format(M, I))

