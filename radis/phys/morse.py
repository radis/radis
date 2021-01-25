# -*- coding: utf-8 -*-
"""Created on Tue Jul 31 16:28:43 2018.

@author: erwan
"""


def morse_increment(Dzero, we, wexe=0, weye=0, weze=0, weae=0, webe=0, wece=0):
    """Get Morse potential energy increment correction for given molecule.

    Parameters
    ----------

    Dzero: cm-1
        dissociation energy

    we, wexe, etc.: cm-1
        spectroscopic parameters

    Notes
    -----

    Internal:

        Dzero = self.Ediss
        d_equil = (Dzero+0.5*we-0.25*wexe+0.125*weye
                   + 0.0625*weze+0.03125*weae+0.015625*webe
                   + 7.8125e-3*wece)

        v_inc = we**2 / (2*d_equil)

    References
    -----------

    Specair
    """

    d_equil = (
        Dzero
        + 0.5 * we
        - 0.25 * wexe
        + 0.125 * weye
        + 0.0625 * weze
        + 0.03125 * weae
        + 0.015625 * webe
        + 7.8125e-3 * wece
    )

    v_inc = we ** 2 / (2 * d_equil)

    return v_inc
