# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests::

    pytest      # (in command line, in project folder)

Run only fast tests (i.e: tests that a  'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------

"""

import numpy as np
import pytest

from radis.los.slabs import MergeSlabs, SerialSlabs


@pytest.mark.fast
def test_merge_slabs(
    verbose=True,
    plot=False,
    close_plots=True,
    warnings=True,
    debug=False,
    *args,
    **kwargs,
):
    """Merge 10 slabs with 1/10 of concentration, and compare the results.
    Ensure error is < 0.1%

    Note that we won't have exactly the same results as the broadening is not
    updated in MergeSlabs
    """

    import matplotlib.pyplot as plt

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    if plot:
        plt.ion()  # dont get stuck with Matplotlib if executing through pytest
        if close_plots:
            plt.close("all")

    for optically_thin in [True, False]:

        # Get Some spectra
        s1 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
        s2 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"), binary=True)
        s1.update("all")
        s2.update("all")

        # Merge 50 times s1 in the same slab
        s50 = [s1] * 50
        s1N = MergeSlabs(*s50, optically_thin=optically_thin, debug=debug)

        if plot:
            for k in ["radiance_noslit"]:  # , 'transmittance_noslit']:
                s2.plot(k, lw=3, label="1x[CO=0.5]")
                s1N.plot(k, nfig="same", label="50x[CO=0.01]")
                plt.legend()
                plt.title("Optically thin: {0}".format(optically_thin))
                plt.tight_layout()

        if verbose:
            print("test_merge_slabs")
            print(
                "... Compare 50x[CO=0.01] vs 1x[CON=0.5] (optically thin: {0})".format(
                    optically_thin
                )
            )
            print(
                "... Difference: {0:.2f}%".format(
                    abs(s1N.get_power() / s2.get_power() - 1) * 100
                )
            )
        assert np.isclose(s2.get_power(), s1N.get_power(), 1.5e-2)

    return True


@pytest.mark.fast
def test_equilibrium_condition():
    """See issue #370

    'thermal_equilibrium' is used to reocmpute some spectral arrays from others,
    it is important to make sure a line-of-sight doesn't assume thermal equilibrium wrongly
    """

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    s1 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s2 = s1.copy()

    s2.conditions["thermal_equilibrium"] = False
    assert s1.conditions["thermal_equilibrium"] == True
    assert s2.conditions["thermal_equilibrium"] == False

    # Test MergeSlabs
    assert (s1 // s1).is_at_equilibrium() == True
    assert (s2 // s2).is_at_equilibrium() == False
    assert (s1 // s2).is_at_equilibrium() == False
    assert (s2 // s1).is_at_equilibrium() == False

    # Test SerialSlabs
    assert (s1 < s1).is_at_equilibrium() == False
    assert (s2 < s2).is_at_equilibrium() == False
    assert (s1 < s2).is_at_equilibrium() == False
    assert (s2 < s1).is_at_equilibrium() == False


@pytest.mark.fast
def test_serial_slabs_transmittance(
    verbose=True, plot=False, warnings=True, debug=False, *args, **kwargs
):
    """Add some slabs in serie, ensure that overall transmittance decreases
    Also check that some quantities are propagated (like path length?)

    """

    import matplotlib.pyplot as plt

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    if plot:
        plt.ion()  # dont get stuck with Matplotlib if executing through pytest

    # Get Some spectra
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"), binary=True)
    #    s.update('all')
    assert "transmittance_noslit" not in s.get_vars()

    s_los = SerialSlabs(s, s, s)

    w, T = s_los.get("transmittance_noslit")
    # note that it was calculated despite transmittance not in the initial
    # spectrum

    s.update("transmittance_noslit")  # get from abscoeff
    w1, T1 = s.get("transmittance_noslit")

    # Check it looks alright
    assert (T <= T1).all()
    assert (T != T1).any()
    assert (w == w1).all()

    # because they were the same slabs, intensive conditions should have been propagated
    assert s.conditions["Tgas"] == s_los.conditions["Tgas"]
    # but extensive shouldnt:
    assert s.conditions["path_length"] * 3 == s_los.conditions["path_length"]

    if verbose:
        print("Tested Serialslabs transmittance decreases s1>s2: OK")

    return True


@pytest.mark.fast
def test_serial_slabs_radiance(
    verbose=True, plot=False, warnings=True, debug=False, *args, **kwargs
):
    """Add some slabs in serie under optically thin conditions, ensures
    that radiance is about the sum of all.
    """

    import matplotlib.pyplot as plt

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec

    if plot:
        plt.ion()  # dont get stuck with Matplotlib if executing through pytest

    # Get Some spectra
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update("radiance_noslit")

    s_los = SerialSlabs(s, s)

    w, I = s_los.get("radiance_noslit")
    w1, I1 = s.get("radiance_noslit")

    if plot:
        s.plot(nfig="test_serial_slabs_radiance", lw=2)
        plt.title("s>s for almost optically thin conditions")
        s_los.plot(nfig="same")

    # Check it looks alright
    assert np.allclose(I, 2 * I1, rtol=1e-3)

    if verbose:
        print(
            "Tested Serialslabs radiance is added when approximately optically thin: OK"
        )

    # Test operators
    s_serial = s > s
    assert s_serial == s_los

    return True


def _run_testcases(
    verbose=True,
    plot=True,
    close_plots=True,
    debug=False,
    warnings=True,
    *args,
    **kwargs,
):

    test_merge_slabs(
        verbose=verbose,
        plot=plot,
        close_plots=close_plots,
        debug=debug,
        warnings=warnings,
        *args,
        **kwargs,
    )

    test_serial_slabs_transmittance(
        verbose=verbose, plot=plot, debug=debug, warnings=warnings, *args, **kwargs
    )

    test_serial_slabs_radiance(
        verbose=verbose, plot=plot, debug=debug, warnings=warnings, *args, **kwargs
    )

    test_equilibrium_condition()

    return True


if __name__ == "__main__":
    print("test_slabs: ", _run_testcases(verbose=True))
