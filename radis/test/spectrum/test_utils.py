# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 02:04:14 2021

@author: erwan
"""

import os
import os.path
import time


def test_perf_profile(*args, **kwargs):
    """Test visual/interactive performance profile

    Examples
    --------
    ::
        s = calc_spectrum(...)
        s.print_perf_profile()

        # output >>
            spectrum_calculation      0.844s ████████████████
                check_line_databank              0.001s
                check_non_eq_param               0.070s █
                fetch_energy_5                   0.027s
                calc_weight_trans                0.016s
                reinitialize                     0.003s
                    copy_database                    0.000s
                    memory_usage_warning             0.003s
                    reset_population                 0.000s
                calc_noneq_population            0.066s █
                    part_function                    0.056s █
                    population                       0.010s
                scaled_non_eq_linestrength       0.004s
                    map_part_func                    0.001s
                    corrected_population_se          0.003s
                calc_emission_integral           0.008s
                applied_linestrength_cutoff      0.003s
                calc_lineshift                   0.001s
                calc_hwhm                        0.008s
                generate_wavenumber_arrays       0.002s
                calc_line_broadening             0.668s ████████████
                    precompute_DLM_lineshapes        0.020s
                    DLM_Initialized_vectors          0.000s
                    DLM_closest_matching_line        0.001s
                    DLM_Distribute_lines             0.002s
                    DLM_convolve                     0.304s █████
                calc_other_spectral_quan         0.005s
                generate_spectrum_obj            0.000s
    ::
        s.generate_perf_profile()

    See output of :py:meth:`~radis.spectrum.spectrum.Spectrum.generate_perf_profile`
    in https://github.com/radis/radis/pull/325

    See Also
    --------
    :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile`,
    :py:meth:`~radis.spectrum.spectrum.Spectrum.generate_perf_profile`
    """
    from radis.lbl.calc import calc_spectrum
    from radis.test.utils import setup_test_line_databases

    setup_test_line_databases()

    s = calc_spectrum(
        2000,
        2300,
        molecule="CO",
        Tvib=3000,
        Trot=2000,
        databank="HITRAN-CO-TEST",
        verbose=3,
    )

    # Test Console profiler
    s.print_perf_profile()

    # Test interactive profiler
    if os.path.exists("spectrum.prof"):
        os.remove("spectrum.prof")
    s.generate_perf_profile()
    time.sleep(
        5
    )  # note @dev: would be better to wait for the subprocess to end; however it won't even if the window is closed
    with open("spectrum.prof", "rb") as f:
        assert "calculation_time" in str(f.read())


def test_perf_profile_from_factory(*args, **kwargs):
    """ See :py:func:`radis.test.spectrum.test_utils.test_perf_profile`"""

    from radis import SpectrumFactory
    from radis.test.utils import setup_test_line_databases

    setup_test_line_databases()

    sf = SpectrumFactory(
        2000,
        2300,
        molecule="CO",
        verbose=3,
    )
    sf.load_databank("HITRAN-CO-TEST")
    sf.non_eq_spectrum(3000, 2000)

    # Test Console profiler
    sf.print_perf_profile()

    # Test interactive profiler
    if os.path.exists("spectrum.prof"):
        os.remove("spectrum.prof")
    sf.generate_perf_profile()
    time.sleep(
        5
    )  # note @dev: would be better to wait for the subprocess to end; however it won't even if the window is closed
    with open("spectrum.prof", "rb") as f:
        assert "calculation_time" in str(f.read())


if __name__ == "__main__":
    test_perf_profile()
    # test_perf_profile_from_factory()
