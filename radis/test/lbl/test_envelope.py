# -*- coding: utf-8 -*-
"""
Test :py:mod:`radis.lbl.envelope` — Gaussian and Lorentzian width envelope
pre-computation.
"""

import json

import numpy as np
import pytest


@pytest.mark.fast
def test_gaussian_envelope_basic():
    """compute_gaussian_envelope returns correct min/max."""
    from radis.lbl.envelope import compute_gaussian_envelope

    log_2vMm = np.array([-3.0, -2.5, -1.0, -2.0, -1.5], dtype=np.float32)
    g_min, g_max = compute_gaussian_envelope(log_2vMm)
    assert np.isclose(g_min, -3.0, atol=1e-6)
    assert np.isclose(g_max, -1.0, atol=1e-6)


@pytest.mark.fast
def test_gaussian_envelope_single_line():
    """Edge case: single line."""
    from radis.lbl.envelope import compute_gaussian_envelope

    log_2vMm = np.array([-2.0], dtype=np.float32)
    g_min, g_max = compute_gaussian_envelope(log_2vMm)
    assert np.isclose(g_min, -2.0, atol=1e-6)
    assert np.isclose(g_max, -2.0, atol=1e-6)


@pytest.mark.fast
def test_lorentzian_envelope_single_na():
    """When all lines have the same na, the envelope is trivial."""
    from radis.lbl.envelope import compute_lorentzian_envelope

    np.random.seed(99)
    N = 100
    na = np.full(N, 0.7, dtype=np.float32)
    gamma_arr = np.random.uniform(0.01, 0.1, size=(2, N)).astype(np.float32)

    envelope = compute_lorentzian_envelope(na, gamma_arr)
    assert len(envelope) == 2  # (min_params, max_params)

    for params in envelope:
        A, B, X = params
        assert len(A) == len(B)
        assert len(A) == len(X)
        assert isinstance(A, np.ndarray)


@pytest.mark.fast
def test_lorentzian_envelope_bounds():
    """Envelope must bound the actual per-line Lorentzian widths."""
    from radis.lbl.envelope import compute_lorentzian_envelope

    np.random.seed(42)
    N = 500
    na = np.random.uniform(0.3, 0.9, size=N).astype(np.float32)
    gamma_arr = np.random.uniform(0.01, 0.15, size=(2, N)).astype(np.float32)

    envelope = compute_lorentzian_envelope(na, gamma_arr)

    for T in [300, 1000, 3000, 5000]:
        for p in [0.01, 1.0, 10.0]:
            log_rT = np.log(296.0 / T)
            log_2p = np.log(2 * p)

            # Evaluate envelope (same logic as params.set_L_params)
            env_vals = []
            for params in envelope:
                A, B, X = params
                i = 0
                while X[i] < log_rT:
                    i += 1
                env_vals.append(log_rT * A[i] + B[i] + log_2p)
            log_wL_min_env, log_wL_max_env = env_vals

            # Compute actual per-line widths
            gamma_min = np.min(gamma_arr, axis=0)
            gamma_max = np.max(gamma_arr, axis=0)
            actual_min = (na * log_rT + np.log(gamma_min) + log_2p).min()
            actual_max = (na * log_rT + np.log(gamma_max) + log_2p).max()

            assert log_wL_min_env <= actual_min + 1e-4, (
                f"Envelope min {log_wL_min_env} > actual min {actual_min} "
                f"at T={T}, p={p}"
            )
            assert log_wL_max_env >= actual_max - 1e-4, (
                f"Envelope max {log_wL_max_env} < actual max {actual_max} "
                f"at T={T}, p={p}"
            )


@pytest.mark.fast
def test_gpu_params_uses_shared_envelope():
    """radis.gpu.params delegates to radis.lbl.envelope."""
    import radis.gpu.params as params_mod
    from radis.gpu.params import init_G_params, init_L_params
    from radis.lbl.envelope import (
        compute_gaussian_envelope,
        compute_lorentzian_envelope,
    )

    np.random.seed(123)
    N = 50
    log_2vMm = np.random.uniform(-3, -1, size=N).astype(np.float32)
    na = np.random.uniform(0.4, 0.8, size=N).astype(np.float32)
    gamma_arr = np.random.uniform(0.02, 0.1, size=(2, N)).astype(np.float32)

    g_direct = compute_gaussian_envelope(log_2vMm)
    l_direct = compute_lorentzian_envelope(na, gamma_arr)

    init_G_params(log_2vMm)
    init_L_params(na, gamma_arr)

    assert np.isclose(g_direct[0], params_mod._G_param_data[0], atol=1e-6)
    assert np.isclose(g_direct[1], params_mod._G_param_data[1], atol=1e-6)

    for i in range(2):
        np.testing.assert_allclose(
            l_direct[i][0], params_mod._L_param_data[i][0], atol=1e-6
        )
        np.testing.assert_allclose(
            l_direct[i][1], params_mod._L_param_data[i][1], atol=1e-6
        )
        np.testing.assert_allclose(
            l_direct[i][2], params_mod._L_param_data[i][2], atol=1e-6
        )


@pytest.mark.fast
def test_envelope_to_metadata():
    """envelope_to_metadata produces correct HDF5-storable dict."""
    from radis.lbl.envelope import (
        compute_gaussian_envelope,
        compute_lorentzian_envelope,
        envelope_to_metadata,
    )

    np.random.seed(77)
    N = 80
    log_2vMm = np.random.uniform(-3, -1, size=N).astype(np.float32)
    na = np.random.uniform(0.3, 0.9, size=N).astype(np.float32)
    gamma_arr = np.random.uniform(0.02, 0.12, size=(2, N)).astype(np.float32)

    g_env = compute_gaussian_envelope(log_2vMm)
    l_env = compute_lorentzian_envelope(na, gamma_arr)

    meta = envelope_to_metadata(g_env, l_env)

    assert np.isclose(meta["envelope_g_min"], g_env[0])
    assert np.isclose(meta["envelope_g_max"], g_env[1])

    l_min_d = json.loads(meta["envelope_l_min"])
    np.testing.assert_allclose(l_min_d["A"], l_env[0][0].tolist(), atol=1e-6)
    np.testing.assert_allclose(l_min_d["B"], l_env[0][1].tolist(), atol=1e-6)

    l_max_d = json.loads(meta["envelope_l_max"])
    np.testing.assert_allclose(l_max_d["A"], l_env[1][0].tolist(), atol=1e-6)
    np.testing.assert_allclose(l_max_d["B"], l_env[1][1].tolist(), atol=1e-6)


if __name__ == "__main__":
    test_gaussian_envelope_basic()
    test_gaussian_envelope_single_line()
    test_lorentzian_envelope_single_na()
    test_lorentzian_envelope_bounds()
    test_gpu_params_uses_shared_envelope()
    test_envelope_to_metadata()
    print("All envelope tests passed!")
