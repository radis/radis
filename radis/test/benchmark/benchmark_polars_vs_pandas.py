"""
Benchmark: Polars vs Pandas for RADIS spectroscopic database loading.

Compares loading and filtering performance to evaluate Polars as a
replacement for Vaex in RADIS (see issue #658).

Run:
    python radis/test/benchmarks/benchmark_polars_vs_pandas.py

Requires:
    pip install polars pyarrow pandas
"""

import time
import os
import tempfile
import numpy as np
import pandas as pd

try:
    import polars as pl

    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False
    print("polars not installed, skipping Polars benchmarks")

try:
    import pyarrow.parquet as pq

    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False
    print("pyarrow not installed, skipping PyArrow benchmarks")


# RADIS column names from radis/api/columns_hitran.py
COLUMNS = ["wav", "int", "A", "airbrd", "selbrd", "El",
           "Tdpair", "Pshft", "gp", "gpp", "iso", "id"]


def generate_hitran_data(n_lines):
    """Generate synthetic spectroscopic data matching HITRAN column format."""
    rng = np.random.default_rng(42)
    data = {
        "wav": np.sort(rng.uniform(1, 25000, n_lines)),
        "int": rng.lognormal(-20, 5, n_lines),
        "A": rng.lognormal(1, 2, n_lines),
        "airbrd": rng.uniform(0.01, 0.15, n_lines),
        "selbrd": rng.uniform(0.01, 0.15, n_lines),
        "El": rng.uniform(0, 5000, n_lines),
        "Tdpair": rng.uniform(0.4, 0.8, n_lines),
        "Pshft": rng.uniform(-0.01, 0.01, n_lines),
        "gp": rng.integers(1, 300, n_lines).astype(float),
        "gpp": rng.integers(1, 300, n_lines).astype(float),
        "iso": rng.integers(1, 4, n_lines).astype(float),
        "id": np.ones(n_lines),
    }
    return data


def write_parquet(data, path):
    df = pd.DataFrame(data)
    df = df.sort_values("wav").reset_index(drop=True)
    df.to_parquet(path, engine="pyarrow", row_group_size=50000)


def write_hdf5(data, path):
    df = pd.DataFrame(data)
    df.to_hdf(path, key="df", format="table", complevel=1)


# --- benchmark functions ---

def bench_pandas_full_load(hdf5_path):
    t0 = time.perf_counter()
    df = pd.read_hdf(hdf5_path)
    t = time.perf_counter() - t0
    mem = df.memory_usage(deep=True).sum() / 1e6
    return t, mem


def bench_pandas_filter(hdf5_path, wmin, wmax):
    t0 = time.perf_counter()
    df = pd.read_hdf(hdf5_path)
    filtered = df[(df["wav"] >= wmin) & (df["wav"] <= wmax)]
    t = time.perf_counter() - t0
    return t, len(filtered)


def bench_polars_full_load(parquet_path):
    t0 = time.perf_counter()
    df = pl.scan_parquet(parquet_path).collect()
    t = time.perf_counter() - t0
    mem = df.estimated_size("mb")
    return t, mem


def bench_polars_filter(parquet_path, wmin, wmax):
    t0 = time.perf_counter()
    df = pl.scan_parquet(parquet_path).filter(
        pl.col("wav").is_between(wmin, wmax)
    ).collect()
    t = time.perf_counter() - t0
    return t, len(df)


def bench_polars_filter_select(parquet_path, wmin, wmax, columns):
    t0 = time.perf_counter()
    df = (
        pl.scan_parquet(parquet_path)
        .select(columns)
        .filter(pl.col("wav").is_between(wmin, wmax))
        .collect()
    )
    t = time.perf_counter() - t0
    return t, len(df)


def bench_pyarrow_filter(parquet_path, wmin, wmax):
    t0 = time.perf_counter()
    table = pq.read_table(
        parquet_path,
        filters=[("wav", ">=", wmin), ("wav", "<=", wmax)],
    )
    t = time.perf_counter() - t0
    return t, len(table)


# --- main ---

def run_benchmarks():
    sizes = {
        "HITRAN CO (160K)": 160_000,
        "HITEMP CO (1.1M)": 1_100_000,
        "HITEMP H2O (5M)": 5_000_000,
    }

    wmin, wmax = 1900, 2300
    calc_columns = ["wav", "int", "A", "airbrd", "selbrd", "El"]

    results = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        for label, n_lines in sizes.items():
            print(f"\n{'='*60}")
            print(f"Database: {label} ({n_lines:,} lines)")
            print(f"{'='*60}")

            data = generate_hitran_data(n_lines)
            parquet_path = os.path.join(tmpdir, f"test_{n_lines}.parquet")
            hdf5_path = os.path.join(tmpdir, f"test_{n_lines}.h5")

            print("Writing test files...")
            write_parquet(data, parquet_path)
            write_hdf5(data, hdf5_path)

            pq_size = os.path.getsize(parquet_path) / 1e6
            h5_size = os.path.getsize(hdf5_path) / 1e6
            print(f"  Parquet: {pq_size:.1f} MB | HDF5: {h5_size:.1f} MB")

            result = {"n_lines": n_lines}

            # pandas + hdf5
            print("\n--- Pandas + HDF5 (current fallback) ---")
            t, mem = bench_pandas_full_load(hdf5_path)
            print(f"  Full load:       {t:.3f}s  ({mem:.1f} MB)")
            result["pandas_full_load"] = t

            t, n = bench_pandas_filter(hdf5_path, wmin, wmax)
            print(f"  Load+filter:     {t:.3f}s  ({n:,} lines matched)")
            result["pandas_filter"] = t

            # polars + parquet
            if HAS_POLARS:
                print("\n--- Polars + Parquet (proposed) ---")
                t, mem = bench_polars_full_load(parquet_path)
                print(f"  Full load:       {t:.3f}s  ({mem:.1f} MB)")
                result["polars_full_load"] = t

                t, n = bench_polars_filter(parquet_path, wmin, wmax)
                print(f"  Pushdown filter: {t:.3f}s  ({n:,} lines)")
                result["polars_filter"] = t

                t, n = bench_polars_filter_select(
                    parquet_path, wmin, wmax, calc_columns
                )
                print(f"  Filter+select:   {t:.3f}s  ({n:,} lines, 6 cols)")
                result["polars_filter_select"] = t

            # pyarrow
            if HAS_PYARROW:
                print("\n--- PyArrow + Parquet (alternative) ---")
                t, n = bench_pyarrow_filter(parquet_path, wmin, wmax)
                print(f"  Filter:          {t:.3f}s  ({n:,} lines)")
                result["pyarrow_filter"] = t

            results[label] = result

    # summary
    print(f"\n\n{'='*60}")
    print("SUMMARY: Wavenumber filter (1900-2300 cm-1)")
    print(f"{'='*60}")
    header = f"{'Database':<25} {'Pandas HDF5':>12} {'Polars Pqt':>12} {'Speedup':>10}"
    print(header)
    print("-" * len(header))
    for label, r in results.items():
        pt = r.get("pandas_filter", 0)
        lt = r.get("polars_filter", 0)
        sp = pt / lt if lt > 0 else 0
        print(f"{label:<25} {pt:>11.3f}s {lt:>11.3f}s {sp:>9.1f}x")

    return results


if __name__ == "__main__":
    run_benchmarks()