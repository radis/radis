"""
Demo: Polars/Parquet lazy-loading backend for RADIS (ref issue #978, #658)

This script demonstrates the DataFrameAdapter pattern working end-to-end
with synthetic HITRAN-format spectroscopic data.

Usage:
    python3.11 radis/io/demo_polars_backend.py

What this demo shows:
    1. Creates a synthetic HITEMP-like database (1M lines, 12 columns)
    2. Saves it as both HDF5 (current format) and Parquet (proposed format)
    3. Runs a typical calc_spectrum loading workflow through both adapters
    4. Compares load time, filter time, and verifies data matches
"""

import time
import os
import tempfile
import numpy as np
import pandas as pd


def generate_hitemp_co_data(n_lines=1_000_000):
    """Generate synthetic HITEMP CO database matching real column format."""
    print(f"Generating synthetic HITEMP CO database ({n_lines:,} lines)...")
    rng = np.random.default_rng(42)
    data = {
        "wav": np.sort(rng.uniform(1, 9000, n_lines)),
        "int": rng.lognormal(-20, 5, n_lines),
        "A": rng.lognormal(1, 2, n_lines),
        "airbrd": rng.uniform(0.01, 0.15, n_lines),
        "selbrd": rng.uniform(0.01, 0.15, n_lines),
        "El": rng.uniform(0, 5000, n_lines),
        "Tdpair": rng.uniform(0.4, 0.8, n_lines),
        "Pshft": rng.uniform(-0.01, 0.01, n_lines),
        "gp": rng.integers(1, 300, n_lines).astype(float),
        "gpp": rng.integers(1, 300, n_lines).astype(float),
        "iso": rng.choice([1.0, 2.0, 3.0], n_lines, p=[0.7, 0.2, 0.1]),
        "id": np.ones(n_lines),
    }
    return pd.DataFrame(data)


def demo_current_workflow(hdf5_path, wmin, wmax):
    """Current RADIS: load everything into RAM, then filter in Python."""
    print("\n--- Current RADIS workflow (Pandas + HDF5) ---")
    t0 = time.perf_counter()
    df = pd.read_hdf(hdf5_path)
    t_load = time.perf_counter() - t0
    mem = df.memory_usage(deep=True).sum() / 1e6

    t1 = time.perf_counter()
    filtered = df[(df["wav"] >= wmin) & (df["wav"] <= wmax)]
    filtered = filtered[filtered["iso"] == 1.0]
    result = filtered[["wav", "int", "A", "airbrd", "selbrd", "El"]]
    t_filter = time.perf_counter() - t1
    t_total = time.perf_counter() - t0

    print(f"  Load entire HDF5:   {t_load:.3f}s ({mem:.1f} MB in RAM)")
    print(f"  Filter + select:    {t_filter:.3f}s")
    print(f"  Total:              {t_total:.3f}s")
    print(f"  Lines loaded:       {len(df):,} (full database)")
    print(f"  Lines after filter: {len(result):,}")
    return result, t_total


def demo_polars_workflow(parquet_path, wmin, wmax):
    """Proposed: lazy load with predicate pushdown via PolarsAdapter."""
    from radis.io.polars_adapter import PolarsAdapter

    print("\n--- Proposed workflow (PolarsAdapter + Parquet) ---")
    t0 = time.perf_counter()
    adapter = PolarsAdapter()
    adapter.load(parquet_path)
    adapter.filter_range("wav", wmin, wmax)
    adapter.filter_equals("iso", 1.0)
    adapter.select_columns(["wav", "int", "A", "airbrd", "selbrd", "El"])
    result = adapter.to_pandas()
    t_total = time.perf_counter() - t0

    print(f"  Total (lazy+exec):  {t_total:.3f}s")
    print(f"  Lines returned:     {len(result):,}")
    print(f"  Memory: only filtered data loaded (predicate pushdown)")
    return result, t_total


def demo_migration(hdf5_path):
    """Show automatic HDF5 to Parquet migration."""
    from radis.io.hdf5_to_parquet import convert_hdf5_to_parquet

    print("\n--- HDF5 to Parquet migration ---")
    parquet_path = hdf5_path.replace(".h5", ".parquet")
    t0 = time.perf_counter()
    convert_hdf5_to_parquet(hdf5_path, parquet_path, sort_by="wav",
                            compression="zstd", validate=True, verbose=False)
    t = time.perf_counter() - t0
    h5 = os.path.getsize(hdf5_path) / 1e6
    pq = os.path.getsize(parquet_path) / 1e6
    print(f"  Time:     {t:.3f}s")
    print(f"  HDF5:     {h5:.1f} MB -> Parquet: {pq:.1f} MB ({(1-pq/h5)*100:.0f}% smaller)")
    print(f"  Sorted by 'wav' for predicate pushdown")
    print(f"  Validated: rtol=1e-14")
    return parquet_path


def demo_adapter_factory():
    """Show engine selection."""
    from radis.io.polars_adapter import get_adapter
    print("\n--- Engine selection ---")
    a1 = get_adapter("polars")
    a2 = get_adapter("pytables")
    print(f"  get_adapter('polars')   -> {type(a1).__name__}")
    print(f"  get_adapter('pytables') -> {type(a2).__name__}")
    print(f"  Production: reads config['DATAFRAME_ENGINE'] from radis.json")


def main():
    print("=" * 65)
    print("DEMO: Polars/Parquet Backend for RADIS")
    print("Issue #978, ref #658")
    print("=" * 65)

    wmin, wmax = 1900, 2300  # CO fundamental band

    with tempfile.TemporaryDirectory() as tmpdir:
        df = generate_hitemp_co_data(1_000_000)
        hdf5_path = os.path.join(tmpdir, "CO-HITEMP.h5")
        df.to_hdf(hdf5_path, key="df", format="table", complevel=1)

        result_pd, t_pd = demo_current_workflow(hdf5_path, wmin, wmax)
        parquet_path = demo_migration(hdf5_path)
        result_pl, t_pl = demo_polars_workflow(parquet_path, wmin, wmax)

        print("\n--- Verification ---")
        wav_pd = np.sort(result_pd["wav"].values)
        wav_pl = np.sort(result_pl["wav"].values)
        match = len(wav_pd) == len(wav_pl) and np.allclose(wav_pd, wav_pl, rtol=1e-14)
        print(f"  Results match: {'YES' if match else 'NO'} ({len(wav_pd):,} lines)")

        demo_adapter_factory()

        speedup = t_pd / t_pl if t_pl > 0 else 0
        print(f"\n{'='*65}")
        print(f"RESULT: Polars is {speedup:.1f}x faster than Pandas for this workflow")
        print(f"{'='*65}")


if __name__ == "__main__":
    main()