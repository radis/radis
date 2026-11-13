"""
Validation case: H2/He-broadened NH3 cross sections from four engine / line-list pairs
=====================================================================================

Computes the cross section of 14NH3 broadened by H2 (0.85) and He (0.15) with

    HITRAN  + HAPI         (hitran-api, absorptionCoefficient_Voigt with Diluent)
    HITRAN  + RADIS        (fetch_databank('hitran', extra_params='all'))
    ExoMol  + PyExoCross   (CoYuTe, .broad files for H2 and He)
    ExoMol  + RADIS        (CoYuTe, fetch_databank('exomol'), same files as PyExoCross)

at 270 K and 500 K, 1 bar and 0.01 bar, and plots them against each other. The two
engines on the same line list should agree closely; HITRAN and ExoMol are different
line lists with different broadening parameters, so they should agree only roughly.

Settings shared by all engines: Voigt profile, line-wing cut-off 25 cm-1 (Gharib-Nezhad
et al. 2024, RASTI 3, 44, recommended value for P <= 200 bar), no pressure shift (ExoMol
has none; HAPI is run with LineShift=False and the RADIS HITRAN shifts are set to 0),
same wavenumber grid. RADIS uses its default LDM optimization ('simple'; pass
``--optimization none`` for line by line, which took > 1 h and ~25 GB of memory per condition for CoYuTe on a laptop).
Cross sections are per 14NH3 molecule (the HITRAN abundance is
divided out of the HITRAN results and RADIS's ExoMol abundance scaling is removed).

Known differences, not errors:
  * RADIS renormalizes each truncated line profile to unit area; HAPI and PyExoCross do
    not (~0.2% at 1 bar with a 25 cm-1 cut-off).
  * The NH3 .broad files mix the a0 (J''), m0 (|m|) and m1 (|m|, K'') recipes. RADIS
    uses m0 over a0 and skips m1 (it does not load K''). PyExoCross would apply m1 but
    reads its first column as a signed m and ignores K'' (3-13% median difference), so
    its copy of the .broad files has the m1 rows removed: both engines use a0 + m0.

Needs network access the first time (HITRAN API; ~300 MB of CoYuTe .trans.bz2 and the
93 MB states file, downloaded once for RADIS and once for PyExoCross, because RADIS may
delete the .bz2 files after caching them); later runs reuse the files in ``--data-dir``::

    python validation_NH3_H2He_hitran_exomol.py --data-dir ~/nh3_validation
"""

import argparse
import os
import time
import shutil
import urllib.request

import numpy as np

MOLECULE = "NH3"
ISOTOPOLOGUE = "14N-1H3"
EXOMOL_DATASET = "CoYuTe"
HITRAN_ID, HITRAN_ISO = 11, 1
DILUENT = {"H2": 0.85, "He": 0.15}
CONDITIONS = [(270.0, 1.0), (270.0, 0.01), (500.0, 1.0), (500.0, 0.01)]  # (K, bar)
TRUNCATION = 25.0  # cm-1
TRANS_SEGMENT = 100  # cm-1, width of the CoYuTe .trans files
BAR_PER_ATM = 1.01325
EXOMOL_DB_URL = "https://www.exomol.com/db"

ENGINES = {  # label: (color, linestyle); fixed categorical order, dashed = RADIS
    "HITRAN, HAPI": ("#2a78d6", "-"),
    "HITRAN, RADIS": ("#eb6834", "--"),
    "ExoMol, PyExoCross": ("#1baf7a", "-"),
    "ExoMol, RADIS": ("#eda100", "--"),
}


def wavenumber_grid(wmin, wmax, wstep):
    return wmin + wstep * np.arange(int(round((wmax - wmin) / wstep)) + 1)


def hitran_abundance():
    from radis.db.molparam import MolParams

    return MolParams().get(MOLECULE, HITRAN_ISO, "abundance")


def exomol_abundance():
    from radis.db.molparam import MOLPARAMS_EXTRA_PATH, MolParams

    return MolParams(extra_file_json=MOLPARAMS_EXTRA_PATH).get(
        MOLECULE, HITRAN_ISO, "abundance"
    )


# %% HAPI
def hapi_xsec(nu, T, P, data_dir):
    import hapi

    hapi.db_begin(os.path.join(data_dir, "hapi"))
    if MOLECULE not in hapi.tableList():
        hapi.fetch(
            MOLECULE,
            HITRAN_ID,
            HITRAN_ISO,
            nu[0] - TRUNCATION,
            nu[-1] + TRUNCATION,
            ParameterGroups=["160-char", "voigt_h2", "voigt_he"],
        )
    _, xs = hapi.absorptionCoefficient_Voigt(
        SourceTables=MOLECULE,
        Components=[(HITRAN_ID, HITRAN_ISO)],
        Environment={"T": T, "p": P / BAR_PER_ATM},
        Diluent=DILUENT,
        WavenumberGrid=nu,
        WavenumberWing=TRUNCATION,
        LineShift=False,
        HITRAN_units=True,
    )
    return xs / hitran_abundance()


# %% RADIS
def radis_xsec(nu, T, P, databank, data_dir, cutoff, optimization="simple"):
    import radis
    from radis import SpectrumFactory

    radis.config["DEFAULT_DOWNLOAD_PATH"] = os.path.join(data_dir, "radisdb")
    x = 1e-6  # trace NH3: self-broadening negligible, as in the other engines
    sf = SpectrumFactory(
        wavenum_min=nu[0],
        wavenum_max=nu[-1],
        wstep=nu[1] - nu[0],
        molecule=MOLECULE,
        isotope="1",
        pressure=P,
        mole_fraction=x,
        diluent={k: v * (1 - x) for k, v in DILUENT.items()},
        truncation=TRUNCATION,
        neighbour_lines=TRUNCATION,
        cutoff=cutoff if databank == "exomol" else 0,
        optimization=optimization,
        broadening_method="convolve",
        verbose=0,
        warnings={"AccuracyError": "ignore", "default": "ignore"},
    )
    if databank == "hitran":
        # "diluent" keeps the gamma_h2 / n_h2 / gamma_he / n_he columns
        sf.fetch_databank(
            "hitran",
            extra_params="all",
            load_columns=["equilibrium", "diluent"],
            db_use_cached=True,
        )
        for col in [c for c in sf.df0.columns if c == "Pshft" or c.startswith("delta_")]:
            sf.df0[col] = 0.0
        abundance = hitran_abundance()
    else:
        sf.fetch_databank("exomol", database=EXOMOL_DATASET)
        abundance = exomol_abundance()
    s = sf.eq_spectrum(Tgas=T)
    w, xs = s.get("xsection", wunit="cm-1", Iunit="cm2")
    return np.interp(nu, w, xs) / abundance


# %% PyExoCross
def exomol_dir(data_dir):
    return os.path.join(data_dir, "pyexocross")


def download(url, path):
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as r, open(path + ".part", "wb") as f:
        shutil.copyfileobj(r, f, length=1 << 20)
    os.replace(path + ".part", path)


def ensure_exomol_files(nu, data_dir):
    """PyExoCross's own copy of the files it reads: .def.json, .pf, states, the .trans
    segments that reach the grid within the cut-off, and the H2 / He .broad files"""
    base = f"{EXOMOL_DB_URL}/{MOLECULE}/{ISOTOPOLOGUE}"
    folder = os.path.join(exomol_dir(data_dir), MOLECULE, ISOTOPOLOGUE, EXOMOL_DATASET)
    prefix = f"{ISOTOPOLOGUE}__{EXOMOL_DATASET}"
    names = [f"{prefix}.def.json", f"{prefix}.pf", f"{prefix}.states.bz2"]
    lo = int((nu[0] - TRUNCATION) // TRANS_SEGMENT) * TRANS_SEGMENT
    for w0 in range(lo, int(nu[-1] + TRUNCATION), TRANS_SEGMENT):
        names.append(f"{prefix}__{w0:05d}-{w0 + TRANS_SEGMENT:05d}.trans.bz2")
    for name in names:
        download(f"{base}/{EXOMOL_DATASET}/{name}", os.path.join(folder, name))
    for broadener in DILUENT:
        name = f"{ISOTOPOLOGUE}__{broadener}.broad"
        path = os.path.join(folder, name)
        download(f"{base}/{name}", path)
        # same recipes as RADIS: without the m1 rows PyExoCross uses m0 there
        with open(path) as f:
            lines = f.readlines()
        keep = [l for l in lines if l.split()[:1] != ["m1"]]
        if len(keep) < len(lines):
            with open(path, "w") as f:
                f.writelines(keep)


def pyexocross_xsec(nu, T, P, data_dir, cutoff, threads):
    import pyexocross as px

    ensure_exomol_files(nu, data_dir)
    result = px.cross_sections(
        database="ExoMol",
        molecule=MOLECULE,
        isotopologue=ISOTOPOLOGUE,
        dataset=EXOMOL_DATASET,
        read_path=exomol_dir(data_dir) + "/",
        temperatures=[T],
        pressures=[P],
        min_range=nu[0],
        max_range=nu[-1],
        n_point=len(nu),
        wn_wl="WN",
        wn_wl_unit="cm-1",
        profile="SciPyVoigt",
        cutoff=TRUNCATION,
        threshold=cutoff,
        broadeners=list(DILUENT),
        ratios=list(DILUENT.values()),
        ncputrans=threads,
        ncpufiles=1,
        output="memory",
        verbose=False,
        log="none",
    )
    w = np.asarray(result.coords["wavenumber"])
    return np.interp(nu, w, np.asarray(result.select(T=T, P=P)))


# %% comparison
def compare(x, ref):
    """median |x/ref - 1| where ref > 1e-3 of its maximum, and the integral ratio"""
    m = ref > ref.max() * 1e-3
    return float(np.median(np.abs(x[m] / ref[m] - 1))), float(x.sum() / ref.sum())


PAIRS = [  # (label, numerator, denominator)
    ("RADIS / HAPI (HITRAN)", "HITRAN, RADIS", "HITRAN, HAPI"),
    ("RADIS / PyExoCross (ExoMol)", "ExoMol, RADIS", "ExoMol, PyExoCross"),
    ("ExoMol / HITRAN (PyExoCross / HAPI)", "ExoMol, PyExoCross", "HITRAN, HAPI"),
]


def plot(
    nu,
    results,
    out,
    title="14NH3 broadened by H2 (0.85) + He (0.15): Voigt, 25 cm$^{-1}$ wing cut-off, no shift",
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        2,
        len(CONDITIONS),
        figsize=(4.2 * len(CONDITIONS), 6.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1.4]},
        constrained_layout=True,
    )
    ratio_style = {PAIRS[0][0]: ENGINES["HITRAN, RADIS"],
                   PAIRS[1][0]: ENGINES["ExoMol, RADIS"],
                   PAIRS[2][0]: ("#6b6a63", "-")}
    for j, (T, P) in enumerate(CONDITIONS):
        res = results[(T, P)]
        ax, axr = axes[0, j], axes[1, j]
        for label, (color, ls) in ENGINES.items():
            ax.plot(nu, res[label], color=color, ls=ls, lw=1.2 if ls == "-" else 1.0,
                    label=label)
        ax.set_yscale("log")
        top = max(r.max() for r in res.values())
        ax.set_ylim(top * 1e-5, top * 3)
        ax.set_title(f"{T:g} K, {P:g} bar", fontsize=11)
        for label, num, den in PAIRS:
            ref = res[den]
            ratio = np.where(ref > ref.max() * 1e-3, res[num] / np.where(ref > 0, ref, 1), np.nan)
            color, ls = ratio_style[label]
            axr.plot(nu, ratio, color=color, ls=ls, lw=0.9, label=label)
        axr.axhline(1, color="#9a998f", lw=0.6)
        axr.set_ylim(0.5, 1.5)
        axr.set_xlabel("Wavenumber (cm$^{-1}$)")
        for a in (ax, axr):
            a.grid(color="#e6e5df", lw=0.5)
            a.tick_params(labelsize=9)
    axes[0, 0].set_ylabel("Cross section (cm$^2$/molecule)")
    axes[1, 0].set_ylabel("Ratio")
    axes[0, 0].legend(fontsize=8, frameon=False, loc="lower left")
    axes[1, 0].legend(fontsize=7, frameon=False, loc="lower left")
    fig.suptitle(title, fontsize=12)
    fig.savefig(out, dpi=150)
    print(f"figure: {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default=os.path.abspath("NH3_validation_data"),
                   help="where HITRAN (HAPI), RADIS and ExoMol files are stored")
    p.add_argument("--wmin", type=float, default=955.0, help="cm-1 (NH3 nu2 band)")
    p.add_argument("--wmax", type=float, default=975.0, help="cm-1")
    p.add_argument("--wstep", type=float, default=2e-4,
                   help="cm-1; resolves the Doppler width (~1.4e-3 cm-1 HWHM at 270 K)")
    p.add_argument("--cutoff", type=float, default=1e-27,
                   help="ExoMol line intensity threshold at T (cm/molecule), RADIS and PyExoCross")
    p.add_argument("--threads", type=int, default=os.cpu_count())
    p.add_argument("--optimization", default="simple",
                   help="RADIS optimization: 'simple' (LDM, default), 'min-RMS' or 'none' (line by line)")
    p.add_argument("--out", default="validation_NH3_H2He_hitran_exomol.png")
    args = p.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    optimization = None if args.optimization.lower() == "none" else args.optimization

    nu = wavenumber_grid(args.wmin, args.wmax, args.wstep)
    results, rows = {}, []
    for T, P in CONDITIONS:
        res, seconds = {}, {}
        for label, fn in [
            ("HITRAN, HAPI", lambda: hapi_xsec(nu, T, P, args.data_dir)),
            ("HITRAN, RADIS", lambda: radis_xsec(nu, T, P, "hitran", args.data_dir, args.cutoff, optimization)),
            ("ExoMol, RADIS", lambda: radis_xsec(nu, T, P, "exomol", args.data_dir, args.cutoff, optimization)),
            ("ExoMol, PyExoCross", lambda: pyexocross_xsec(nu, T, P, args.data_dir, args.cutoff, args.threads)),
        ]:
            t0 = time.time()
            res[label] = np.asarray(fn(), dtype=float)
            seconds[label] = time.time() - t0
        results[(T, P)] = res
        print(f"\n{T:g} K, {P:g} bar  (" + ", ".join(f"{k} {v:.0f} s" for k, v in seconds.items()) + ")")
        for label, num, den in PAIRS:
            med, integral = compare(res[num], res[den])
            rows.append((T, P, label, med, integral))
            print(f"  {label:38s} median |ratio - 1| {med:8.2e}   integral ratio {integral:.4f}")

    out_csv = os.path.splitext(args.out)[0] + ".csv"
    with open(out_csv, "w") as f:
        f.write("temperature_K,pressure_bar,pair,median_rel_diff,integral_ratio\n")
        for r in rows:
            f.write(f"{r[0]:g},{r[1]:g},{r[2]},{r[3]:.4e},{r[4]:.5f}\n")
    print(f"table: {out_csv}")
    plot(nu, results, args.out)


if __name__ == "__main__":
    main()
