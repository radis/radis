#!/usr/bin/env python3
"""
SIMD HITRAN Parser Integration Module for RADIS

Provides high-performance SIMD-accelerated parsing of HITRAN/HITEMP files
as an optional alternative to the standard Python parser.
"""

import os
import platform
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd


def get_simd_parser_path():
    """Get the path to the SIMD parser executable"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if platform.system().lower() == "windows":
        return os.path.join(current_dir, "hitran_simd_parser.exe")
    else:
        return os.path.join(current_dir, "hitran_simd_parser")


def is_simd_parser_available():
    """Check if the SIMD parser executable is available"""
    parser_path = get_simd_parser_path()
    return os.path.exists(parser_path) and os.access(parser_path, os.X_OK)


def compile_simd_parser_if_needed(force_recompile=False):
    """Compile the SIMD parser if it doesn't exist or if forced"""
    if force_recompile or not is_simd_parser_available():
        compile_script = os.path.join(
            os.path.dirname(get_simd_parser_path()), "compile_simd_parser.py"
        )

        if not os.path.exists(compile_script):
            warnings.warn(f"Compilation script not found: {compile_script}")
            return False

        try:
            subprocess.run(
                [sys.executable, compile_script],
                capture_output=True,
                text=True,
                check=True,
                cwd=os.path.dirname(compile_script),
            )
            return True
        except subprocess.CalledProcessError as e:
            warnings.warn(f"Failed to compile SIMD parser: {e.stderr}")
            return False
    return True


def read_simd_binary_output(output_dir, verbose=True):
    """Read binary data files created by the SIMD C++ HITRAN parser"""
    if verbose:
        print(f"Reading SIMD binary data from {output_dir}/")

    if not os.path.exists(output_dir):
        raise FileNotFoundError(f"SIMD output directory not found: {output_dir}")

    def read_numeric_file(filename, dtype):
        filepath = os.path.join(output_dir, filename)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"SIMD output file not found: {filepath}")
        with open(filepath, "rb") as f:
            return np.frombuffer(f.read(), dtype=dtype)

    def read_string_file(filename):
        filepath = os.path.join(output_dir, filename)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"SIMD output file not found: {filepath}")

        strings = []
        with open(filepath, "rb") as f:
            while True:
                len_bytes = f.read(8)
                if len(len_bytes) < 8:
                    break

                str_len = int.from_bytes(len_bytes, byteorder="little")
                if str_len == 0:
                    strings.append("")
                    continue

                str_bytes = f.read(str_len)
                if len(str_bytes) < str_len:
                    break

                strings.append(str_bytes.decode("utf-8"))

        return strings

    # Load all numerical columns
    id_data = read_numeric_file("id.dat", np.int32)
    iso_data = read_numeric_file("iso.dat", np.int32)
    wav_data = read_numeric_file("wav.dat", np.float32)
    int_data = read_numeric_file("int.dat", np.float32)
    A_data = read_numeric_file("A.dat", np.float32)
    airbrd_data = read_numeric_file("airbrd.dat", np.float32)
    selbrd_data = read_numeric_file("selbrd.dat", np.float32)
    El_data = read_numeric_file("El.dat", np.float32)
    Tdpair_data = read_numeric_file("Tdpair.dat", np.float32)
    Pshft_data = read_numeric_file("Pshft.dat", np.float32)
    gp_data = read_numeric_file("gp.dat", np.float32)
    gpp_data = read_numeric_file("gpp.dat", np.float32)

    # Load all string columns
    globu_data = read_string_file("globu.dat")
    globl_data = read_string_file("globl.dat")
    locu_data = read_string_file("locu.dat")
    locl_data = read_string_file("locl.dat")
    ierr_data = read_string_file("ierr.dat")
    iref_data = read_string_file("iref.dat")
    lmix_data = read_string_file("lmix.dat")

    # Verify all arrays have the same length
    record_count = len(id_data)
    arrays_to_check = [
        ("iso", iso_data),
        ("wav", wav_data),
        ("int", int_data),
        ("A", A_data),
        ("airbrd", airbrd_data),
        ("selbrd", selbrd_data),
        ("El", El_data),
        ("Tdpair", Tdpair_data),
        ("Pshft", Pshft_data),
        ("gp", gp_data),
        ("gpp", gpp_data),
        ("globu", globu_data),
        ("globl", globl_data),
        ("locu", locu_data),
        ("locl", locl_data),
        ("ierr", ierr_data),
        ("iref", iref_data),
        ("lmix", lmix_data),
    ]

    for name, data in arrays_to_check:
        if len(data) != record_count:
            raise ValueError(
                f"SIMD parser array length mismatch: {name} has {len(data)} records, expected {record_count}"
            )

    # Create DataFrame with exact RADIS column names and order
    df = pd.DataFrame(
        {
            "id": id_data,
            "iso": iso_data,
            "wav": wav_data.astype(
                np.float64
            ),  # Convert to double precision for RADIS compatibility
            "int": int_data.astype(np.float64),
            "A": A_data.astype(np.float64),
            "airbrd": airbrd_data.astype(np.float64),
            "selbrd": selbrd_data.astype(np.float64),
            "El": El_data.astype(np.float64),
            "Tdpair": Tdpair_data.astype(np.float64),
            "Pshft": Pshft_data.astype(np.float64),
            "globu": globu_data,
            "globl": globl_data,
            "locu": locu_data,
            "locl": locl_data,
            "ierr": ierr_data,
            "iref": iref_data,
            "lmix": lmix_data,
            "gp": gp_data.astype(np.float64),
            "gpp": gpp_data.astype(np.float64),
        }
    )

    if verbose:
        print(f"Successfully loaded {len(df)} HITRAN records using SIMD parser")

    return df


def run_simd_parser(input_file, output_dir=None, verbose=True):
    """Run the SIMD parser on a HITRAN/HITEMP .par file"""
    if not is_simd_parser_available():
        if not compile_simd_parser_if_needed():
            raise RuntimeError("SIMD parser is not available and compilation failed")

    parser_path = get_simd_parser_path()

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Create temporary output directory if not specified
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix="simd_hitran_")
    else:
        os.makedirs(output_dir, exist_ok=True)

    # The C++ parser writes to "hitemp_output" in the current directory
    original_cwd = os.getcwd()

    try:
        os.chdir(output_dir)

        if verbose:
            print(f"Running SIMD parser on {input_file}")
            print(f"Output directory: {output_dir}")

        result = subprocess.run(
            [parser_path, input_file],
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
        )

        if result.returncode != 0:
            raise RuntimeError(f"SIMD parser failed: {result.stderr}")

        if verbose and result.stdout:
            print("SIMD parser output:")
            print(result.stdout)

        # Move files from hitemp_output to main directory
        hitemp_output_dir = os.path.join(output_dir, "hitemp_output")
        if os.path.exists(hitemp_output_dir):
            for filename in os.listdir(hitemp_output_dir):
                src = os.path.join(hitemp_output_dir, filename)
                dst = os.path.join(output_dir, filename)
                shutil.move(src, dst)
            os.rmdir(hitemp_output_dir)

        # Also clean up any hitemp_output in the current working directory
        cwd_hitemp_output = os.path.join(os.getcwd(), "hitemp_output")
        if os.path.exists(cwd_hitemp_output):
            shutil.rmtree(cwd_hitemp_output)
            if verbose:
                print("Cleaned up hitemp_output directory")

        # Verify output files exist
        expected_files = [
            "id.dat",
            "iso.dat",
            "wav.dat",
            "int.dat",
            "A.dat",
            "airbrd.dat",
            "selbrd.dat",
            "El.dat",
            "Tdpair.dat",
            "Pshft.dat",
            "globu.dat",
            "globl.dat",
            "locu.dat",
            "locl.dat",
            "ierr.dat",
            "iref.dat",
            "lmix.dat",
            "gp.dat",
            "gpp.dat",
        ]

        missing_files = [
            f for f in expected_files if not os.path.exists(os.path.join(output_dir, f))
        ]
        if missing_files:
            raise RuntimeError(
                f"SIMD parser did not generate expected files: {missing_files}"
            )

        return output_dir

    finally:
        os.chdir(original_cwd)

        # Always clean up any hitemp_output directory in original working directory
        cwd_hitemp_output = os.path.join(original_cwd, "hitemp_output")
        if os.path.exists(cwd_hitemp_output):
            try:
                shutil.rmtree(cwd_hitemp_output)
            except:
                pass  # Ignore cleanup errors


def parse_hitran_simd(input_file, output_dir=None, verbose=True, cleanup=True):
    """Parse a HITRAN file using the SIMD parser and return a DataFrame"""
    temp_dir_created = output_dir is None

    try:
        # Run SIMD parser
        output_dir = run_simd_parser(input_file, output_dir, verbose)

        # Read the binary output
        df = read_simd_binary_output(output_dir, verbose)

        return df

    finally:
        # Clean up temporary directory if we created it
        if cleanup and temp_dir_created and output_dir and os.path.exists(output_dir):
            try:
                shutil.rmtree(output_dir)
                if verbose:
                    print(f"Cleaned up temporary directory: {output_dir}")
            except:
                warnings.warn(f"Could not clean up temporary directory: {output_dir}")
