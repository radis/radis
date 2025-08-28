#!/usr/bin/env python3
"""
Compilation script for the SIMD HITRAN parser
"""

import os
import platform
import subprocess
import sys


def get_compiler_flags():
    """Get optimized compiler flags for the current platform"""
    system = platform.system().lower()

    if system == "windows":
        # Windows with MSVC or MinGW
        msvc_flags = [
            "/O2",  # Maximum optimization (MSVC)
            "/std:c++17",  # C++17 standard (MSVC)
            "/openmp",  # OpenMP support (MSVC)
            "/arch:AVX2",  # AVX2 support (MSVC)
        ]
        mingw_flags = [
            "-O3",  # Maximum optimization (MinGW)
            "-std=c++17",  # C++17 standard
            "-fopenmp",  # OpenMP support
            "-mavx2",  # AVX2 support
            "-msse4.2",  # SSE4.2 support
            "-mfma",  # FMA support
        ]
        return {"msvc": msvc_flags, "mingw": mingw_flags}
    else:
        # Linux/macOS with GCC/Clang
        flags = [
            "-O3",  # Maximum optimization
            "-std=c++17",  # C++17 standard
            "-fopenmp",  # OpenMP support
            "-march=native",  # Optimize for current CPU
            "-mavx2",  # AVX2 support
            "-msse4.2",  # SSE4.2 support
            "-mfma",  # FMA support
        ]
        return {"gcc": flags}


def compile_simd_parser(source_file="corrected_hitran_parser.cpp", output_file=None):
    """Compile the SIMD HITRAN parser"""

    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    source_path = os.path.join(script_dir, source_file)

    # Set platform-specific output file
    if output_file is None:
        if platform.system().lower() == "windows":
            output_file = "hitran_simd_parser.exe"
        else:
            output_file = "hitran_simd_parser"

    output_path = os.path.join(script_dir, output_file)

    if not os.path.exists(source_path):
        raise FileNotFoundError(f"Source file not found: {source_path}")

    # Detect compiler based on platform
    system = platform.system().lower()
    compilers = []

    if system == "windows":
        # Try MSVC first, then MinGW
        compilers = ["cl", "g++", "clang++"]
    else:
        # Linux/macOS
        compilers = ["g++", "clang++", "c++"]

    compiler = None
    is_msvc = False
    is_mingw = False

    for comp in compilers:
        try:
            subprocess.run(
                [comp, "--version"], capture_output=True, check=True, text=True
            )
            compiler = comp
            is_msvc = comp == "cl"
            is_mingw = comp == "g++" and system == "windows"
            break
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue

    if not compiler:
        raise RuntimeError(
            "No suitable C++ compiler found. Please install MSVC, MinGW, g++, or clang++"
        )

    # Get appropriate flags for the detected compiler
    flag_sets = get_compiler_flags()

    if is_msvc:
        flags = flag_sets["msvc"]
    elif is_mingw:
        flags = flag_sets["mingw"]
    else:
        flags = flag_sets["gcc"]

    # Build compilation command
    if is_msvc:
        # MSVC specific command structure
        cmd = [compiler] + flags + [source_path, f"/Fe:{output_path}"]
    else:
        # GCC/Clang structure
        cmd = [compiler] + flags + [source_path, "-o", output_path]

    print(f"Compiling SIMD HITRAN parser...")
    print(f"Platform: {platform.system()}")
    print(f"Compiler: {compiler}")
    print(f"Command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"✓ Successfully compiled: {output_path}")

        # Make executable on Unix-like systems
        if system != "windows":
            os.chmod(output_path, 0o755)

        return output_path

    except subprocess.CalledProcessError as e:
        print(f"✗ Compilation failed!")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise


def check_dependencies():
    """Check if required dependencies are available"""
    print("Checking dependencies...")
    system = platform.system().lower()

    # Check for OpenMP
    try:
        if system == "windows":
            # Windows OpenMP check
            print("✓ OpenMP support assumed (check compiler documentation)")
        else:
            result = subprocess.run(
                ["echo", "#include <omp.h>\nint main(){return 0;}"],
                shell=True,
                capture_output=True,
            )
            print("✓ OpenMP headers available")
    except:
        print("⚠ OpenMP may not be available")

    # Check CPU features
    try:
        if system == "linux":
            with open("/proc/cpuinfo", "r") as f:
                cpuinfo = f.read()
                if "avx2" in cpuinfo:
                    print("✓ AVX2 support detected")
                if "sse4_2" in cpuinfo:
                    print("✓ SSE4.2 support detected")
                if "fma" in cpuinfo:
                    print("✓ FMA support detected")
        elif system == "windows":
            # Windows CPU feature detection
            try:
                import cpuinfo

                info = cpuinfo.get_cpu_info()
                flags = info.get("flags", [])
                if "avx2" in flags:
                    print("✓ AVX2 support detected")
                if "sse4_2" in flags:
                    print("✓ SSE4.2 support detected")
                if "fma" in flags:
                    print("✓ FMA support detected")
            except ImportError:
                print(
                    "⚠ Could not detect CPU features (install py-cpuinfo for details)"
                )
            except Exception:
                print("⚠ Could not detect CPU features")
        elif system == "darwin":
            # macOS CPU feature detection
            try:
                result = subprocess.run(
                    ["sysctl", "-n", "machdep.cpu.features"],
                    capture_output=True,
                    text=True,
                )
                features = result.stdout.lower()
                if "avx2" in features:
                    print("✓ AVX2 support detected")
                if "sse4.2" in features:
                    print("✓ SSE4.2 support detected")
            except:
                print("⚠ Could not detect CPU features")
    except:
        print("⚠ Could not detect CPU features")


if __name__ == "__main__":
    try:
        check_dependencies()
        executable_path = compile_simd_parser()

        print(f"\n=== SIMD Parser Ready ===")
        print(f"Executable: {executable_path}")
        print(f"Usage: {executable_path} <input.par>")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
