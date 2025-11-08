# SIMD Parser for RADIS

High-performance SIMD-accelerated parser for HITRAN/HITEMP data files.

## Overview

This module provides a C++ SIMD-accelerated parser that significantly speeds up parsing of large HITRAN/HITEMP `.par` files. It uses AVX2/SSE4.2 instructions and OpenMP for parallel processing.

## Features

- **SIMD acceleration**: Uses AVX2/SSE4.2 instructions for vectorized parsing
- **Multi-threading**: OpenMP parallelization for multi-core CPUs
- **Automatic fallback**: Falls back to standard Python parser if SIMD is unavailable
- **Seamless integration**: Automatically used by `fetch_hitemp()` when available

## Compilation

The parser is automatically compiled when first used. To manually compile:

```bash
cd radis/simd_parser
python compile_simd_parser.py
```

### Requirements

- C++ compiler with C++17 support (g++, clang++, or MSVC)
- OpenMP support
- AVX2/SSE4.2 CPU instructions (most modern CPUs)

## Usage

The SIMD parser is automatically used by `fetch_hitemp()` for CO2 molecules:

```python
from radis import fetch_hitemp

# SIMD parser will be used automatically for CO2
df = fetch_hitemp("CO2", load_wavenum_min=2000, load_wavenum_max=2300)

# Disable SIMD parser if needed
df = fetch_hitemp("CO2", load_wavenum_min=2000, load_wavenum_max=2300, use_simd=False)
```

## Files

- `simd_parser.py`: Python interface to the SIMD parser
- `corrected_hitran_parser.cpp`: C++ SIMD parser implementation
- `compile_simd_parser.py`: Compilation script
- `hitran_simd_parser` or `hitran_simd_parser.exe`: Compiled executable (generated)

## Performance

The SIMD parser can be 5-10x faster than the standard Python parser for large files, depending on your CPU and data size.

## Troubleshooting

If compilation fails:
1. Ensure you have a C++ compiler installed (g++ on Linux/Mac, MSVC or MinGW on Windows)
2. Check that your CPU supports AVX2 instructions
3. The system will automatically fall back to the standard parser if SIMD is unavailable
