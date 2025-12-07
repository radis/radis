"""
SIMD Parser Module for RADIS

High-performance SIMD-accelerated parsing of HITRAN/HITEMP files.
"""

from .simd_parser import (
    compile_simd_parser_if_needed,
    is_simd_parser_available,
    parse_hitran_simd,
)

__all__ = [
    "compile_simd_parser_if_needed",
    "is_simd_parser_available",
    "parse_hitran_simd",
]
