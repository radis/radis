"""
Example: ExoJax-style usage of radis.api
==========================================
Shows how ExoJax can subclass radis.api classes directly.
Ref: https://github.com/radis/radis/issues/474
"""
from radis.api.exomolapi import MdbExomol, apply_jax_array_conversion


class ExoJaxMdbExomol(MdbExomol):
    """ExoJax-style subclass — uses radis.api directly"""

    def load_as_jax(self, wavenum_min, wavenum_max):
        """Load molecular data as JAX arrays (GPU-ready)"""
        return self.load(
            output="jax",
            load_wavenum_min=wavenum_min,
            load_wavenum_max=wavenum_max,
        )

# GPU arrays:  nu_lines, logsij0, elower  → jax arrays
# DRAM arrays: n_Texp, alpha_ref          → numpy arrays
