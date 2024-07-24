.. _label_radis_gpu:

==============================
RADIS-GPU Spectrum Calculation
==============================

RADIS provides GPU acceleration to massively speedup spectral computations.
Currently only Nvidia GPU's are supported, but this will likely change in the future
by using Vulkan as a backend.

Generally GPU computations are memory bandwidth limited, meaning the computation time of
a single spectrum is determined by the time it takes to move the database data from host
(=CPU) to device (=GPU) memory. Because of this, GPU computations take place in two steps:
An initialization step :py:func:`~radis.gpu.gpu.gpu_init` where, among other things, the database is
uploaded to the GPU, and an iteration step :py:func:`~radis.gpu.gpu.gpu_iterate`, where a new spectrum
with different parameters but the same database is computed. The latter could be repeated indefinitely
as long as the same database and spectral axis is used, resulting in extremely fast spectrum generation.

RADIS implements two functions that expose GPU functionality:

- :py:func:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu` computes a single spectrum and returns.

- :py:func:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu_interactive` computes a single
  spectrum and allows to interactively adjust parameters such as temperature, pressure, and
  composition, updating the spectrum in real time. Because the database only has to be transferred
  once, the updated spectra are calculated extremely fast.

By default both functions will be ran on a GPU if available. The CUDA code can also be compiled as pure
C++, which means it can be compiled for CPU in addition to GPU.
As a result, it ispossible to use the same GPU functions without an actual GPU by passing the
keyword ``backend='cpu-cuda'``, which forces use of the CPU targeted compiled code. This feature is
mostly for developers to check for errors in the CUDA code, but it can also be used for interactive
plotting on the CPU for small spectra.

GPU computation is currently only supported for equilibrium spectra. It is likely that
non-equilibrium spectra will be supported at some point in the future.


Single Spectrum
---------------

As mentioned above, the function :py:func:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu`
produces a single equilibrium spectrum using GPU acceleration. Below is a usage example::

    sf = SpectrumFactory(
        2100,
        2400,  # cm-1
        molecule="CO2",
        isotope="1,2,3",
        wstep=0.002,
    )

    sf.fetch_databank("hitemp")

    s = sf.eq_spectrum_gpu(
        Tgas=1100.0,  # K
        pressure=1,  # bar
        mole_fraction=0.8,
        path_length=0.2,  # cm
    )

    s.apply_slit(0.5)#cm-1
    s.plot("radiance", show=True)


.. minigallery:: radis.lbl.SpectrumFactory.eq_spectrum_gpu
    :add-heading:


Interactive Spectrum
--------------------

As mentioned before, computing the first GPU spectrum in a session takes a comparatively long time because the
entire database must be transferred to the GPU. The real power of GPU acceleration
becomes evident when computation times are not limited by data-transfer, i.e., when multiple
consecutive spectra are synthesized. One obvious use case would be the fitting of a spectrum.
Another one is interactive plotting, which can be done by calling
:py:func:`~radis.lbl.factory.eq_spectrum_gpu_interactive()`. A usage example is shown below::

    from radis.tools.plot_tools import ParamRange
    sf = SpectrumFactory(
        2100,
        2400,  # cm-1
        molecule="CO2",
        isotope="1,2,3",
        wstep=0.002,
    )

    sf.fetch_databank("hitemp")

    s = sf.eq_spectrum_gpu_interactive(
        var="radiance",
        Tgas=ParamRange(300.0, 2500.0, 1100.0),  # K
        pressure=ParamRange(0.1, 2, 1),  # bar
        mole_fraction=ParamRange(0, 1, 0.8),
        path_length=ParamRange(0, 1, 0.2),  # cm
        slit_FWHM=ParamRange(0, 1.5, 0.24),  # cm-1
        emulate=False,  # runs on GPU
        plotkwargs={"nfig": "same", "wunit": "nm"},
    )


.. minigallery:: radis.lbl.SpectrumFactory.eq_spectrum_gpu_interactive
    :add-heading:


Note that `eq_spectrum_gpu_interactive()` replaces all of `eq_spectrum_gpu()`,
`s.apply_slit()`, and `s.plot()` seen in the earlier example, and for this reason the
syntax is a little bit different. For example, we directly pass the `var` keyword to
`eq_spectrum_gpu_interactive()` to specify which spectrum should be plotted, and keyword arguments to `s.plot()`
are passed through `plotkwargs`.

Quantities that are to be varied must be initialized by a
:py:func:`~radis.tools.plot_tools.ParamRange` (valmin, valmax, valinit) object, which
takes the minimum value, maximum value, and init values of the scan range. Each `ParamRange()`
object will spawn a slider widget in the plot window with which the parameter can be
 interactively adjusted. The algorithm is extremely fast for a large number of lines (>100M)
 and will update with very low latency (<200ms typically). The code is not currently optimized
 for large wavenumber ranges (>500cm-1) however, which may take a bit longer (up to a couple seconds),
 provided the GPU didn't run out of memory.

At this moment the application of the instrumental function is done on the CPU to benefit from all features
already implemented in :py:func:`~radis.spectrum.Spectrum.apply_slit`. It is expected that these computations
will also move to the GPU at some point in the future.

Did you miss any feature implemented on GPU? or support for your particular system? The GPU code is heavily under development, so drop us a visit on [our Githup](https://github.com/radis/radis/issues/616) and let us know what you're looking for!



