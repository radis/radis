# -*- coding: utf-8 -*-
"""
Created on Sat May  1 13:50:36 2021

@author: erwan
"""

import numpy as np

import radis

KNOWN_CONTEXT = ["paper", "notebook", "talk", "poster"]


def _get_defaults(plotlib, context, style):

    expected_format = {
        "plot": {
            "plotlib": "ANY OF "
            + "/".join(['"publib"', '"seaborn"', '"matplotlib"', '"none"']),
            "context": "ANY OF " + "/".join(f'"{c}"' for c in KNOWN_CONTEXT),
            "style": "ANY OF the publib / seaborn / matplotlib themes",
        }
    }

    # Get defaults
    if plotlib == "config-default":
        try:
            plotlib = radis.config["plot"]["plotlib"]
        except KeyError:
            raise ValueError(
                "Missing key in your ~/radis.json user config file. Expected format :\n\n {0}".format(
                    expected_format
                )
            )

    if context == "config-default":
        try:
            context = radis.config["plot"]["context"]
        except KeyError:
            raise ValueError(
                "Missing key `context` in your ~/radis.json user config file. Expected format :\n\n {0}".format(
                    expected_format
                )
            )

    if style == "config-default":
        try:
            style = radis.config["plot"]["style"]
        except KeyError:
            raise ValueError(
                "Missing key `style` in your ~/radis.json user config file. Expected format :\n\n {0}".format(
                    expected_format
                )
            )

    if context and context not in KNOWN_CONTEXT:
        raise ValueError(f"context must be in {KNOWN_CONTEXT}. Got {context}")

    if plotlib == "publib":
        if not isinstance(style, list):
            style = [style]
        if context == "paper":
            style = style + ["article"]
        elif context in ["poster", "notebook", "talk"]:
            style = style + [context]

    return plotlib, context, style


def set_style(
    plotlib="config-default", context="config-default", style="config-default"
):
    """Set styles used by Radis plot functions (:py:func:`~radis.spectrum.spectrum.Spectrum.plot`,
    :py:func:`~radis.spectrum.compare.plot_diff`, :py:func:`~radis.tools.slit.plot_slit`, etc.)

    Parameters
    ----------
    plotlib: ``"publib"``, ``"seaborn"``, ``"matplotlib"``
        `publib <https://github.com/erwanp/publib>`__ , :py:mod:`seaborn`, :py:mod:`matplotlib`
    context: ``"poster"``, ``"talk"``, ``"paper"``, ``"notebook"``
        See :py:func:`seaborn.set_theme`
    style: a ``publib`` or ``matplotlib`` style, or a ``seaborn`` theme
        if ``"oring"``, use your current defaults.
        if ``"none"``, use your current defaults.

    Examples
    --------

    .. minigallery:: radis.misc.plot.set_style

    See also
    --------
    :py:func:`seaborn.set_theme`,
    :py:func:`publib.main.set_style`
    """

    plotlib, context, style = _get_defaults(plotlib, context, style)

    if plotlib != "none":
        import matplotlib

        matplotlib.rc_file_defaults()

    # Apply style
    if plotlib == "publib":
        import publib

        publib.set_style(style)
    elif plotlib == "seaborn":
        import seaborn as sns

        sns.set_theme(context=context, style=style)
    elif plotlib == "matplotlib":
        import matplotlib.pyplot as plt

        plt.style.use(style)
    elif plotlib == "none":
        pass
    else:
        raise ValueError(f"plot library: {plotlib}")


def fix_style(
    plotlib="config-default",
    context="config-default",
    style="config-default",
    *args,
    **kwargs,
):
    """A posteriori improvement of the last figure generated

    Effective when ``plotlib`` is `publib <https://github.com/erwanp/publib>`__

    Parameters
    ----------
    plotlib: ``"publib"``, ``"seaborn"``, ``"matplotlib"``
        `publib <https://github.com/erwanp/publib>`__ , :py:mod:`seaborn`, :py:mod:`matplotlib`
    context: ``"poster"``, ``"talk"``, ``"paper"``, ``"notebook"``
        See :py:func:`seaborn.set_theme`
    style: a ``publib`` or ``matplotlib`` style, or a ``seaborn`` theme
        if ``"oring"``, use your current defaults.
        if ``"none"``, use your current defaults.

    Examples
    --------

    .. minigallery:: radis.misc.plot.fix_style

    See also
    --------
    :py:func:`publib.main.fix_style`
    """

    plotlib, context, style = _get_defaults(plotlib, context, style)

    # Apply style
    if plotlib == "publib":
        import publib

        publib.fix_style(style, *args, **kwargs)
    elif plotlib == "seaborn":
        pass
    elif plotlib == "matplotlib":
        pass
    elif plotlib == "none":
        pass
    else:
        raise ValueError(f"plot library: {plotlib}")


def split_and_plot_by_parts(w, I, *args, **kwargs):
    """Plot two discontinued arrays (typically a spectrum) without showing
    junctions: first identify junctions then split and plot separately.

    Useful for plotting an experimental spectrum defined on different, non overlapping
    ranges without showing connecting lines between the ranges, or to plot an
    experimental spectrum defined on overlapping ranges, without showing connecting
    lines neither.

    Parameters
    ----------
    w, I: arrays
        typically output of :py:func:`~numpy.hstack`.

    Other Parameters
    ----------------
    split_threshold: int
        number of standard deviation for threshold. Default 10
    ax: matplotlib axe
        plot on a particular axe
    kwargs: dict
        forwarded to :func:`~matplotlib.pyplot.plot`
    cutwings: int
        discard elements on the side. Default 0

    Returns
    -------
    ax.plot
    """

    import matplotlib.pyplot as plt
    from publib.tools import keep_color

    # Get defaults
    ax = kwargs.pop("ax", None)
    if ax is None:
        ax = plt.gca()
    split_threshold = kwargs.pop("split_threshold", 10)  # type: int
    cutwings = kwargs.pop("cutwings", 0)  # type: int
    label = kwargs.pop("label", None)  # type: str

    # identify joints
    dw = np.diff(w)
    dwmean = dw.mean()
    joints = np.argwhere((abs(dw - dwmean) > split_threshold * dw.std())) + 1

    # Split
    if len(joints) > 0:
        ws = np.split(w, joints.flatten())
        Is = np.split(I, joints.flatten())

        # Plot separately
        out = []
        for i, (wi, Ii) in enumerate(zip(ws, Is)):
            if cutwings:
                wi = wi[cutwings:-cutwings]
                Ii = Ii[cutwings:-cutwings]
            if i == 0:  # label once only
                out.append(
                    ax.plot(
                        wi, Ii, *args, **dict(list(kwargs.items()) + [("label", label)])
                    )
                )
            else:
                keep_color()
                out.append(ax.plot(wi, Ii, *args, **kwargs))

        return list(zip(*out))

    else:
        if cutwings:
            w = w[cutwings:-cutwings]
            I = I[cutwings:-cutwings]
        return ax.plot(w, I, *args, **dict(list(kwargs.items()) + [("label", label)]))
