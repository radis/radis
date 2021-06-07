# -*- coding: utf-8 -*-
"""
Created on Sat May  1 13:50:36 2021

@author: erwan
"""

import matplotlib.pyplot as plt
import numpy as np


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
