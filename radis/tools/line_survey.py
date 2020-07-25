# -*- coding: utf-8 -*-
"""
Functions to plot line surveys 


-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from radis.io.hitran import (
    get_molecule,
    HITRAN_CLASS1,
    HITRAN_CLASS2,
    HITRAN_CLASS3,
    HITRAN_CLASS4,
    HITRAN_CLASS5,
    HITRAN_CLASS6,
    HITRAN_CLASS7,
    HITRAN_CLASS8,
    HITRAN_CLASS9,
    HITRAN_CLASS10,
)
from radis.io.hitran import columns_2004 as hitrancolumns
from radis.io.cdsd import columns_hitemp as cdsdcolumns
from radis.io.cdsd import columns_4000 as cdsd4000columns
from radis.phys.convert import cm2nm
from radis.phys.air import vacuum2air
from radis.phys.constants import k_b
from radis.misc.basics import is_float
from radis.misc.utils import NotInstalled
from warnings import warn

try:
    import plotly.offline as py
    import plotly.graph_objs as go
except ImportError:
    py = NotInstalled(
        name="plotly.offline",
        info="LineSurvey requires Plotly. Please install it manually",
    )
    go = NotInstalled(
        name="plotly.graph_objs",
        info="LineSurvey requires Plotly. Please install it manually",
    )


# dictionary to convert branch in numeric (-1, 0, 1) or (P, Q, R) format to (P, Q, R) format
_fix_branch_format = {
    -1: "P",
    "P": "P",
    0: "Q",
    "Q": "Q",
    1: "R",
    "R": "R",
}


def LineSurvey(
    spec,
    overlay=None,
    wunit="cm-1",
    Iunit="hitran",
    medium="air",
    cutoff=None,
    plot="S",
    lineinfo=["int", "A", "El"],
    barwidth=0.07,
    yscale="log",
    display=True,
    filename="line-survey.html",
    xunit=None,
    yunit=None,  # deprecated
):
    """ Plot Line Survey (all linestrengths above cutoff criteria) in Plotly (html)


    Parameters    
    ----------

    spec: Spectrum 
        result from SpectrumFactory calculation (see spectrum.py)

    overlay: tuple (w, I, [name], [units]), or list or tuples
        plot (w, I) on a secondary axis. Useful to compare linestrength with
        calculated / measured data 

    wunit: ``'nm'``, ``'cm-1'``
        wavelength / wavenumber units

    Iunit: ``'hitran'``, ``'splot'`` 
        Linestrength output units:

        - ``'hitran'``: (cm-1/(molecule/cm-2))
        - ``'splot'`` : (cm-1/atm)   (Spectraplot units [2]_)

        Note: if not None, cutoff criteria is applied in this unit.
        Not used if plot is not 'S'

    medium: ``'air'``, ``'vacuum'``
        show wavelength either in air or vacuum. Default ``'air'``

    plot: str
        what to plot. Default ``'S'`` (scaled line intensity). But it can be 
        any key in the lines, such as population (``'nu'``), or Einstein coefficient (``'Aul'``)

    lineinfo: list, or ``'all'``
        extra line information to plot. Should be a column name in the databank
        (s.lines). For instance: ``'int'``, ``'selbrd'``, etc... Default [``'int'``]

    Other Parameters
    ----------------

    display: boolean
        if True, open the image in a browser. Default ``True``. 

    filename: str
        filename to save .html 

    yscale: ``'log'``, ``'linear'``
        Default ``'log'``

    Returns
    -------

    Plot in Plotly. See Output in [1]_


    Examples
    --------

    An example using the :class:`~radis.lbl.factory.SpectrumFactory` to generate a spectrum::

        from radis import SpectrumFactory
        sf = SpectrumFactory(
                             wavenum_min=2380,
                             wavenum_max=2400,
                             mole_fraction=400e-6,
                             path_length=100,  # cm
                             isotope=[1],
                             db_use_cached=True) 
        sf.load_databank('HITRAN-CO2-TEST')
        s = sf.eq_spectrum(Tgas=1500)
        s.apply_slit(0.5)
        s.line_survey(overlay='radiance_noslit', barwidth=0.01)
        
    See the output in :ref:`Examples <label_examples>`
    
    References
    ----------

    .. [1] `RADIS Online Documentation (LineSurvey) <https://radis.readthedocs.io/en/latest/tools/line_survey.html>`__

    .. [2] `SpectraPlot <http://www.spectraplot.com/survey>`__


    """

    # Check inputs
    if xunit is not None:
        warn(DeprecationWarning("xunit replaced with wunit"))
        wunit = xunit
    if yunit is not None:
        warn(DeprecationWarning("yunit replaced with Iunit"))
        Iunit = yunit
    assert yscale in ["log", "linear"]

    try:
        spec.lines
        assert spec.lines is not None
    except (AttributeError, AssertionError):
        raise AttributeError(
            "Spectrum has no `lines`. Cant run line survey. If your code allows "
            + "it, recompute the spectrum with 'export_lines=True'"
        )

    # Get input
    T = spec.conditions["Tgas"]
    P = spec.conditions["pressure_mbar"] / 1000  # bar1013.25 # atm
    Xi = spec.conditions["mole_fraction"]
    sp = spec.lines.copy()
    dbformat = spec.conditions["dbformat"]

    if not plot in list(sp.keys()):
        raise KeyError(
            "Key {0} is not in line database: {1}".format(plot, list(sp.keys()))
        )

    def hitran2splot(S):
        """ convert Linestrength in HITRAN units (cm-1/(molecules.cm-2)) to
        SpectraPlot units (cm-2/atm) """

        return S / (k_b * T) / 10

    # Parsers to get units and more details
    if dbformat == "hitran":
        columndescriptor = hitrancolumns
    elif dbformat == "cdsd-hitemp":
        columndescriptor = cdsdcolumns
    elif dbformat == "cdsd-4000":
        columndescriptor = cdsd4000columns

    # Apply cutoff, get ylabel
    if plot == "S":
        if Iunit == "hitran":
            if cutoff is None:
                cutoff = spec.conditions["cutoff"]
            sp = sp[(sp.S > cutoff)]
            Iunit_str = "cm-1/(molecule/cm-2)"
        elif Iunit == "splot":
            if cutoff is None:
                cutoff = spec.conditions["cutoff"]
                # if None, use default cutoff expressed in Hitran units
                sp = sp[(sp.S > cutoff)]
            else:
                sp["S"] = hitran2splot(sp.S)
                sp = sp[(sp.S > cutoff)]
            Iunit_str = "cm-1/atm"
        else:
            raise ValueError("Unknown Iunit: {0}".format(Iunit))
        ylabel = "Linestrength ({0})".format(Iunit_str)
    else:
        cutoff = 0
        try:  # to find units and real name (if exists in initial databank)
            _, _, name, unit = columndescriptor[plot]
            ylabel = "{0} - {1} [{2}]".format(name, plot, unit)
        except KeyError:
            ylabel = plot

    print(len(sp), "lines")
    if len(sp) == 0:
        raise ValueError("0 lines. Change your cutoff criteria?")

    # %% Plot - Plotly version

    def get_x(w):
        """ w (input) is supposed to be in vacuum wavenumbers in the 
        lines database """

        # Convert wavelength / wavenumber
        if wunit == "cm-1":
            x = w
        elif wunit == "nm":
            x = cm2nm(w)

            # Correct if requested medium is air
            if medium == "air":
                x = vacuum2air(x)
            elif medium == "vacuum":
                pass
            else:
                raise ValueError("Unknown medium: {0}".format(medium))
        else:
            raise ValueError("Unknown wunit: {0}".format(wunit))
        return x

    # Parse databank to get relevant information on each line
    # (one function per databank format)

    def get_label_hitran(row, details):
        """ 
        Todo
        -------

        replace with simple astype(str) statements and str operations

        ex:        
        > '['+df[locl].astype(str)+']('+df[globl].astype(str)+'->'+
        >     df[globu].astype(str)'+)'

        will be much faster! 
        """

        molecule = get_molecule(row.id)

        # Get global labels
        if molecule in HITRAN_CLASS1:
            label = "{molec}[iso{iso:.0f}] [{branch}{jl:.0f}]({vl:.0f})->({vu:.0f})".format(
                **dict(
                    [(k, row[k]) for k in ["vu", "vl", "jl", "iso"]]
                    + [
                        ("molec", molecule),
                        ("branch", _fix_branch_format[row["branch"]]),
                    ]
                )
            )
        elif molecule in HITRAN_CLASS4:
            label = "{molec}[iso{iso:.0f}] [{branch}{jl:.0f}]({v1l:.0f}{v2l:.0f}`{l2l:.0f}`{v3l:.0f})->({v1u:.0f}{v2u:.0f}`{l2u:.0f}`{v3u:.0f})".format(
                **dict(
                    [
                        (k, row[k])
                        for k in [
                            "v1u",
                            "v2u",
                            "l2u",
                            "v3u",
                            "v1l",
                            "v2l",
                            "l2l",
                            "v3l",
                            "jl",
                            "iso",
                        ]
                    ]
                    + [
                        ("molec", molecule),
                        ("branch", _fix_branch_format[row["branch"]]),
                    ]
                )
            )
        elif molecule in HITRAN_CLASS5:
            label = "{molec}[iso{iso:.0f}] [{branch}{jl:.0f}]({v1l:.0f}{v2l:.0f}`{l2l:.0f}`{v3l:.0f} {rl:.0f})->({v1u:.0f}{v2u:.0f}`{l2u:.0f}`{v3u:.0f} {ru:.0f})".format(
                **dict(
                    [
                        (k, row[k])
                        for k in [
                            "v1u",
                            "v2u",
                            "l2u",
                            "v3u",
                            "v1l",
                            "v2l",
                            "l2l",
                            "v3l",
                            "rl",
                            "ru",
                            "jl",
                            "iso",
                        ]
                    ]
                    + [
                        ("molec", molecule),
                        ("branch", _fix_branch_format[row["branch"]]),
                    ]
                )
            )
        else:
            raise NotImplementedError(
                "No label for {0}. Please add it!".format(molecule)
            )

        # Add details about some line properties
        for k in details:
            name, _, unit = details[k]
            if is_float(row[k]):
                label += "<br>{0} {1}: {2:.3g} {3}".format(k, name, row[k], unit)
            else:
                label += "<br>{0} {1}: {2} {3}".format(k, name, row[k], unit)

        return label

    def get_label_cdsd(row, details):
        label = "CO2[iso{iso}] [{branch}{jl:.0f}](p{polyl:.0f}c{wangl:.0f}n{rankl:.0f})->(p{polyu:.0f}c{wangu:.0f}n{ranku:.0f})".format(
            **dict(
                [
                    (k, row[k])
                    for k in [
                        "polyl",
                        "wangl",
                        "rankl",
                        "polyu",
                        "wangu",
                        "ranku",
                        "jl",
                        "iso",
                    ]
                ]
                + [("branch", _fix_branch_format[row["branch"]])]
            )
        )

        for k in details:
            name, _, unit = details[k]
            if is_float(row[k]):
                label += "<br>{0} {1}: {2:.3g} {3}".format(k, name, row[k], unit)
            else:
                label += "<br>{0} {1}: {2} {3}".format(name, k, row[k], unit)

        return label

    def get_label_cdsd_hitran(row, details):
        label = "CO2[iso{iso}] [{branch}{jl:.0f}]({v1l:.0f}{v2l:.0f}`{l2l:.0f}`{v3l:.0f})->({v1u:.0f}{v2u:.0f}`{l2u:.0f}`{v3u:.0f})".format(
            **dict(
                [
                    (k, row[k])
                    for k in [
                        "v1u",
                        "v2u",
                        "l2u",
                        "v3u",
                        "v1l",
                        "v2l",
                        "l2l",
                        "v3l",
                        "jl",
                        "iso",
                    ]
                ]
                + [("branch", _fix_branch_format[row["branch"]])]
            )
        )

        for k in details:
            name, _, unit = details[k]
            if is_float(row[k]):
                label += "<br>{0} {1}: {2:.3g} {3}".format(k, name, row[k], unit)
            else:
                label += "<br>{0} {1}: {2} {3}".format(name, k, row[k], unit)

        return label

    def get_label_none(row):
        return "unknown databank format. \ndetails cant be read"

    # add extra info to label
    details = {}
    if columndescriptor:
        for k in lineinfo:
            if not k in sp.columns:
                raise KeyError(
                    "{0} not a {1} databank entry ({2} format)".format(
                        k, columndescriptor, dbformat.upper()
                    )
                )
            try:  # to find units and real name (if exists in initial databank)
                _, ktype, name, unit = columndescriptor[k]
                details[k] = (name, ktype, " [{0}]".format(unit))
            except:
                details[k] = ("", None, "")  # keep short name

    # Get label
    if dbformat == "hitran":
        sp["label"] = sp.apply(lambda r: get_label_hitran(r, details), axis=1)
    elif dbformat in ["cdsd-hitemp", "cdsd-4000"]:
        try:
            sp["label"] = sp.apply(lambda r: get_label_cdsd_hitran(r, details), axis=1)
        except KeyError:
            sp["label"] = sp.apply(lambda r: get_label_cdsd(r, details), axis=1)
    else:
        sp["label"] = sp.apply(get_label_none, axis=1)

    # from plotly.graph_objs import Scatter, Figure, Layout
    #
    l = [
        go.Bar(
            x=get_x(sp.shiftwav),
            y=sp[plot],
            text=sp.label,
            width=barwidth,
            name="linestrength",
        )
    ]

    if wunit == "nm":
        xlabel = "Wavelength (nm) [{0}]".format(medium)
    elif wunit == "cm-1":
        xlabel = "Wavenumber (cm-1)"
    else:
        raise ValueError("unknown wunit: {0}".format(wunit))

    if yscale == "log":
        plot_range = (
            int(np.log10(max(cutoff, sp[plot].min()))),
            int(np.log10(sp[plot].max())) + 1,
        )
    else:
        plot_range = (max(cutoff, sp[plot].min()), sp[plot].max() + 1)

    layout = go.Layout(
        title="Line Survey ({T}K, {P:.3f}bar, Mfrac={Xi:.3f})".format(
            **{"T": T, "P": P, "Xi": Xi}
        ),
        hovermode="closest",
        xaxis=dict(title=xlabel,),
        yaxis=dict(
            title=ylabel,
            # note: LaTeX doesnt seem to work in Offline mode yet.
            type=yscale,
            range=plot_range,
            titlefont=dict(color="#1f77b4"),
            tickfont=dict(color="#1f77b4"),
        ),
        showlegend=False,
    )

    # %% Add overlay
    def add_overlay(overlay):

        over_w = overlay[0]
        over_I = overlay[1]
        args = overlay[2:]
        over_name = args[0] if len(args) > 0 else "overlay"
        over_units = args[1] if len(args) > 1 else "a.u"

        l.append(
            go.Scatter(
                x=over_w,
                #                        y=spec.out[overlay],
                y=over_I,
                yaxis="y2",
                name=over_name,
            )
        )

        layout["yaxis2"] = dict(
            title="{0} ({1})".format(
                over_name.capitalize(), over_units
            ),  # note: LaTeX doesnt seem to
            # work in Offline mode yet.
            overlaying="y",
            type="log",
            side="right",
            titlefont=dict(color="#ff7f0e"),
            tickfont=dict(color="#ff7f0e"),
        )

    if overlay is not None:

        if type(overlay) is not list:
            overlay = [overlay]

        for ov in overlay:
            add_overlay(ov)

    # %% Plot

    fig = go.Figure(data=l, layout=layout)
    py.plot(fig, filename=filename, auto_open=display)


# %% Test
if __name__ == "__main__":

    from radis.test.tools.test_line_survey import _run_testcases

    print("Testing line survey functions:", _run_testcases(plot=True))
