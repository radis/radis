# -*- coding: utf-8 -*-
"""
Functions to plot line surveys


-------------------------------------------------------------------------------


"""


from warnings import warn

import numpy as np

from radis.api.cdsdapi import columns_4000 as cdsd4000columns
from radis.api.cdsdapi import columns_cdsdhitemp as cdsdhitempcolumns
from radis.api.hitranapi import (  # HITRAN_CLASS2,; HITRAN_CLASS3,; HITRAN_CLASS6,; HITRAN_CLASS7,; HITRAN_CLASS8,; HITRAN_CLASS9,; HITRAN_CLASS10,
    HITRAN_CLASS1,
    HITRAN_CLASS4,
    HITRAN_CLASS5,
)
from radis.api.hitranapi import columns_2004 as hitrancolumns
from radis.db.classes import get_molecule, get_molecule_identifier
from radis.misc.basics import is_float
from radis.misc.printer import superscript_map
from radis.misc.utils import NotInstalled
from radis.phys.air import vacuum2air
from radis.phys.constants import k_b
from radis.phys.convert import cm2nm

try:
    import plotly.graph_objs as go
    import plotly.offline as py
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
    barwidth="hwhm_voigt",  # 0.01,
    yscale="log",
    writefile=None,
    xunit=None,
    yunit=None,  # deprecated
):
    """Plot Line Survey (all linestrengths above cutoff criteria) in Plotly
    (html)

    Parameters
    ----------
    spec: Spectrum
        result from SpectrumFactory calculation (see spectrum.py)
    overlay: tuple (w, I, [name], [units]), or list or tuples
        plot (w, I) on a secondary axis. Useful to compare linestrength with
        calculated / measured data::

            LineSurvey(overlay='abscoeff')
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
        any key in the lines, such as population (``'nu'``), or Einstein coefficient (``'A'``)
    lineinfo: list, or ``'all'``
        extra line information to plot. Should be a column name in the databank
        (s.lines). For instance: ``'int'``, ``'selbrd'``, etc... Default [``'int'``]

    Other Parameters
    ----------------
    writefile: str
        if not ``None``, a valid filename to save the plot under .html format.
        If ``None``, use the ``fig`` object returned to show the plot.
    yscale: ``'log'``, ``'linear'``
        Default ``'log'``
    barwidth: float or str
        if float, width of bars, in ``wunit``, as a fraction of full-range; i.e. ::

            barwidth=0.01

        makes bars span 1% of the full range.
        if ``str``, uses the column as width. Example ::

            barwidth = 'hwhm_voigt'

    Returns
    -------
    fig: a Plotly figure.
        If using a Jupyter Notebook, the plot will appear. Else, use ``writefile``
        to export to an html file.

    See typical output in [1]_

    Examples
    --------
    An example using the :class:`~radis.lbl.factory.SpectrumFactory` to generate a spectrum,
    using the Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey` method
    directly ::

        from radis import SpectrumFactory
        sf = SpectrumFactory(
                             wavenum_min=2380,
                             wavenum_max=2400,
                             mole_fraction=400e-6,
                             path_length=100,  # cm
                             isotope=[1],
                             export_lines=True,    # required for LineSurvey!
                             db_use_cached=True)
        sf.load_databank('HITRAN-CO2-TEST')
        s = sf.eq_spectrum(Tgas=1500)
        s.apply_slit(0.5)
        s.line_survey(overlay='radiance_noslit', barwidth=0.01, lineinfo="all")

    See the output in :ref:`Examples <label_examples>`

    .. raw:: html

        <iframe id="igraph" src="//plotly.com/~erwanp/6.embed" width="650" height="420" seamless="seamless" scrolling="no"></iframe>

    .. minigallery:: radis.spectrum.spectrum.Spectrum.line_survey
        :add-heading: Examples using Line Survey

    References
    ----------

    .. [1] `RADIS Online Documentation (LineSurvey) <https://radis.readthedocs.io/en/latest/tools/line_survey.html>`__

    .. [2] `SpectraPlot <http://www.spectraplot.com/survey>`__

    See Also
    --------
    :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`
    """

    # Check inputs
    if xunit is not None:
        warn(DeprecationWarning("xunit replaced with wunit"))
        wunit = xunit
    if yunit is not None:
        warn(DeprecationWarning("yunit replaced with Iunit"))
        Iunit = yunit
    assert yscale in ["log", "linear"]

    # auto plot in Spyder IDE
    from os import environ

    if writefile is None and any("SPYDER" in name for name in environ):
        writefile = "line_survey.html"

    try:
        spec.lines
        assert spec.lines is not None
    except (AttributeError, AssertionError):
        raise AttributeError(
            "Spectrum has no `lines`. Cant run line survey. If your code allows "
            + "it, recompute the spectrum with 'export_lines=True'"
        )

    # Get input
    # T = spec.conditions["Tgas"]
    # P = spec.conditions["pressure"]  # bar1.01325 # atm
    # Xi = spec.conditions["mole_fraction"]
    sp = spec.lines.copy()

    if not plot in list(sp.keys()):
        raise KeyError(
            f"Key {plot} is not in line database: {list(sp.keys())}. Change `plot=` parameter of line_survey"
        )

    def hitran2splot(S, T):
        r"""convert Linestrength in HITRAN units (cm-1/(molecules.cm-2)) to
        SpectraPlot units (cm-2/atm)"""

        return S / (k_b * T) / 10

    # Known parsers to get units and more details
    dbformat = None
    if "dbformat" in spec.conditions:
        dbformat = spec.conditions["dbformat"]
        if dbformat in ["hitran", "hitemp", "hitemp-radisdb", "geisa"]:
            columndescriptor = hitrancolumns
        elif dbformat in ["nist"]:
            columndescriptor = hitrancolumns  # Placeholder. TODO replace with NIST (for the moment we copy the names anyway)
        elif dbformat == "cdsd-hitemp":
            columndescriptor = cdsdhitempcolumns
        elif dbformat == "cdsd-4000":
            columndescriptor = cdsd4000columns
        else:
            warn(f"unknown dbformat {dbformat}")
            columndescriptor = {}
    else:
        columndescriptor = {}

    def _fmt_val(k, v):
        """Format a value: wavenumber with 2 decimals, others with 3 significant figures."""
        if not is_float(v):
            return f"{v}"
        if k == "wav":
            return f"{v:.2f}"
        return f"{v:.3g}"

    def _superscript_unit(unit):
        """Convert plain-text unit exponents to Unicode superscripts.

        e.g. 'cm-1/(molecule/cm-2)' -> 'cm⁻¹/(molecule/cm⁻²)'
        """
        repl = {"-1": "\u207b\u00b9", "-2": "\u207b\u00b2", "-3": "\u207b\u00b3"}
        for old, new in repl.items():
            unit = unit.replace(old, new)
        return unit

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
                Tgas = spec.conditions["Tgas"]
                sp["S"] = hitran2splot(sp.S, Tgas)
                sp = sp[(sp.S > cutoff)]
            Iunit_str = "cm-1/atm"
        else:
            raise ValueError(f"Unknown Iunit: {Iunit}")
        ylabel = f"Linestrength ({_superscript_unit(Iunit_str)})"
    else:
        cutoff = 0
        try:  # to find units and real name (if columndescriptor is known exists in initial databank)
            _, _, name, unit = columndescriptor[plot]
            ylabel = f"{name} - {plot} ({_superscript_unit(unit)})"
        except KeyError:
            ylabel = plot

    print(len(sp), "lines")
    if len(sp) == 0:
        raise ValueError("0 lines. Change your cutoff criteria?")

    # %% Plot - Plotly version

    def get_x(w):
        r"""w (input) is supposed to be in vacuum wavenumbers in the lines
        database."""

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
                raise ValueError(f"Unknown medium: {medium}")
        else:
            raise ValueError(f"Unknown wunit: {wunit}")
        return x

    # Parse databank to get relevant information on each line
    # (one function per databank format)

    # Display name mapping for common spectroscopic quantities
    _display_names = {
        "int": "S",
        "El": "E\u2032\u2032",
        "wav": "\U0001d708",
    }

    # Units for common columns (used for unknown/ExoMol formats)
    _display_units = {
        "wav": "cm\u207b\u00b9",
        "int": "cm\u207b\u00b9/(molecule/cm\u207b\u00b2)",
        "El": "cm\u207b\u00b9",
        "A": "s\u207b\u00b9",
        "airbrd": "cm\u207b\u00b9/atm",
        "selbrd": "cm\u207b\u00b9/atm",
        "S": "cm\u207b\u00b9/(molecule/cm\u207b\u00b2)",
    }

    def add_details(row, details):
        r"""Add details in string; add on 2 columns if "upper" and "lower" value follow"

        Example :  keep "gu", "gl" on the same line

        Also add unit if known
        """
        label = ""
        for k in details:
            _, _, unit = details[k]
            if k[-1] == "u" and k[:-1] + "l" in details:
                continue  #  don't show "upper" value. Will be shown on same line as lower

            display = _display_names.get(k, k)
            val = _fmt_val(k, row[k])
            label += f"<br><b>{display}</b>: {val} {unit}"
            # If lower value, also show upper value on same line
            if k[-1] == "l" and k[:-1] + "u" in details:
                k_up = k[:-1] + "u"
                display_up = _display_names.get(k_up, k_up)
                val_up = _fmt_val(k_up, row[k_up])
                label += (
                    "&nbsp;" * (30 - len(display) - len(val))
                    + f"<b>{display_up}</b>: {val_up} {unit}"
                )

        return label

    def get_label_hitran(row, details, attrs):
        r"""
        Todo
        -------

        replace with simple astype(str) statements and str operations

        ex:
        > '['+df[locl].astype(str)+']('+df[globl].astype(str)+'->'+
        >     df[globu].astype(str)'+)'

        will be much faster!
        """

        if "id" in attrs:
            molecule = get_molecule(attrs["id"])

        else:
            molecule = spec.conditions["molecule"]
            if isinstance(molecule, str):
                id = get_molecule_identifier(molecule)
                molecule = get_molecule(id)

        try:
            # Try customized input

            # Get global labels
            if (
                molecule in HITRAN_CLASS1
            ):  # ["CO", "HF", "HCl", "HBr", "HI", "N2", "NO+"]

                add = ["vu", "vl", "jl"]

                if "iso" in row:
                    iso = int(row["iso"])
                elif "isotope" in spec.conditions:
                    iso = int(spec.conditions["isotope"])
                else:
                    iso = "?"

                label = "{molec}[iso{iso}] {branch}({vl:.0f},{jl:.0f})".format(
                    **dict(
                        [(k, row[k]) for k in add]
                        + [
                            ("iso", iso),
                            ("molec", molecule),
                            ("branch", _fix_branch_format[row["branch"]]),
                        ]
                    )
                )
            elif molecule in HITRAN_CLASS4:  #  ["N2O", "OCS", "HCN"]

                add = [
                    "v1u",
                    "v2u",
                    "l2u",
                    "v3u",
                    "v1l",
                    "v2l",
                    "l2l",
                    "v3l",
                    "jl",
                ]

                if "iso" in row:
                    iso = int(row["iso"])
                elif "isotope" in spec.conditions:
                    iso = int(spec.conditions["isotope"])
                else:
                    iso = "?"

                label = "{molec}[iso{iso}] {branch}({v1l:.0f}{v2l:.0f}`{l2l:.0f}`{v3l:.0f},{jl:.0f})->({v1u:.0f}{v2u:.0f}`{l2u:.0f}`{v3u:.0f})".format(
                    **dict(
                        [(k, row[k]) for k in add]
                        + [
                            ("iso", iso),
                            ("molec", molecule),
                            ("branch", _fix_branch_format[row["branch"]]),
                        ]
                    )
                )
            elif molecule in HITRAN_CLASS5:  # ["CO2"]

                add = [
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
                ]

                if "iso" in row:
                    iso = int(row["iso"])
                elif "isotope" in spec.conditions:
                    iso = int(spec.conditions["isotope"])
                else:
                    iso = "?"

                label = "{molec}[iso{iso}] {branch}({v1l:.0f}{v2l:.0f}{l2l}{v3l:.0f},{jl:.0f} r={rl})->({v1u:.0f}{v2u:.0f}{l2u}{v3u:.0f} r={ru})"
                label = label.format(
                    **dict(
                        [
                            (k, row[k])
                            for k in add
                            if k not in ["l2l", "l2u", "ru", "rl"]
                        ]
                        +
                        # Try to get l number as UTF-8 superscript (works if in superscript_map)
                        [
                            (
                                k,
                                superscript_map.get(
                                    str(int(row[k])), f"`{row[k]:.0f}`"
                                ),
                            )
                            for k in add
                            if k in ["l2l", "l2u"]
                        ]  # ignore missing params (like rl, ru ? )
                        + [
                            ("iso", iso),
                            ("molec", molecule),
                            ("branch", _fix_branch_format[row["branch"]]),
                            # 'ru', 'rl' may be missing if not all columns loaded
                            ("ru", f"{row['ru']:.0f}" if "ru" in row else "?"),
                            ("rl", f"{row['rl']:.0f}" if "rl" in row else "?"),
                        ]
                    )
                )

            else:
                raise NotImplementedError(
                    f"No customized label for {molecule}. Please add it!"
                )

            # Add details about some line properties
            label += add_details(row, details)

        except KeyError as err:
            print(
                f"Error during customized Line survey labelling. Printing everything : \n{str(err)}"
            )

            # Add everything we know
            # raise
            label = row.to_string().replace("\n", "<br>")
            if "id" not in row:
                label = "{molecule}<br>" + label
            else:
                label = get_molecule(row["id"]) + "<br>" + label

        return label

    def get_label_cdsd(row, details):

        add = [
            "polyl",
            "wangl",
            "rankl",
            "polyu",
            "wangu",
            "ranku",
            "jl",
        ]

        if "iso" in row:
            iso = int(row["iso"])
        elif "isotope" in spec.conditions:
            iso = int(spec.conditions["isotope"])
        else:
            iso = "?"

        label = "CO2[iso{iso}] {branch}(p{polyl:.0f}c{wangl:.0f}n{rankl:.0f},{jl:.0f})->(p{polyu:.0f}c{wangu:.0f}n{ranku:.0f})".format(
            **dict(
                [(k, row[k]) for k in add]
                + [("iso", iso), ("branch", _fix_branch_format[row["branch"]])]
            )
        )

        label += add_details(row, details)

        return label

    def get_label_cdsd_hitran(row, details):
        if "iso" in row:
            iso = int(row["iso"])
        elif "isotope" in spec.conditions:
            iso = int(spec.conditions["isotope"])
        else:
            iso = "?"

        label = "CO2[iso{iso}] {branch}({v1l:.0f}{v2l:.0f}`{l2l:.0f}`{v3l:.0f},{jl:.0f})->({v1u:.0f}{v2u:.0f}`{l2u:.0f}`{v3u:.0f})".format(
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
                    ]
                ]
                + [("iso", iso), ("branch", _fix_branch_format[row["branch"]])]
            )
        )

        label += add_details(row, details)

        return label

    def get_label_nist(row, attrs):
        label = attrs[
            "molecule"
        ] + " [{Lower level}] ({El:.2f} eV) -> [{Upper level}] ({Eu:.2f} eV)".format(
            **dict(
                [
                    (k, row[k])
                    for k in [
                        "Lower level",  # TODO: have nicer Term Symbol appear
                        "Upper level",  # TODO: have nicer Term Symbol appear
                        "El",
                        "Eu",
                    ]
                ]
            )
        )

        return label

    def get_label_all(row, columns=None):
        r"""Print line details, optionally filtered to certain columns.

        Parameters
        ----------
        row : pandas.Series
            Line data row
        columns : list, optional
            If provided, only show these columns. Otherwise show all.
        """
        # Build a spectroscopic header if branch/vl/jl columns are available
        header = ""
        if "branch" in row and "jl" in row:
            try:
                branch = _fix_branch_format.get(row["branch"], row["branch"])
                if "vl" in row:
                    header = f"{branch}({row['vl']:.0f},{row['jl']:.0f})"
                else:
                    header = f"{branch}({row['jl']:.0f})"
            except (ValueError, TypeError):
                pass
        if header:
            label = header + "<br>"
        else:
            label = ""
        if columns is not None:
            items = [(k, row[k]) for k in columns if k in row.index]
        else:
            items = row.items()
        label += "<br>".join(
            [
                f"<b>{_display_names.get(k, k)}</b>: {_fmt_val(k, v)}"
                + (f" {_display_units[k]}" if k in _display_units else "")
                for k, v in items
            ]
        )
        return label

    def get_label_none(row):
        return "unknown databank format. \ndetails cant be read"

    # add extra info to label
    if lineinfo == "all":
        lineinfo = sp.columns
    details = {}
    if columndescriptor:
        for k in lineinfo:
            if not k in sp.columns:
                raise KeyError(
                    f"{k} not a {columndescriptor} databank entry ({dbformat.upper()} format)"
                )
            try:  # to find units and real name (if exists in initial databank)
                _, ktype, name, unit = columndescriptor[k]
            except KeyError:
                fallback_unit = _display_units.get(k, "")
                details[k] = ("", None, f" {fallback_unit}" if fallback_unit else "")
            else:
                details[k] = (name, ktype, f" {_superscript_unit(unit)}")

    # Get label
    if dbformat in ["hitran", "hitemp", "hitemp-radisdb", "radisdb-hitemp", "geisa"]:
        sp["label"] = sp.apply(lambda r: get_label_hitran(r, details, sp.attrs), axis=1)
    elif dbformat in ["cdsd-hitemp", "cdsd-4000"]:
        try:
            sp["label"] = sp.apply(lambda r: get_label_cdsd_hitran(r, details), axis=1)
        except KeyError:
            sp["label"] = sp.apply(lambda r: get_label_cdsd(r, details), axis=1)
    elif dbformat in ["nist"]:
        sp["label"] = sp.apply(lambda r: get_label_nist(r, sp.attrs), axis=1)
    else:
        sp["label"] = sp.apply(lambda r: get_label_all(r, lineinfo), axis=1)

        # sp["label"] = sp.apply(get_label_none, axis=1)

    # from plotly.graph_objs import Scatter, Figure, Layout

    if wunit == "nm":
        xlabel = f"Wavelength (nm) [{medium}]"
        x_range = spec.get_wavelength()
    elif wunit == "cm-1":
        xlabel = "Wavenumber (cm-1)"
        x_range = spec.get_wavenumber()
    else:
        raise ValueError(f"unknown wunit: {wunit}")

    # TODO : implement barwidth -->  if ``hwhm_voigt``, bars are one half-width at half-maximum (in cm-1)
    if isinstance(barwidth, str):
        if not barwidth in sp.columns:
            raise ValueError(
                f"`{barwidth}` not in the lines Dataframe columns. Use one of the column names ({sp.columns}) or set `barwidth=[a 0-1 number]` to plot bars with widths as a fraction of the total range, ex : `barwidth=0.01`"
            )
        barwidth = sp[barwidth]
    else:
        barwidth = (x_range.max() - x_range.min()) * barwidth
    l = [
        go.Bar(
            x=get_x(sp.shiftwav),
            y=sp[plot],
            text=sp.label,
            hoverinfo="text",
            width=barwidth,
            name="linestrength",
        )
    ]

    if yscale == "log":
        plot_range = (
            int(np.log10(max(cutoff, sp[plot].min()))) - 1,
            int(np.log10(sp[plot].max())) + 1,
        )
    else:
        plot_range = (max(cutoff, sp[plot].min()), sp[plot].max() + 1)

    layout = go.Layout(
        title="Line Survey",
        # TODO : re-add some parameters in the title (if they exist)
        # ({T}K, {P:.3f}bar, Mfrac={Xi:.3f})".format(
        #    **{"T": T, "P": P, "Xi": Xi}
        # ),
        hovermode="closest",
        hoverlabel=dict(
            bgcolor="white",
            font=dict(color="#333", size=13, family="Arial, sans-serif"),
            bordercolor="#ccc",
        ),
        xaxis=dict(
            title=xlabel,
            range=(x_range.min(), x_range.max()),
        ),
        yaxis=dict(
            title={"text": ylabel, "font": {"color": "#1f77b4"}},
            # note: LaTeX doesnt seem to work in Offline mode yet.
            type=yscale,
            range=plot_range,
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
            title={
                "text": f"{over_name.capitalize()} ({over_units})",
                "font": {"color": "#ff7f0e"},
            },
            # note: LaTeX doesnt seem to
            # work in Offline mode yet.
            overlaying="y",
            type="log",
            side="right",
            tickfont=dict(color="#ff7f0e"),
        )

    if overlay is not None:

        if type(overlay) is not list:
            overlay = [overlay]

        for ov in overlay:
            add_overlay(ov)

    # %% Plot

    fig = go.Figure(data=l, layout=layout)

    if writefile:
        py.plot(fig, filename=writefile, auto_open=True)

        # Inject copy-to-clipboard JS + hint text into the HTML file
        _clipboard_js = (
            "\n<script>\n"
            "document.addEventListener('DOMContentLoaded', function() {\n"
            "    var plot = document.querySelector('.plotly-graph-div');\n"
            "    if (plot) {\n"
            "        plot.on('plotly_click', function(data) {\n"
            "            var text = data.points[0].text;\n"
            "            if (text) {\n"
            "                var plain = text.replace(/<br>/g, '\\n').replace(/<[^>]*>/g, '');\n"
            "                navigator.clipboard.writeText(plain).then(function() {\n"
            "                    var el = document.createElement('div');\n"
            "                    el.textContent = 'Copied to clipboard!';\n"
            "                    el.style.cssText = 'position:fixed;bottom:20px;left:50%;"
            "transform:translateX(-50%);background:#4CAF50;color:white;"
            "padding:10px 20px;border-radius:5px;z-index:9999;"
            "font-family:sans-serif;';\n"
            "                    document.body.appendChild(el);\n"
            "                    setTimeout(function(){ el.remove(); }, 1500);\n"
            "                });\n"
            "            }\n"
            "        });\n"
            "    }\n"
            "});\n"
            "</script>\n"
            "<p style='text-align:center;color:#888;font-family:sans-serif;"
            "font-size:12px;margin-top:5px;'>"
            "Click on plot to copy data to clipboard.</p>\n"
        )
        # Insert before </html> tag
        with open(writefile, "r", encoding="utf-8") as f:
            html = f.read()
        html = html.replace("</html>", _clipboard_js + "</html>")
        with open(writefile, "w", encoding="utf-8") as f:
            f.write(html)

    return fig


# %% Test
if __name__ == "__main__":

    from radis.test.tools.test_line_survey import _run_testcases

    print("Testing line survey functions:", _run_testcases(plot=True))
