# -*- coding: utf-8 -*-
"""
=========================
Use different plot themes
=========================

Shows how to customize plot styles, using :py:mod:`seaborn` , :py:mod:`matplotlib`
or `publib <https://github.com/erwanp/publib>`__

Examples
--------
To change it in your user script, set the keys of the ``"plot"`` bloc in :py:attr:`radis.config` :

::

    import radis
    radis.config["plot"]["plotlib"] = "seaborn"
    radis.config["plot"]["context"] = "paper"
    radis.config["plot"]["style"] = "darkgrid"

To change your default settings, edit the
``~/radis.json`` :ref:`Configuration file <label_lbl_config_file>`

See Also
--------
:py:func:`~radis.misc.plot.set_style`, :py:func:`~radis.misc.plot.fix_style`,

"""


import matplotlib.pyplot as plt

import radis

plt.close("all")
s = radis.spectrum_test()

for plotlib, context, style in [
    ("matplotlib", "", "default"),
    ("seaborn", "paper", "darkgrid"),
    ("publib", "paper", "origin"),
    ("seaborn", "poster", "darkgrid"),
    ("publib", "poster", "origin"),
]:
    radis.config["plot"]["plotlib"] = plotlib
    radis.config["plot"]["context"] = context
    radis.config["plot"]["style"] = style
    s.plot()
    plt.title(", ".join(f"{v}" for v in radis.config["plot"].values()))
    plt.tight_layout()

#%%
# You can also choose to not show the figure immediately, and instead return
# the Matplotlib or Plotly object. This allows for **advanced customization**
# (e.g., modifying the matplotlib figure).

fig, ax = s.plot("radiance", plotting_library="matplotlib", show=False)

ax.set_xlabel("Updated X-axis Label")
ax.set_ylabel("Updated Y-axis Label")
ax.grid(True, linestyle=":", linewidth=1)
ax.set_facecolor("beige")

fig.gcf().set_size_inches(8.6, 5)  # Resize
fig.gcf().set_dpi(150)  # Change DPI
# Then Show
fig.show()
