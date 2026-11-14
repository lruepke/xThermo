"""Gallery: the H2O-NaCl phase diagram at four pressures (depth_sections.png).

Writes to docs/assets/img, or to the directory given as the first argument or in $XTHERMO_FIGURE_DIR.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# --8<-- [start:figure]
import matplotlib.pyplot as plt
import xthermo
from xthermo import H2ONaCl

import _style as st                              # fonts and colours of the documentation figures (docs/gallery)

st.apply()
sw = H2ONaCl.cH2ONaCl("IAPS84")
pressures = [100e5, 250e5, 400e5, 1000e5]        # Pa

fig, axes = plt.subplots(1, 4, figsize=(13, 4.4), sharey=True, layout="constrained")
handles = {}
for ax, p in zip(axes, pressures):
    ax.set_xscale("log")
    h = xthermo.plot_slice(ax, sw.Slice_constP(p), x="X", y="T", points=False)
    handles.update(h)
    ax.set_xlim(1e-3, 100)                       # wt% NaCl
    ax.set_ylim(0, 1000)                         # °C
    ax.set_title(f"{p / 1e5:.0f} bar", loc="center")
    ax.axvline(3.2, color=st.CORAL, lw=1.2, ls=(0, (5, 3)))   # seawater
    if ax is not axes[0]:
        ax.set_ylabel("")
# one legend entry for the three-phase assemblage (its temperature differs between the panels)
vlh = [k for k in handles if k.startswith("V+L+H")]
handles = {k: v for k, v in handles.items() if k not in vlh} | {"V+L+H (a line in T–X)": handles[vlh[0]]}
fig.legend(handles.values(), handles.keys(), loc="outside lower center", ncol=7, fontsize="small")
fig.suptitle("The H₂O–NaCl phase diagram at four pressures (dashed: seawater, 3.2 wt% NaCl)", x=0.01, ha="left",
             color=st.NAVY, fontweight="bold")
plt.show()
# --8<-- [end:figure]
fig.savefig(os.path.join(st.outdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img")),
                         "depth_sections.png"))
print("depth_sections.png written")
