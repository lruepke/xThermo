"""Python guide, "Property maps" page: the figures of docs/python/property_maps.md.

The page shows each section between the snippet markers verbatim; everything outside the markers only saves the
figures for the documentation (docs/assets/img, or the directory given as the first argument or in
$XTHERMO_FIGURE_DIR). Writes map_density_300bar.png, map_temperature_hx.png, map_water_tp.png.
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
warnings.filterwarnings("ignore", message=".*non-interactive.*")   # plt.show() under the Agg backend
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style  # noqa: E402  (output directory only; the figures use Matplotlib's defaults, as in the shown code)

OUT = _style.outdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img"))

# --8<-- [start:density]
import matplotlib.pyplot as plt
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.linspace(1, 999, 250) + 273.15              # K
X = np.logspace(-5, 0, 250)                        # kg/kg
s = xthermo.tpx(sw, T[:, None], 300e5, X[None, :], errors="nan")   # a 250 x 250 grid at 300 bar

fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
ax.set_xscale("log")
xthermo.plot_property(ax, s, "Rho", x="X", y="T", cmap="viridis", vmin=0, vmax=1300,   # cap the colour scale
                      contours=np.arange(200, 1300, 200))
xthermo.plot_slice(ax, sw.Slice_constP(300e5), x="X", y="T", regions=False, points=False)  # boundaries on top
ax.set_title("Bulk density at 300 bar")
plt.show()
# --8<-- [end:density]
fig.savefig(os.path.join(OUT, "map_density_300bar.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# --8<-- [start:temperature]
import matplotlib.pyplot as plt
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
H = np.linspace(0.1e6, 4.4e6, 200)                 # J/kg
X = np.logspace(-4, 0, 200)                        # kg/kg
s = xthermo.hpx(sw, H[:, None], 250e5, X[None, :], errors="nan")

fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
ax.set_xscale("log")
xthermo.plot_property(ax, s, "T", x="X", y=("H", H), cmap="magma",       # H: the input values (see below)
                      contours=np.arange(100, 1000, 100) + 273.15, label="temperature (K)")
xthermo.plot_slice(ax, sw.Slice_constP(250e5), x="X", y="H", regions=False, points=False)
ax.set_xlim(1e-2, 100)                             # wt% NaCl: the grid (boundary lines reach X = 0)
ax.set_ylim(0, 4.6)                                # MJ/kg
ax.set_title("Temperature from H, p, X at 250 bar (grey: not evaluated)")
plt.show()
# --8<-- [end:temperature]
fig.savefig(os.path.join(OUT, "map_temperature_hx.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# --8<-- [start:water]
import matplotlib.pyplot as plt
import numpy as np
import xthermo
from matplotlib.colors import LogNorm
from xthermo import H2O

water = H2O.cIAPS84()
T = np.linspace(280, 1000, 250)                    # K
p = np.geomspace(1e5, 1000e5, 250)                 # Pa
s = xthermo.tpx(water, T[None, :], p[:, None], errors="nan")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
for ax in (ax1, ax2):
    ax.set_yscale("log")
xthermo.plot_phases(ax1, s, x="T", y="P")
xthermo.plot_property(ax2, s, "Rho", x="T", y="P", norm=LogNorm(),
                      contours=[1, 3, 10, 30, 100, 300, 600, 800, 1000])
ax1.set_title("Phase (IAPS-84 water)")
ax2.set_title("Density")
plt.show()
# --8<-- [end:water]
fig.savefig(os.path.join(OUT, "map_water_tp.png"), dpi=150, bbox_inches="tight")
plt.close("all")
print("map_density_300bar.png map_temperature_hx.png map_water_tp.png written to", os.path.abspath(OUT))
