"""Gallery: seawater rising and decompressing at constant specific enthalpy (upflow_seawater.png).

Writes to docs/assets/img, or to the directory given as the first argument or in $XTHERMO_FIGURE_DIR.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# --8<-- [start:figure]
import matplotlib.pyplot as plt
import numpy as np
import xthermo
from xthermo import H2ONaCl

import _style as st                              # fonts and colours of the documentation figures (docs/gallery)

st.apply()
sw = H2ONaCl.cH2ONaCl("IAPS84")
X = 0.032                                        # seawater salinity, kg/kg
p = np.geomspace(500e5, 50e5, 300)               # decompression from 500 to 50 bar, Pa

# where seawater is a single phase or vapour + liquid, in temperature and pressure
T_grid = np.linspace(200, 500, 301) + 273.15
P_grid = np.geomspace(40e5, 550e5, 300)
grid = xthermo.tpx(sw, T_grid[None, :], P_grid[:, None], X, errors="nan")

fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
ax.set_yscale("log")
xthermo.plot_phases(ax, grid, x="T", y="P", legend=False)
ax.invert_yaxis()                                # pressure decreases upwards, as on the way to the seafloor

for T0, colour in [(380, st.NAVY), (420, st.TEAL), (450, st.CORAL)]:
    start = xthermo.tpx(sw, T0 + 273.15, p[0], X)            # the fluid at depth
    path = xthermo.hpx(sw, start.H, p, X, errors="nan")       # same specific enthalpy at every pressure
    ax.plot(path.T - 273.15, p / 1e5, color=colour, lw=2.2, label=f"starts at {T0} °C, 500 bar")
    boiling = path.has("vapour") & path.has("liquid")
    w_vapour = np.where(boiling, path.S_v * path.Rho_v / path.Rho, np.nan)   # vapour mass fraction (two phases only)
    bx.plot(p / 1e5, 100 * w_vapour, color=colour, lw=2.2)
    i = np.flatnonzero(boiling)[0]
    ax.plot(path.T[i] - 273.15, p[i] / 1e5, "o", ms=7, mfc=colour, mec="white", mew=1.5)

ticks = [50, 100, 200, 300, 500]                 # bar
ax.set_yticks(ticks, [str(t) for t in ticks])
ax.set_title("Path in temperature and pressure")
ax.legend(loc="lower left", fontsize="small")
ax.text(300, 300, "single-phase fluid", ha="center", fontsize="small")
ax.text(300, 75, "vapour + liquid", ha="center", fontsize="small")
ax.text(470, 70, "vapour\n+ halite", ha="center", fontsize="small")

bx.set_xscale("log")
bx.set_xticks(ticks, [str(t) for t in ticks])
bx.invert_xaxis()
bx.annotate("starts as a vapour-like fluid:\na brine condenses", (325, 70), xytext=(240, 88), fontsize="small",
            color=st.CORAL, arrowprops=dict(arrowstyle="-", color=st.CORAL))
bx.text(150, 8, "start as liquid-like fluids: they boil", fontsize="small", color=st.INK)
bx.set_xlabel("Pressure (bar)")
bx.set_ylabel("vapour mass fraction (%)")
bx.set_ylim(0, 100)
bx.set_title("Vapour mass fraction once two phases coexist")
fig.suptitle("Seawater rising from 500 bar at constant specific enthalpy (no heat exchange)", x=0.01, ha="left",
             color=st.NAVY, fontweight="bold")
plt.show()
# --8<-- [end:figure]
fig.savefig(os.path.join(st.outdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img")),
                         "upflow_seawater.png"))
print("upflow_seawater.png written")
