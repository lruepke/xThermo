"""Gallery: density and salinity of coexisting brine and vapour (coexisting_phases.png).

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
T = np.linspace(380, 800, 421) + 273.15          # K; up to 800 °C (see the text)
colours = {250e5: st.CORAL, 400e5: st.TEAL, 600e5: "#1f6f9a", 1000e5: st.NAVY}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
for p, colour in colours.items():
    # compositions of the coexisting liquid and vapour; a bulk salinity between them lies in the two-phase field
    x_l = np.array(sw.XL_VL(T, np.full_like(T, p)))
    x_v = np.array(sw.XV_VL(T, np.full_like(T, p)))
    s = xthermo.tpx(sw, T, p, 0.5 * (x_l + x_v), errors="nan", absent="nan")
    two_phase = s.ok & s.has("liquid") & s.has("vapour") & ~s.has("halite")
    rho_l = np.where(two_phase, s.Rho_l, np.nan)
    rho_v = np.where(two_phase, s.Rho_v, np.nan)
    label = f"{p / 1e5:.0f} bar"
    ax.plot(T - 273.15, rho_l, color=colour, lw=2, label=label)
    ax.plot(T - 273.15, rho_v, color=colour, lw=2, ls="--")
    bx.semilogy(T - 273.15, np.where(two_phase, 100 * s.X_l, np.nan), color=colour, lw=2, label=label)
    bx.semilogy(T - 273.15, np.where(two_phase, 100 * s.X_v, np.nan), color=colour, lw=2, ls="--")

ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("density (kg/m³)")
ax.set_title("Density: brine (solid) and vapour (dashed)")
ax.legend(fontsize="small")
bx.set_xlabel("Temperature (°C)")
bx.set_ylabel("salinity (wt% NaCl)")
bx.set_title("Salinity: brine (solid) and vapour (dashed)")
bx.legend(fontsize="small", loc="lower right")
fig.suptitle("Coexisting brine and vapour in the vapour + liquid field of H₂O–NaCl", x=0.01, ha="left",
             color=st.NAVY, fontweight="bold")
plt.show()
# --8<-- [end:figure]
fig.savefig(os.path.join(st.outdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img")),
                         "coexisting_phases.png"))
print("coexisting_phases.png written")
