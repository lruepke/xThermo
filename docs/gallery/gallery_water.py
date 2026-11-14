"""Gallery and "Pure water and NaCl" page: water near its critical point (water_near_critical.png) and the boiling
curve of water (boiling_curve_H2O.png).

Writes to docs/assets/img, or to the directory given as the first argument or in $XTHERMO_FIGURE_DIR.
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
warnings.filterwarnings("ignore", message=".*non-interactive.*")   # plt.show() under the Agg backend
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style  # noqa: E402

OUT = _style.outdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "img"))

# --8<-- [start:critical]
import matplotlib.pyplot as plt
import numpy as np
import xthermo
from xthermo import H2O

import _style as st                              # fonts and colours of the documentation figures (docs/gallery)

st.apply()
models = {"IAPS-84": H2O.cIAPS84(), "IAPWS-95": H2O.cIAPWS95()}
T = np.linspace(330, 480, 1501) + 273.15         # K, steps of 0.1 K
colours = {230e5: st.CORAL, 250e5: st.TEAL, 300e5: "#1f6f9a", 400e5: st.NAVY}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
for p, colour in colours.items():
    for (name, water), ls in zip(models.items(), ["-", ":"]):
        s = xthermo.tpx(water, T, p, errors="nan")
        label = f"{p / 1e5:.0f} bar" if name == "IAPS-84" else None
        ax.semilogy(T - 273.15, s.Cp / 1e3, color=colour, ls=ls, lw=1.8, label=label)
        bx.plot(T - 273.15, s.Rho, color=colour, ls=ls, lw=1.8, label=label)
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("isobaric heat capacity (kJ/(kg K))")
ax.set_title("Heat capacity above the critical pressure")
ax.legend(title="solid: IAPS-84, dotted: IAPWS-95", fontsize="small", title_fontsize="small")
bx.set_xlabel("Temperature (°C)")
bx.set_ylabel("density (kg/m³)")
bx.set_title("Density")
fig.suptitle("Pure water near its critical point (374 °C; 220.5 bar IAPS-84, 220.6 bar IAPWS-95)", x=0.01, ha="left", color=st.NAVY,
             fontweight="bold")
plt.show()
# --8<-- [end:critical]
fig.savefig(os.path.join(OUT, "water_near_critical.png"))
plt.close("all")
matplotlib.rcdefaults()                          # the next figure uses Matplotlib's defaults, as its shown code

# --8<-- [start:boiling]
import matplotlib.pyplot as plt
import numpy as np
from xthermo import H2O

water = H2O.cIAPS84()

T = np.linspace(water.Tmin(), water.T_critical(), 100)   # K
p_boil = np.array([water.Boiling_p(Ti) for Ti in T])     # Pa

fig, ax = plt.subplots()
ax.plot(T - 273.15, p_boil / 1e5)
ax.plot(water.T_critical() - 273.15, water.p_critical() / 1e5,
        'ro', label=f"Critical point ({water.T_critical()-273.15:.1f} °C)")
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Saturation pressure (bar)")
ax.set_yscale("log")
ax.legend()
plt.show()
# --8<-- [end:boiling]
fig.savefig(os.path.join(OUT, "boiling_curve_H2O.png"), dpi=150, bbox_inches="tight")
plt.close("all")
print("water_near_critical.png boiling_curve_H2O.png written to", os.path.abspath(OUT))
