"""Gallery and "Phase boundaries" page: the phase boundaries of H2O-NaCl in T-p-X (phase_diagram_3d.png) and the
three-phase curve V+L+H (vlh_curve.png).

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

# --8<-- [start:diagram3d]
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_hex
from xthermo import H2ONaCl

warnings.filterwarnings("ignore", category=RuntimeWarning)
sw = H2ONaCl.cH2ONaCl("IAPS84")

# ── VLH surface (V→L side and L→H side) ──────────────────────────────────
nT = 80
nT_high = nT // 3
T_high  = H2ONaCl.T_MIN_VLH + (sw.Tmax_VLH() - H2ONaCl.T_MIN_VLH) * 0.95
T_vlh   = np.append(np.linspace(H2ONaCl.T_MIN_VLH, T_high, nT - nT_high),
                    np.linspace(T_high, sw.Tmax_VLH(), nT_high))
P_vlh   = np.array(sw.P_VLH(T_vlh))
Xl_vlh, Xv_vlh = np.array(sw.X_VLH(T_vlh, P_vlh))

n_log, n_lin = 20, 40
nX = n_log + n_lin
TT_vlh = np.tile(T_vlh, (nX, 1))
PP_vlh = np.tile(P_vlh, (nX, 1))
XX_V2L = np.zeros_like(TT_vlh)
XX_L2H = np.zeros_like(TT_vlh)
for j in range(len(T_vlh)):
    XX_V2L[:, j] = np.append(
        10 ** np.linspace(np.log10(max(Xv_vlh[j], 1e-10)), np.log10(0.01), n_log),
        np.linspace(0.01, Xl_vlh[j], n_lin))
    XX_L2H[:, j] = np.linspace(Xl_vlh[j], 1.0, nX)

# ── VH surface ─────────────────────────────────────────────────────────────
nP = 60
T_vh = np.linspace(H2ONaCl.T_MIN_VLH, sw.Tmax_VLH(), nT)
P_vlh2 = np.array(sw.P_VLH(T_vh))
TT_vh = np.tile(T_vh, (nP, 1))
PP_vh = np.zeros_like(TT_vh)
np_low = nP // 3
for i in range(len(T_vh)):
    p_low = sw.pmin() + (P_vlh2[i] - sw.pmin()) * 0.1
    PP_vh[:, i] = np.append(np.linspace(sw.pmin(), p_low, np_low),
                            np.linspace(p_low, P_vlh2[i], nP - np_low))
XX_VH = np.array(sw.X_VH(TT_vh.ravel(), PP_vh.ravel())).reshape(TT_vh.shape)

# ── Halite liquidus ────────────────────────────────────────────────────────
TT_lh = np.tile(T_vh, (nP, 1))
PP_lh = np.zeros_like(TT_lh)
for i in range(len(T_vh)):
    PP_lh[:, i] = np.linspace(P_vlh2[i], 2500e5, nP)
XL_LH = np.array(sw.X_HaliteLiquidus(TT_lh.ravel(), PP_lh.ravel())).reshape(TT_lh.shape)

# ── VL surfaces ────────────────────────────────────────────────────────────
pb_l = sw.PhaseBoundary_VL_DeformLinear(H2ONaCl.Liquid)
pb_v = sw.PhaseBoundary_VL_DeformLinear(H2ONaCl.Vapor)

# ── Critical curve ─────────────────────────────────────────────────────────
T_crit = np.linspace(sw.get_pWater().T_critical(), sw.Tmax(), 80)
p_crit, X_crit = np.array(sw.P_X_Critical(T_crit))

# ── Plot ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(111, projection='3d', facecolor='None')
legend = []   # 3-D surfaces give no usable legend handles: one coloured swatch per surface

def surf(TT, PP, XX, color, label):
    ax.plot_surface(XX * 100, TT - 273.15, PP / 1e5,
                    color=color, alpha=0.45, linewidth=0)
    ax.plot_wireframe(XX * 100, TT - 273.15, PP / 1e5,
                      color=color, lw=0.4)
    legend.append(Patch(fc=color, ec=color, alpha=0.6, label=label))

surf(TT_vlh, PP_vlh, XX_V2L,  'steelblue', 'VLH (V side)')
surf(TT_vlh, PP_vlh, XX_L2H,  'orange',    'VLH (L+H side)')
surf(TT_vh,  PP_vh,  XX_VH,   'purple',    'VH surface')
surf(TT_lh,  PP_lh,  XL_LH,   'limegreen', 'Halite liquidus')

for pb, color, label in [(pb_l, 'black', 'VL liquid'), (pb_v, 'gray', 'VL vapor')]:
    T2 = np.array(pb.T); P2 = np.array(pb.p); X2 = np.array(pb.X)
    ax.plot_surface(X2 * 100, T2 - 273.15, P2 / 1e5, color=color, alpha=0.3, linewidth=0)
    ax.plot_wireframe(X2 * 100, T2 - 273.15, P2 / 1e5, color=color, lw=0.4)
    legend.append(Patch(fc=color, ec=color, alpha=0.4, label=label))

ax.plot(X_crit * 100, T_crit - 273.15, p_crit / 1e5, color='red', lw=2)
legend.append(Line2D([], [], color='red', lw=2, label='Critical curve'))
ax.plot(Xv_vlh * 100, T_vlh - 273.15, P_vlh / 1e5, color='red',   lw=1)
ax.plot(Xl_vlh * 100, T_vlh - 273.15, P_vlh / 1e5, color='green', lw=1)

ax.set_xlabel('Salinity (wt% NaCl)', fontsize=10, labelpad=8)
ax.set_ylabel('Temperature (°C)',    fontsize=10, labelpad=8)
ax.set_zlabel('Pressure (bar)',      fontsize=10, labelpad=8)
ax.set_xlim(0, 100); ax.set_ylim(0, 1000); ax.set_zlim(0, 2500)
ax.view_init(elev=25, azim=-145)
ax.set_title(r'H$_2$O–NaCl phase boundaries in T–p–X space', fontsize=13, pad=12)
ax.legend(handles=legend, loc='upper left', fontsize=9)
ax.set_box_aspect(None, zoom=0.88)   # leaves room for the axis labels
plt.show()
# --8<-- [end:diagram3d]
fig.savefig(os.path.join(OUT, "phase_diagram_3d.png"), dpi=150, bbox_inches="tight", pad_inches=0.3)
plt.close("all")

# --8<-- [start:vlh]
import matplotlib.pyplot as plt
import numpy as np
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

T = np.linspace(H2ONaCl.T_MIN_VLH, H2ONaCl.T_MAX_VLH, 200)   # K
P_vlh = np.array(sw.P_VLH(T))                                  # Pa
XL_vlh, XV_vlh = np.array(sw.X_VLH(T, P_vlh))                 # kg/kg

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].plot(T - 273.15, P_vlh / 1e5)
axes[0].set_xlabel("Temperature (°C)")
axes[0].set_ylabel("Pressure (bar)")
axes[0].set_title("VLH coexistence pressure")

axes[1].plot(T - 273.15, XL_vlh * 100, label="Liquid (halite liquidus)")
axes[1].plot(T - 273.15, XV_vlh * 100, label="Vapour")
axes[1].set_yscale("log")
axes[1].set_xlabel("Temperature (°C)")
axes[1].set_ylabel("Salinity (wt% NaCl)")
axes[1].legend()
plt.tight_layout()
plt.show()
# --8<-- [end:vlh]
fig.savefig(os.path.join(OUT, "vlh_curve.png"), dpi=150, bbox_inches="tight")
plt.close("all")
print("phase_diagram_3d.png vlh_curve.png written to", os.path.abspath(OUT))
