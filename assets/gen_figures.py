"""Generate figures for MkDocs documentation site.
Run from repo root: python docs/assets/gen_figures.py
Requires the xThermal Python package to be installed.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import os

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'cm'
OUT = os.path.join(os.path.dirname(__file__), 'img')

# ---------------------------------------------------------------------------
# 1. H2O-NaCl p-T-X phase diagram (isobaric slice, linear scale)
# ---------------------------------------------------------------------------
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

def phase_diagram_isobaric(P0=300e5, out=OUT):
    T = np.linspace(sw.Tmin(), sw.Tmax(), 300)
    X = np.linspace(1e-4, 1.0, 300)
    TT, XX = np.meshgrid(T, X)
    PP = np.full_like(TT, P0)

    state = sw.UpdateState_TPX(TT.ravel(), PP.ravel(), XX.ravel())
    Phase = np.array(state.phase).reshape(TT.shape)

    # colour map
    cmap = plt.get_cmap("Dark2")
    colors = list(copy.deepcopy(cmap.colors))
    colors[:8] = ['#AED6F1','#A9DFBF','#F9E79F','#F5CBA7','#D7BDE2','#FDFEFE','#85C1E9','#A9CCE3']
    cmap = mpl.colors.ListedColormap(colors)

    phase_unique = np.sort(np.unique(Phase))
    phase_names  = [sw.phase_name(int(p)) for p in phase_unique]
    # remap to 0..N for clean colorbar
    Phase_plot = np.zeros_like(Phase)
    for i, p in enumerate(phase_unique):
        Phase_plot[Phase == p] = i

    fig, ax = plt.subplots(figsize=(8, 5))
    cs = ax.contourf(TT - 273.15, XX * 100, Phase_plot,
                     levels=np.arange(-0.5, len(phase_unique)),
                     cmap=cmap, vmin=-0.5, vmax=len(phase_unique) - 0.5)
    cb = fig.colorbar(cs, ax=ax, ticks=np.arange(len(phase_unique)), pad=0.02)
    cb.ax.set_yticklabels(phase_names, fontsize=10)

    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Salinity (wt% NaCl)', fontsize=12)
    ax.set_title(r'H$_2$O–NaCl phase diagram  |  p = ' + f'{P0/1e5:.0f} bar', fontsize=13)
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    ax.tick_params(which='both', direction='in')

    fig.tight_layout()
    path = os.path.join(out, 'phase_diagram_PTX.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f'saved {path}')
    plt.close(fig)

# ---------------------------------------------------------------------------
# 2. H2O boiling curve
# ---------------------------------------------------------------------------
from xThermal import H2O

def boiling_curve(out=OUT):
    water = H2O.cIAPS84()
    T = np.linspace(water.Tmin(), water.T_critical(), 120)
    p = np.array([water.Boiling_p(Ti) for Ti in T])

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(T - 273.15, p / 1e5, color='steelblue', lw=2, label='Boiling curve')
    ax.plot(water.T_critical() - 273.15, water.p_critical() / 1e5,
            'o', color='crimson', ms=8, label=f'Critical point\n({water.T_critical()-273.15:.1f} °C, {water.p_critical()/1e5:.0f} bar)')
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Saturation pressure (bar)', fontsize=12)
    ax.set_title(r'Pure H$_2$O boiling curve (IAPS-84)', fontsize=13)
    ax.set_yscale('log')
    ax.grid(which='major', color='lightgray', lw=0.5)
    ax.grid(which='minor', color='#f0f0f0', lw=0.3)
    ax.legend(fontsize=10)
    ax.tick_params(which='both', direction='in')
    fig.tight_layout()
    path = os.path.join(out, 'boiling_curve_H2O.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f'saved {path}')
    plt.close(fig)

# ---------------------------------------------------------------------------
# 3. VLH pressure curve
# ---------------------------------------------------------------------------
def vlh_curve(out=OUT):
    T = np.linspace(H2ONaCl.T_MIN_VLH, H2ONaCl.T_MAX_VLH, 200)
    P = np.array(sw.P_VLH(T))
    XL, XV = np.array(sw.X_VLH(T, P))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={'wspace': 0.35})

    ax = axes[0]
    ax.plot(T - 273.15, P / 1e5, color='navy', lw=2)
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Pressure (bar)', fontsize=12)
    ax.set_title('VLH coexistence pressure', fontsize=12)
    ax.grid(color='lightgray', lw=0.5)
    ax.tick_params(which='both', direction='in')

    ax = axes[1]
    ax.plot(T - 273.15, XL * 100, color='steelblue', lw=2, label='Liquid (halite liquidus)')
    ax.plot(T - 273.15, XV * 100, color='darkorange', lw=2, label='Vapour')
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Salinity (wt% NaCl)', fontsize=12)
    ax.set_title('VLH phase compositions', fontsize=12)
    ax.set_yscale('log')
    ax.legend(fontsize=10)
    ax.grid(color='lightgray', lw=0.5)
    ax.tick_params(which='both', direction='in')

    fig.tight_layout()
    path = os.path.join(out, 'vlh_curve.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f'saved {path}')
    plt.close(fig)

# ---------------------------------------------------------------------------
# 4. 3-D T-p-X phase boundary diagram
# ---------------------------------------------------------------------------
def phase_diagram_3d(out=OUT):
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    from matplotlib.patches import Patch
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    # --- compute phase boundary surfaces ---
    nT = 80

    # VLH surface
    nT_high = int(nT / 3)
    nT_low  = nT - nT_high
    T_high  = H2ONaCl.T_MIN_VLH + (sw.Tmax_VLH() - H2ONaCl.T_MIN_VLH) * 0.95
    T_vlh   = np.append(np.linspace(H2ONaCl.T_MIN_VLH, T_high, nT_low),
                        np.linspace(T_high, sw.Tmax_VLH(), nT_high))
    P_vlh   = np.array(sw.P_VLH(T_vlh))
    Xl_vlh, Xv_vlh = np.array(sw.X_VLH(T_vlh, P_vlh))

    n_log, n_lin = 20, 40
    X_dummy = np.zeros(n_log + n_lin)
    TT_vlh, _ = np.meshgrid(T_vlh, X_dummy)
    PP_vlh    = np.zeros_like(TT_vlh)
    XX_V2L    = np.zeros_like(TT_vlh)
    XX_L2H    = np.zeros_like(TT_vlh)
    for i in range(PP_vlh.shape[0]):
        PP_vlh[i, :] = P_vlh
    for j in range(TT_vlh.shape[1]):
        XX_V2L[:, j] = np.append(
            10 ** np.linspace(np.log10(max(Xv_vlh[j], 1e-10)),
                              np.log10(0.01), n_log),
            np.linspace(0.01, Xl_vlh[j], n_lin))
        XX_L2H[:, j] = np.linspace(Xl_vlh[j], 1.0, n_log + n_lin)

    # VH surface
    nP = 60
    T_vh = np.linspace(H2ONaCl.T_MIN_VLH, sw.Tmax_VLH(), nT)
    P_vlh2 = np.array(sw.P_VLH(T_vh))
    TT_vh, PP_vh_dummy = np.meshgrid(T_vh, np.zeros(nP))
    PP_vh = np.zeros_like(TT_vh)
    np_low = int(nP / 3)
    np_high = nP - np_low
    for i in range(len(T_vh)):
        p_low = sw.pmin() + (P_vlh2[i] - sw.pmin()) * 0.1
        PP_vh[:, i] = np.append(np.linspace(sw.pmin(), p_low, np_low),
                                np.linspace(p_low, P_vlh2[i], np_high))
    XX_VH = np.array(sw.X_VH(TT_vh.ravel(), PP_vh.ravel())).reshape(TT_vh.shape)

    # Halite liquidus
    TT_lh, _ = np.meshgrid(T_vh, np.zeros(nP))
    PP_lh = np.zeros_like(TT_lh)
    for i in range(len(T_vh)):
        PP_lh[:, i] = np.linspace(P_vlh2[i], 2500e5, nP)
    XL_LH = np.array(sw.X_HaliteLiquidus(TT_lh.ravel(), PP_lh.ravel())).reshape(TT_lh.shape)

    # VL surfaces (liquid and vapor branches)
    pb_l = sw.PhaseBoundary_VL_DeformLinear(H2ONaCl.Liquid)
    pb_v = sw.PhaseBoundary_VL_DeformLinear(H2ONaCl.Vapor)

    # Critical curve
    T_crit = np.linspace(sw.get_pWater().T_critical(), sw.Tmax(), 80)
    p_crit, X_crit = np.array(sw.P_X_Critical(T_crit))

    # --- plot ---
    fig = plt.figure(figsize=(12, 9))
    ax  = fig.add_subplot(111, projection='3d', facecolor='None')

    def surf(TT, PP, XX, color, label):
        ax.plot_surface(XX * 100, TT - 273.15, PP / 1e5,
                        color=color, alpha=0.45, linewidth=0)
        ax.plot_wireframe(XX * 100, TT - 273.15, PP / 1e5,
                          color=color, lw=0.4, label=label)

    surf(TT_vlh, PP_vlh, XX_V2L,  'steelblue', 'VLH (V side)')
    surf(TT_vlh, PP_vlh, XX_L2H,  'orange',    'VLH (L+H side)')
    surf(TT_vh,  PP_vh,  XX_VH,   'purple',    'VH surface')
    surf(TT_lh,  PP_lh,  XL_LH,   'limegreen', 'Halite liquidus')

    for pb, color, label in [(pb_l, 'black', 'VL liquid'), (pb_v, 'gray', 'VL vapor')]:
        TT_vl = np.array(pb.T)
        PP_vl = np.array(pb.p)
        XX_vl = np.array(pb.X)
        ax.plot_surface(XX_vl * 100, TT_vl - 273.15, PP_vl / 1e5,
                        color=color, alpha=0.3, linewidth=0)
        ax.plot_wireframe(XX_vl * 100, TT_vl - 273.15, PP_vl / 1e5,
                          color=color, lw=0.4, label=label)

    ax.plot(X_crit * 100, T_crit - 273.15, p_crit / 1e5,
            color='red', lw=2, label='Critical curve', zorder=10)
    ax.plot(Xv_vlh * 100, T_vlh - 273.15, P_vlh / 1e5, color='red',    lw=1)
    ax.plot(Xl_vlh * 100, T_vlh - 273.15, P_vlh / 1e5, color='green',  lw=1)

    ax.set_xlabel('Salinity (wt% NaCl)', fontsize=10, labelpad=8)
    ax.set_ylabel('Temperature (°C)',    fontsize=10, labelpad=8)
    ax.set_zlabel('Pressure (bar)',      fontsize=10, labelpad=8)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1000)
    ax.set_zlim(0, 2500)
    ax.view_init(elev=25, azim=-145)

    # clean legend with patch handles
    from matplotlib.colors import to_hex
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        if l not in seen:
            color = None
            for attr in ('get_edgecolor', 'get_facecolor', 'get_color'):
                if hasattr(h, attr):
                    try:
                        c = getattr(h, attr)()
                        if hasattr(c, '__len__') and len(c) >= 3:
                            color = to_hex(c[0] if hasattr(c[0], '__len__') else c)
                            break
                    except Exception:
                        pass
            seen[l] = Patch(facecolor='white', edgecolor=color or 'black',
                            linewidth=1.5, label=l)
    ax.legend(handles=list(seen.values()), loc='upper left',
              fontsize=8, ncol=2, framealpha=0.8)

    ax.set_title(r'H$_2$O–NaCl phase boundaries in T–p–X space', fontsize=13, pad=12)
    fig.tight_layout()
    path = os.path.join(out, 'phase_diagram_3d_PTX.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f'saved {path}')
    plt.close(fig)


if __name__ == '__main__':
    phase_diagram_isobaric()
    boiling_curve()
    vlh_curve()
    phase_diagram_3d()
    print('done')
