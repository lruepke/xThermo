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

if __name__ == '__main__':
    phase_diagram_isobaric()
    boiling_curve()
    vlh_curve()
    print('done')
