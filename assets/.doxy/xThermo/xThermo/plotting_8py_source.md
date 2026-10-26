

# File plotting.py

[**File List**](files.md) **>** [**API**](dir_4b630c9df1aec2cd2fe9acaf75644364.md) **>** [**python**](dir_a5d5f7c295cec99c8e29471acd4df0e8.md) **>** [**xthermo**](dir_c7c9957b02953ad15ede12b0595b5fff.md) **>** [**plotting.py**](plotting_8py.md)

[Go to the documentation of this file](plotting_8py.md)


```Python
"""Plotting helpers (require matplotlib): phase-diagram sections from the library's phase-boundary geometry and property
maps from xthermo.tpx / xthermo.hpx results. Imported on first use (xthermo.plot_slice etc.), so "import xthermo"
does not import matplotlib.

    import matplotlib.pyplot as plt, numpy as np, xthermo
    from xthermo import H2ONaCl
    sw = H2ONaCl.cH2ONaCl("IAPS84")
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    xthermo.plot_slice(ax, sw.Slice_constP(300e5), x="X", y="T", legend=True)   # isobaric T-X section

    T, X = np.linspace(280, 990, 300) + 273.15, np.logspace(-4, 0, 300)
    s = xthermo.tpx(sw, T[:, None], 300e5, X[None, :], errors="nan")
    xthermo.plot_property(ax, s, "Rho", x="X", y="T", contours=10)              # property map, contours per phase

Display units (axes): T in deg C, P in bar, X in wt% NaCl, H in MJ/kg, as in the former plot_Pslice. Property values
(colour scales) are the library's SI values.

Sections: plot_slice draws a Slice_constP (isobaric: T-X, H-X) or Slice_constT (isothermal: P-X) result: the library's
multiphase region polygons with their colours, boundary lines and invariant points. Unshaded area is single-phase fluid.
plot_Pslice / plot_Tslice build complete figures on log-linear salinity axes (makeLogLinearAxes, createWaterAxes).

Grids: plot_property and plot_phases take a 2-D xthermo.States (inputs broadcast to a grid, e.g. T[:, None], X[None, :])
and draw one cell per point, without interpolation. Contours of plot_property are drawn separately in each phase
assemblage (liquid, vapour, L+H, V+L, ...), so no contour line is interpolated across a phase boundary; failed points
(not ok) are masked and shown in grey.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.ticker import LogLocator, MultipleLocator

from . import H2ONaCl
from .states import _PRESENT

__all__ = ["plot_slice", "plot_Pslice", "plot_Tslice", "plot_property", "plot_phases", "makeLogLinearAxes",
           "createWaterAxes", "UNITS"]

# display conversion of the coordinates: SI value -> axis value, axis label
UNITS = {"T": (lambda v: v - 273.15, "Temperature (°C)"), "P": (lambda v: v / 1e5, "Pressure (bar)"),
         "X": (lambda v: v * 100, "Bulk salinity (wt% NaCl)"), "H": (lambda v: v / 1e6, "Specific enthalpy (MJ/kg)")}
_STATE_FIELD = {"T": "T", "P": "p", "X": "X", "H": "H"}
# SI units of the property fields (DataStructures.h), for colour-bar labels
_FIELD_UNITS = {"T": "K", "p": "Pa", "X": "kg/kg", "H": "J/kg", "Rho": "kg/m³", "Cp": "J/(kg K)", "Mu": "Pa s",
                "dRhodP": "kg/(m³ Pa)", "dRhodT": "kg/(m³ K)", "IsothermalCompressibility": "1/Pa",
                "IsobaricExpansivity": "1/K", "S": "–"}
# phase codes (PhaseRegion, DataStructures.h) -> colour for plot_phases; multiphase colours as the slice regions
_PHASE_COLOURS = {0: "#b8e186", 1: "#fdbf6f", 2: "#8c8c8c", 3: "#e0e0f8", 4: "#fee6c8", 5: "#e3f2c4", 6: "#e31a1c",
                  7: "#f0e442", 8: (0.68, 0.85, 0.90), 9: (0.66, 0.66, 0.66), 10: (0.87, 0.63, 0.87),
                  11: (1.0, 0.87, 0.68), 12: "#ffffff", 13: "#ffffff"}
_SINGLE_PHASE_LABEL = "Single-phase fluid"
# phase assemblage (phases present, states._PRESENT) -> name; single-phase fluid states above the critical point are
# reported in the liquid fields (library convention), so they belong to "L"
_ASSEMBLAGE_NAME = {"l": "L", "v": "V", "h": "H", "lv": "V+L", "lh": "L+H", "vh": "V+H", "lvh": "V+L+H"}


def _coords(obj, k):
    conv, _ = UNITS[k]
    return conv(np.asarray(getattr(obj, k), dtype=float))


def _axes_list(ax):
    return list(ax) if isinstance(ax, (list, tuple)) else [ax]


def slice_kind(sl):
    """("isobaric", p [Pa]) for a Slice_constP result, ("isothermal", T [K]) for a Slice_constT result."""
    for name in sl.regions.keys():
        r = sl.regions[name][0]
        P, T = np.asarray(r.P, float), np.asarray(r.T, float)
        if np.all(P == P[0]):
            return "isobaric", float(P[0])
        if np.all(T == T[0]):
            return "isothermal", float(T[0])
    raise ValueError("cannot tell whether the slice is isobaric or isothermal")


def _in_fixed_limits(a, px, py):
    """False if a point lies outside an axis range that was set explicitly (autoscaling off) or on a log axis at <= 0."""
    for v, scale, lim, auto in ((px, a.get_xscale(), a.get_xlim(), a.get_autoscalex_on()),
                                (py, a.get_yscale(), a.get_ylim(), a.get_autoscaley_on())):
        if scale == "log" and v <= 0:
            return False
        if not auto and not min(lim) <= v <= max(lim):
            return False
    return True


def plot_slice(ax, sl, x="X", y=None, *, regions=True, lines=True, points=True, legend=False):
    """Draw a phase-region section (the result of cH2ONaCl.Slice_constP(p) or Slice_constT(T)) on ax.

    ax: an Axes or a list of Axes (e.g. the two halves of makeLogLinearAxes); everything is drawn on each.
    x, y: coordinates "T", "P", "X" or "H" (display units, see UNITS). y defaults to "T" for an isobaric and "P" for an
    isothermal section. Isobaric sections: T-X or H-X; isothermal: P-X (or H-X).
    The region polygons, colours, boundary lines and invariant points are the library's. Unshaded area is single-phase
    fluid. The three-phase V+L+H region of an isobaric section is a line in T-X and a triangle in H-X.
    legend: add a legend to the (last) Axes. Returns {label: artist} for building a legend.
    """
    if y is None:
        y = "T" if slice_kind(sl)[0] == "isobaric" else "P"
    for k in (x, y):
        if k not in UNITS:
            raise ValueError(f"coordinates must be in {sorted(UNITS)}, got {k!r}")
    axes = _axes_list(ax)
    handles = {}
    for a in axes:
        if regions:
            # L+H first: the L+H polygon of Slice_constT follows the halite liquidus below the V+L+H pressure, where
            # it overlaps V+H (V+H is the stable assemblage there), so V+H must be drawn over it
            for name in sorted(sl.regions.keys(), key=lambda n: n != "L+H"):
                for r in sl.regions[name]:
                    vlh = name.startswith("V+L+H")
                    h = a.fill(_coords(r, x), _coords(r, y), fc=tuple(r.fc), ec=tuple(r.ec) if vlh else "none",
                               lw=0.8 if vlh else 0, zorder=1)[0]
                    handles.setdefault(name, h)
        if lines:
            for name in sl.lines.keys():
                for ln in sl.lines[name]:
                    a.plot(_coords(ln, x), _coords(ln, y), color=tuple(ln.color), ls=ln.ls, lw=ln.lw, zorder=2)
        if points:
            for name in sl.points.keys():
                for pt in sl.points[name]:
                    px, py = float(_coords(pt, x)), float(_coords(pt, y))
                    if not _in_fixed_limits(a, px, py):   # e.g. X = 0 on the linear half of a log-linear axis
                        continue
                    h = a.scatter([px], [py], marker=pt.marker, fc=[tuple(pt.mfc)],
                                  ec=[tuple(pt.mec)], lw=0.5, zorder=4)
                    handles.setdefault(name, h)
        a.set_xlabel(UNITS[x][1])
        a.set_ylabel(UNITS[y][1])
    if regions:
        handles[_SINGLE_PHASE_LABEL] = patches.Patch(fc="white", ec="0.5", lw=0.5)
    if legend:
        axes[-1].legend(handles.values(), handles.keys(), fontsize="small")
    return handles


def makeLogLinearAxes(ax, width_logaxis=0.4, xlim_log=(2E-3, 1), xmax=100, xloc_major=20, xloc_minor=4, yloc='left',
                      color_log=(127/255, 0, 1), xlabel='Bulk salinity (wt % NaCl)', bkcolor='None', ylim=None,
                      yloc_major=None, yloc_minor=None, ylabel=None):
    """Replace ax by a salinity axis that is logarithmic up to xlim_log[1] (width fraction width_logaxis) and linear up
    to xmax (wt% NaCl). Returns (ax_log, ax_linear), which share the y axis; y tick labels on the side yloc."""
    # log segment
    ax_log = ax.inset_axes([0, 0, width_logaxis, 1], facecolor=bkcolor)
    ax_log.set_xscale('log')
    ax_log.set_xlim(xlim_log)
    ax_log.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    ax_log.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    ax_log.spines['right'].set_visible(False)
    ax_log.spines['bottom'].set_color(color_log)
    ax_log.spines['top'].set_color(color_log)
    ax_log.tick_params(axis='x', which='both', colors=color_log)
    ax_log.annotate("", xy=(0.03, -0.14), xycoords='axes fraction', xytext=(1, -0.14), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color=color_log), annotation_clip=False,
                    zorder=0)
    ax_log.set_xlabel('Log-scale', color=color_log, bbox={'fc': 'w', 'ec': 'None'})

    # linear segment, sharing the y axis of the log segment
    ax_linear = ax.inset_axes([width_logaxis, 0, 1 - width_logaxis, 1], facecolor=bkcolor, sharey=ax_log)
    ax_linear.set_xlim(xlim_log[1], xmax)
    ax_linear.xaxis.set_major_locator(MultipleLocator(xloc_major))
    ax_linear.xaxis.set_minor_locator(MultipleLocator(xloc_minor))
    ax_linear.spines['left'].set_linewidth(0.5)
    ax_linear.spines['left'].set_color(color_log)
    ax_linear.annotate("", xy=(-0.03, -0.14), xycoords='axes fraction', xytext=(1, -0.14), textcoords='axes fraction',
                       arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"), annotation_clip=False, zorder=0)
    ax_linear.set_xlabel('Linear-scale', labelpad=10, bbox={'fc': 'w', 'ec': 'None'})
    # y ticks: labels on one side only (tick_params, not set_ticklabels: the tickers are shared)
    ax_linear.yaxis.set_ticks_position('right')
    if yloc == 'left':
        ax_log.tick_params(axis='y', which='both', left=True, labelleft=True)
        ax_linear.tick_params(axis='y', which='both', labelright=False)
    else:
        ax_log.tick_params(axis='y', which='both', labelleft=False)
        ax_linear.tick_params(axis='y', which='both', right=True, labelright=True)
        ax_linear.yaxis.set_label_position('right')
    if ylim is not None:
        ax_log.set_ylim(ylim)
    if yloc_major is not None:
        ax_log.yaxis.set_major_locator(MultipleLocator(yloc_major))
    if yloc_minor is not None:
        ax_log.yaxis.set_minor_locator(MultipleLocator(yloc_minor))
    if ylabel is not None:
        ax_log.set_ylabel(ylabel) if yloc == 'left' else ax_linear.set_ylabel(ylabel)
    # hide the base axes, keep its x label as the common salinity label
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if xlabel is not None:
        ax.set_xlabel(xlabel, labelpad=40)
    return ax_log, ax_linear


def createWaterAxes(ax, fc='whitesmoke'):
    """A narrow strip at the left of ax (the log salinity axis) for pure water, X = 0. Returns the strip's Axes,
    sharing the y axis of ax."""
    ax_water = ax.inset_axes([0.004, 0, 0.1, 1], facecolor=fc, sharey=ax)
    ax_water.spines['left'].set_visible(False)
    ax_water.spines['right'].set_visible(False)
    ax_water.set_xlim(-0.5, 1)
    ax_water.set_xticks([0])
    ax_water.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
    ax_water.axvspan(0.7, 1, hatch='////', fc=fc, ec='k', lw=0)
    return ax_water


def _water_strip(ax_water, sl, y, P0, sw):
    """Pure water (X = 0) on the strip of an isobaric section: boiling point (T) or liquid / V+L / vapour (H)."""
    two_phase = [ln for name in sl.lines.keys() if name == "Two phase water" for ln in sl.lines[name]]
    if P0 > sw.get_pWater().p_critical():
        ax_water.text(0.1, 0.5, 'Supercritical water', transform=ax_water.transAxes, ha='center', va='center',
                      rotation=90, fontsize='small')
        return
    for ln in two_phase:
        v = _coords(ln, y)
        if y == "T":
            ax_water.axhline(v[0], xmax=0.7, color=tuple(ln.color), lw=2)
            ax_water.text(0.1, v[0], 'boiling', ha='center', va='bottom', fontsize='x-small')
        else:
            lo, hi = sorted(v[:2])
            y0, y1 = ax_water.get_ylim()
            for a, b, c, t in ((y0, lo, 'lightgreen', 'Liquid'), (lo, hi, 'gold', 'V+L'), (hi, y1, 'lightcyan', 'Vapour')):
                ax_water.axhspan(a, b, 0.1, 0.6, fc=c, ec='none')
                ax_water.text(0.1, (a + b) / 2, t, ha='center', va='center', rotation=90, fontsize='x-small')


def plot_Pslice(P0, x0_log_linear=0.01, sw=None):
    """Isobaric section at pressure P0 [Pa]: T-X (a) and H-X (b) on log-linear salinity axes (log below
    x0_log_linear [kg/kg], linear above), with a strip for pure water. sw: a cH2ONaCl model (default IAPS84).
    Returns {'fig', 'data' (the Slice_constP result), 'ax_T': [log, linear], 'ax_H': [log, linear], 'ax_water': [T, H],
    'handles'}."""
    if sw is None:
        sw = H2ONaCl.cH2ONaCl("IAPS84")
    sl = sw.Slice_constP(P0)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), gridspec_kw={'wspace': 0.08})
    xlim_log = (2E-4, x0_log_linear * 100)   # wt% NaCl
    Hmax = max(_coords(r, "H").max() for name in sl.regions.keys() for r in sl.regions[name])
    Hmax = np.ceil(Hmax * 5) / 5
    ax_T = makeLogLinearAxes(axes[0], xlim_log=xlim_log, ylim=(0, 1000), yloc_major=200, yloc_minor=40)
    ax_H = makeLogLinearAxes(axes[1], xlim_log=xlim_log, ylim=(0, Hmax), yloc='right', yloc_major=1, yloc_minor=0.2)
    handles = plot_slice(list(ax_T), sl, x="X", y="T")
    plot_slice(list(ax_H), sl, x="X", y="H")
    for a in ax_T + ax_H:   # plot_slice labels each axes; the hybrid axes carry their own labels
        a.set_xlabel('Log-scale' if a.get_xscale() == 'log' else 'Linear-scale')
        a.set_ylabel('')
    ax_T[0].set_ylabel('Temperature (°C)')
    ax_H[1].set_ylabel('Specific enthalpy (MJ/kg)')
    ax_water = [createWaterAxes(ax_T[0]), createWaterAxes(ax_H[0])]
    _water_strip(ax_water[0], sl, "T", P0, sw)
    _water_strip(ax_water[1], sl, "H", P0, sw)
    for a, tag in ((ax_T[0], '(a)'), (ax_H[0], '(b)')):
        a.text(0.97, 0.98, tag, transform=a.transAxes, fontsize=12, fontweight='bold', va='top', ha='right',
               bbox={'fc': 'w', 'ec': 'None', 'alpha': 0.9}, zorder=9)
        a.text(0.97, 0.02, 'P = %.0f bar' % (P0 / 1E5), transform=a.transAxes, ha='right', va='bottom', zorder=9)
    fig.legend(handles.values(), handles.keys(), loc='lower center', bbox_to_anchor=(0.5, 0.93), ncol=6,
               fontsize='small', frameon=False, columnspacing=0.8, handletextpad=0.3, handlelength=1.2)
    return {'fig': fig, 'data': sl, 'ax_T': list(ax_T), 'ax_H': list(ax_H), 'ax_water': ax_water, 'handles': handles}


def plot_Tslice(T0, x0_log_linear=0.01, sw=None, plog=True):
    """Isothermal P-X section at temperature T0 [K] on a log-linear salinity axis (log below x0_log_linear [kg/kg]),
    pressure on a log axis (plog=True) or linear. sw: a cH2ONaCl model (default IAPS84).
    Returns {'fig', 'data' (the Slice_constT result), 'ax_P': [log, linear], 'handles'}."""
    if sw is None:
        sw = H2ONaCl.cH2ONaCl("IAPS84")
    sl = sw.Slice_constT(T0)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax_P = makeLogLinearAxes(ax, xlim_log=(2E-4, x0_log_linear * 100), ylabel='Pressure (bar)')
    if plog:
        ax_P[0].set_yscale('log')
    ax_P[0].set_ylim(sw.pmin() / 1e5, sw.pmax() / 1e5)
    handles = plot_slice(list(ax_P), sl, x="X", y="P")
    for a in ax_P:
        a.set_xlabel('Log-scale' if a.get_xscale() == 'log' else 'Linear-scale')
        a.set_ylabel('')
    ax_P[0].set_ylabel('Pressure (bar)')
    ax_P[1].text(0.97, 0.98, 'T = %.0f °C' % (T0 - 273.15), transform=ax_P[1].transAxes, ha='right', va='top')
    fig.legend(handles.values(), handles.keys(), loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=5,
               fontsize='small', frameon=False, columnspacing=0.8, handletextpad=0.3, handlelength=1.2)
    return {'fig': fig, 'data': sl, 'ax_P': list(ax_P), 'handles': handles}


def _name(spec):
    name = spec[0] if isinstance(spec, tuple) else spec
    if name not in UNITS:
        raise ValueError(f"coordinates must be in {sorted(UNITS)}, got {name!r}")
    return name


def _grid(states, x, y):
    """1-D display coordinates of a 2-D States grid along its two axes, and whether x varies along axis 0.

    x, y: a coordinate name, taken from the result fields, or (name, values) with the 1-D SI input values. Result fields
    must form an exact grid: T, p and X of tpx and p and X of hpx do (they are the inputs); H of hpx is the evaluated
    enthalpy, which can differ from the input, so H grids need (name, values)."""
    if len(states.shape) != 2:
        raise ValueError(f"a 2-D grid of states is required (inputs broadcast to two dimensions), got shape {states.shape}")
    n0, n1 = states.shape

    def axis_values(spec):   # -> (values, axis) with axis 0, 1 or None (either; square grid with given values)
        name = _name(spec)
        if isinstance(spec, tuple):
            v = UNITS[name][0](np.asarray(spec[1], float).ravel())
            axes = [a for a, n in ((1, n1), (0, n0)) if v.size == n]
            if not axes:
                raise ValueError(f"{len(v)} {name} values do not match the grid shape {states.shape}")
            return v, (axes[0] if len(axes) == 1 else None)
        c = UNITS[name][0](np.asarray(states[_STATE_FIELD[name]], float))
        if np.all(c == c[:1, :]):
            return c[0, :], 1
        if np.all(c == c[:, :1]):
            return c[:, 0], 0
        hint = f"; pass {name} as ({name!r}, values) with the input values" if name == "H" else ""
        raise ValueError(f"the result field {name} does not vary along one axis of the grid{hint}")

    xv, ax_x = axis_values(x)
    yv, ax_y = axis_values(y)
    if ax_x is None:
        ax_x = 1 - ax_y if ax_y is not None else 1
    if ax_y is None:
        ax_y = 1 - ax_x
    if ax_x == ax_y:
        raise ValueError(f"{_name(x)} and {_name(y)} vary along the same axis of the grid")
    return xv, yv, ax_x == 0


def _edges(v):
    """Cell edges around the points v: midpoints in log space for logarithmically spaced positive values."""
    v = np.asarray(v, float)
    if v.size < 2:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    if np.all(v > 0):
        d = np.diff(np.log(v))
        if np.allclose(d, d[0], rtol=1e-6, atol=0):
            lv = np.log(v)
            m = (lv[1:] + lv[:-1]) / 2
            return np.exp(np.r_[2 * lv[0] - m[0], m, 2 * lv[-1] - m[-1]])
    m = (v[1:] + v[:-1]) / 2
    return np.r_[2 * v[0] - m[0], m, 2 * v[-1] - m[-1]]


def _assemblage(states, grid_t):
    phase = np.asarray(states["phase"]).astype(int)
    keys = np.array([_PRESENT.get(int(c), "") for c in np.unique(phase)])
    lut = {int(c): k for c, k in zip(np.unique(phase), keys)}
    a = np.vectorize(lambda c: lut[int(c)], otypes=[object])(phase)
    return a.T if grid_t else a


def plot_property(ax, states, field, x="X", y="T", *, cmap="viridis", norm=None, vmin=None, vmax=None, contours=None,
                  contour_kw=None, colorbar=True, failed="0.75", label=None):
    """Map of a property field of a 2-D xthermo.States grid, one cell per point (no interpolation).

    x, y: the grid's coordinates "T", "P", "X" or "H" (display units), taken from the result, or (name, values) with the
    1-D SI input values, needed for H of hpx grids (e.g. y=("H", H)). field: a States field ("Rho", "Cp", "Mu_l", ...),
    in SI units. Points that were not evaluated (not states.ok) are masked and drawn in the colour failed (None: not
    drawn); NaN values of evaluated points (e.g. an absent phase with absent="nan") are left blank.
    contours: number of levels or a sequence of levels; contour lines are drawn separately in each phase assemblage
    (the set of phases present, from the phase code), so they end at phase boundaries instead of being interpolated
    across them. norm: a matplotlib norm, e.g. matplotlib.colors.LogNorm() for viscosity.
    Returns {'mesh', 'contours' ({assemblage: ContourSet}, e.g. "L", "V", "V+L"), 'colorbar', 'failed' (the mesh of
    failed points or None)}.
    """
    xv, yv, grid_t = _grid(states, x, y)
    Z = np.asarray(states[field], float)
    ok = np.asarray(states.ok)
    if grid_t:
        Z, ok = Z.T, ok.T
    Zm = np.ma.masked_where(~ok | ~np.isfinite(Z), Z)
    xe, ye = _edges(xv), _edges(yv)
    mesh = ax.pcolormesh(xe, ye, Zm, cmap=cmap, norm=norm, vmin=None if norm else vmin, vmax=None if norm else vmax,
                         shading="flat", zorder=0.5)
    fmesh = None
    if failed is not None and not ok.all():
        fmesh = ax.pcolormesh(xe, ye, np.ma.masked_where(ok, np.ones_like(Z)),
                              cmap=mpl.colors.ListedColormap([failed]), shading="flat", zorder=0.6)
    sets = {}
    if contours is not None and Zm.count():
        levels = contours
        if np.ndim(contours) == 0:
            lo, hi = float(Zm.min()), float(Zm.max())
            if isinstance(norm, mpl.colors.LogNorm) and lo > 0:
                levels = np.geomspace(lo, hi, int(contours) + 2)[1:-1]
            else:
                levels = mpl.ticker.MaxNLocator(int(contours)).tick_values(lo, hi)
        kw = dict(colors="k", linewidths=0.6)
        kw.update(contour_kw or {})
        groups = _assemblage(states, grid_t)
        for g in np.unique(groups[ok]):
            Zg = np.ma.masked_where((groups != g) | ~ok | ~np.isfinite(Z), Z)
            if Zg.count() >= 4:
                sets[_ASSEMBLAGE_NAME.get(g, g)] = ax.contour(xv, yv, Zg, levels=levels, **kw)
    cb = None
    if colorbar:
        unit = _FIELD_UNITS.get(field.split("_")[0] if field[:2] != "S_" else "S", "")
        cb = ax.figure.colorbar(mesh, ax=ax, label=label or (f"{field} ({unit})" if unit else field))
    ax.set_xlabel(UNITS[_name(x)][1])
    ax.set_ylabel(UNITS[_name(y)][1])
    return {"mesh": mesh, "contours": sets, "colorbar": cb, "failed": fmesh}


def plot_phases(ax, states, x="X", y="T", *, legend=True, failed="0.75"):
    """Phase map of a 2-D xthermo.States grid (x, y as for plot_property): one cell per point coloured by its phase code (no interpolation between
    phase codes); failed points in the colour failed. For sections of H2O-NaCl at constant p or T, plot_slice (the
    library's region geometry) is the more accurate diagram; plot_phases also covers other grids, e.g. T-P at fixed X.
    Returns {'mesh', 'handles' ({phase name: artist})}."""
    xv, yv, grid_t = _grid(states, x, y)
    phase = np.asarray(states["phase"]).astype(int)
    ok = np.asarray(states.ok)
    if grid_t:
        phase, ok = phase.T, ok.T
    codes = sorted(set(np.unique(phase[ok]).tolist()))
    colours = [_PHASE_COLOURS.get(c, "white") for c in codes]
    index = np.searchsorted(codes, np.where(ok, phase, codes[0] if codes else 0))
    data = np.ma.masked_where(~ok, index)
    mesh = None
    if codes:
        mesh = ax.pcolormesh(_edges(xv), _edges(yv), data, cmap=mpl.colors.ListedColormap(colours),
                             norm=mpl.colors.BoundaryNorm(np.arange(len(codes) + 1) - 0.5, len(codes)), shading="flat")
    model = states._model
    handles = {model.phase_name(c): patches.Patch(fc=col, ec="0.5", lw=0.5) for c, col in zip(codes, colours)}
    if failed is not None and not ok.all():
        ax.pcolormesh(_edges(xv), _edges(yv), np.ma.masked_where(ok, np.ones(ok.shape)),
                      cmap=mpl.colors.ListedColormap([failed]), shading="flat")
        handles["not evaluated"] = patches.Patch(fc=failed, ec="0.5", lw=0.5)
    ax.set_xlabel(UNITS[_name(x)][1])
    ax.set_ylabel(UNITS[_name(y)][1])
    if legend:
        ax.legend(handles.values(), handles.keys(), fontsize="small")
    return {"mesh": mesh, "handles": handles}
```


