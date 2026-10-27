"""Python guide, "Phase diagrams" page: the figures of docs/python/phase_diagrams.md.

The page shows each section between the snippet markers verbatim; everything outside the markers only saves the
figures for the documentation (docs/assets/img, or the directory given as the first argument or in
$XTHERMO_FIGURE_DIR). Writes pslice_300bar.png, tslice_450C.png, slice_25bar.png.
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

# --8<-- [start:pslice]
import matplotlib.pyplot as plt
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
r = xthermo.plot_Pslice(300e5, sw=sw)        # isobaric section at 300 bar: T-X (a) and H-X (b)
plt.show()
# --8<-- [end:pslice]
r["fig"].savefig(os.path.join(OUT, "pslice_300bar.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# --8<-- [start:tslice]
import matplotlib.pyplot as plt
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
r = xthermo.plot_Tslice(450 + 273.15, sw=sw)   # isothermal section at 450 °C: P-X
plt.show()
# --8<-- [end:tslice]
r["fig"].savefig(os.path.join(OUT, "tslice_450C.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# --8<-- [start:slice]
import matplotlib.pyplot as plt
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
section = sw.Slice_constP(25e5)                  # the phase-region geometry at 25 bar

fig, (ax_t, ax_h) = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
for ax in (ax_t, ax_h):
    ax.set_xscale("log")                         # salinity on a log axis
    ax.set_xlim(1e-4, 100)                       # wt% NaCl
handles = xthermo.plot_slice(ax_t, section, x="X", y="T")
xthermo.plot_slice(ax_h, section, x="X", y="H")
ax_t.set_ylim(0, 1000)
fig.legend(handles.values(), handles.keys(), loc="outside lower center", ncol=6, fontsize="small")
plt.show()
# --8<-- [end:slice]
fig.savefig(os.path.join(OUT, "slice_25bar.png"), dpi=150, bbox_inches="tight")
plt.close("all")
print("pslice_300bar.png tslice_450C.png slice_25bar.png written to", os.path.abspath(OUT))
