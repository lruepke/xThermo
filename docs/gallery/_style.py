"""Shared look of the documentation figures (docs/gallery/*.py). Not part of the xthermo package.

Phase regions keep the colours that xthermo.plot_slice draws. This module only sets fonts, axes, a sequential
colour map for property fields and the accent colours used for annotations.
"""
import os

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

NAVY = "#0b2545"
INK = "#1d2d44"
MUTED = "#5c677d"
CORAL = "#e2583e"
TEAL = "#16808a"

# light aqua to deep navy, for densities and other fields that increase with depth
OCEAN = LinearSegmentedColormap.from_list("xthermo_ocean", ["#f1f9fb", "#a9d6e5", "#4ba3c3", "#1f6f9a", "#13315c",
                                                            "#0b2545"])

RC = {"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": INK, "axes.labelcolor": INK,
      "axes.titlecolor": NAVY, "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.linewidth": 0.8,
      "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
      "legend.frameon": False, "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 150}

CREDIT = "Computed with xThermo · Driesner & Heinrich (2007), Driesner (2007) · IAPS-84 water"


def apply():
    matplotlib.rcParams.update(RC)


def outdir(default):
    """Output directory: the first command-line argument, else $XTHERMO_FIGURE_DIR, else `default`."""
    import sys
    d = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("XTHERMO_FIGURE_DIR", default)
    os.makedirs(d, exist_ok=True)
    return d
