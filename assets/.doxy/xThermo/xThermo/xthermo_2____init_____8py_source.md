

# File \_\_init\_\_.py

[**File List**](files.md) **>** [**API**](dir_4b630c9df1aec2cd2fe9acaf75644364.md) **>** [**python**](dir_a5d5f7c295cec99c8e29471acd4df0e8.md) **>** [**xthermo**](dir_c7c9957b02953ad15ede12b0595b5fff.md) **>** [**\_\_init\_\_.py**](xthermo_2____init_____8py.md)

[Go to the documentation of this file](xthermo_2____init_____8py.md)


```Python
"""xthermo: equation of state of H2O-NaCl (Driesner & Heinrich 2007; Driesner 2007), with pure-water (IAPS-84, IAPWS-95)
and pure-NaCl models.

Modules: xthermo.H2ONaCl (cH2ONaCl), xthermo.H2O (cIAPS84, cIAPWS95), xthermo.NaCl (cNaCl). Importing the package
loads these compiled modules and nothing else: no model is constructed and matplotlib is not imported. The plotting
helpers are in xthermo.plotting (requires matplotlib); xthermo.plot_Pslice etc. import it on first use.

NumPy interface: xthermo.tpx(model, T, p, X) and xthermo.hpx(model, H, p, X) evaluate states on broadcast arrays and
return xthermo.States (see xthermo.states); xthermo.Status and the EvaluationError classes belong to it. They are
imported on first use, so "import xthermo" does not import numpy.

The former package name xThermal remains available: "from xThermal import H2ONaCl" gives the same modules.
"""
from ._version import __version__
# H2ONaCl must be loaded first. The three SWIG modules share one run-time type table, and if H2O or NaCl registers its
# types first, element access of H2ONaCl's nested vectors (e.g. PhaseBoundary_VLH_DeformLinear().X[3]) fails with a
# TypeError. The former package (xThermal) also loaded H2ONaCl first. The package is always initialised before any of
# its modules, so this order holds whatever the user imports.
from . import H2ONaCl
from . import H2O, NaCl

_PLOTTING = ("makeLogLinearAxes", "createWaterAxes", "plot_Pslice")
_STATES = ("tpx", "hpx", "States", "Status", "EvaluationError", "InvalidInputError", "OutOfRangeError",
           "NotConvergedError", "NotAvailableError")


def __getattr__(name):   # imported on first use (PEP 562): the plotting helpers (matplotlib) and the NumPy interface
    if name in _PLOTTING:
        from . import plotting
        return getattr(plotting, name)
    if name in _STATES:
        from . import states
        return getattr(states, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
```


