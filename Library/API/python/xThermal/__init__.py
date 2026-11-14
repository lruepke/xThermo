"""Compatibility package: the former name of xthermo.

"import xThermal", "from xThermal import H2ONaCl" and "import xThermal.H2ONaCl" give the modules of xthermo itself (the
same module objects, not copies). New code should use "import xthermo".
"""
import sys as _sys

import xthermo as _xthermo
from xthermo import __version__, H2O, NaCl, H2ONaCl

for _name in ("H2O", "NaCl", "H2ONaCl"):
    _sys.modules[__name__ + "." + _name] = getattr(_xthermo, _name)


def __getattr__(name):   # everything else, e.g. the plotting helpers, from xthermo
    return getattr(_xthermo, name)
