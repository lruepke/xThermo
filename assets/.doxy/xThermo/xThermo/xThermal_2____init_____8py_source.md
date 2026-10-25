

# File \_\_init\_\_.py

[**File List**](files.md) **>** [**API**](dir_4b630c9df1aec2cd2fe9acaf75644364.md) **>** [**python**](dir_a5d5f7c295cec99c8e29471acd4df0e8.md) **>** [**xThermal**](dir_c28a5a2dd7a8abb7b693848948dd4e9c.md) **>** [**\_\_init\_\_.py**](xThermal_2____init_____8py.md)

[Go to the documentation of this file](xThermal_2____init_____8py.md)


```Python
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
```


