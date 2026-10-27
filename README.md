# xThermo

Phase relations and thermodynamic properties of H₂O–NaCl fluids from 0 to 1000 °C and 1 to 5000 bar, after
Driesner & Heinrich (2007) and Driesner (2007): a C++ library for simulation codes, with a Python interface for NumPy
and Matplotlib. Developed at the GEOMAR Helmholtz Centre for Ocean Research Kiel.

![The H₂O–NaCl system at 250 bar, computed with xThermo](docs/assets/img/hero_desktop.png)

**Documentation, examples and the model's scope and limitations: https://lruepke.github.io/xThermo/**

## Install (Python, from source)

Requires Python ≥ 3.10, a C/C++ compiler and `make`:

```bash
git clone https://github.com/lruepke/xThermo.git
cd xThermo
pip install .                 # or: pip install ".[plotting]" to add matplotlib
```

A package on PyPI (`pip install xthermo`) is planned but **not yet published**.

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")                              # water backend: IAPS-84
s = xthermo.tpx(sw, np.array([573.15, 673.15]), 250e5, 0.032)  # T [K], p [Pa], X [kg/kg]
print(s.phase_name, s.Rho)
```

## C++

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j
cmake --install build
```

Then `find_package(xThermo)` and link `xThermo::xThermal`. Dependencies (GSL, optional SWIG) and options:
[installation guide](https://lruepke.github.io/xThermo/installation/).

## Citation

Please cite Driesner & Heinrich (2007), *Geochim. Cosmochim. Acta* 71, 4880–4901,
[doi:10.1016/j.gca.2006.01.033](https://doi.org/10.1016/j.gca.2006.01.033), and Driesner (2007), *Geochim.
Cosmochim. Acta* 71, 4902–4919, [doi:10.1016/j.gca.2007.05.026](https://doi.org/10.1016/j.gca.2007.05.026). Authors and citation of the software: [About](https://lruepke.github.io/xThermo/about/).

## Licence

xThermo's own source code is MIT-licensed (`LICENSE`). Builds also contain third-party components under their own
licences: PROST (GNU LGPL v2.0), GSL (GNU GPL v3 or later, linked statically) and Triangle (non-commercial use only).
A compiled xThermo is therefore not covered by the MIT licence alone. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
