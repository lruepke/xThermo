# Python Bindings — Installation and Usage

xThermo exposes its C++ equation-of-state library to Python via [SWIG](https://www.swig.org).
The resulting `xThermal` package provides three modules:

| Module | Class | Description |
|---|---|---|
| `H2O` | `cIAPS84`, `cIAPWS95_CoolProp` | Pure water (IAPWS-84 / IAPWS-95) |
| `NaCl` | `cNaCl` | Pure NaCl |
| `H2ONaCl` | `cH2ONaCl` | H₂O–NaCl mixture (Driesner 2007) |

---

## Prerequisites

| Dependency | macOS | Linux |
|---|---|---|
| C++ compiler | Xcode CLT (`xcode-select --install`) | `gcc` / `g++` |
| CMake ≥ 3.18 | `brew install cmake` | `apt install cmake` |
| SWIG ≥ 4 | `brew install swig` | `apt install swig` |
| conda / miniconda | [miniconda](https://docs.conda.io/en/latest/miniconda.html) | same |

Python itself and all Python packages are managed through conda + poetry as described below — do not install them separately.

---

## Setting up the Python environment

The project uses a two-layer approach:

- **conda** (`environment.yaml`) creates a base environment with Python 3.13, pip, and poetry
- **poetry** (`pyproject.toml` + `poetry.lock`) installs all Python package dependencies at pinned versions

```bash
# 1. Create the conda environment (only needed once)
conda env create -f environment.yaml

# 2. Activate it
conda activate py313_xthermo

# 3. Install Python package dependencies via poetry
poetry install
```

The conda environment is named `py313_xthermo` (defined in `environment.yaml`). The `poetry install` step reads `poetry.lock` to install exact pinned versions of all packages (NumPy, Matplotlib, Sphinx, etc.) into the active environment.

> **Updating dependencies:** if you change `pyproject.toml`, run `poetry update` to regenerate `poetry.lock`, then commit both files.

---

## Build and Install

With the conda environment active:

```bash
# 1. Clone and enter the repo
git clone <repo-url>
cd xThermo

# 2. Activate the environment
conda activate py313_xthermo

# 3. Configure
mkdir build && cd build
cmake .. \
  -DBuild_API_Python=ON \
  -DUSE_NUMPY=ON \
  -DCMAKE_INSTALL_PREFIX=../install

# 4. Compile and install
make -j$(nproc)
make install
```

After a successful build the install tree contains a ready-to-use Python package:

```
install/API/python/cp3.13/   # directory name matches your Python version
├── setup.py
├── setup.cfg
├── xThermal/
│   ├── __init__.py
│   ├── H2ONaCl/
│   │   ├── H2ONaCl.py      # SWIG-generated wrapper
│   │   └── _H2ONaCl.so     # compiled extension
│   ├── H2O/
│   │   ├── H2O.py
│   │   └── _H2O.so
│   └── NaCl/
│       ├── NaCl.py
│       └── _NaCl.so
```

### Option A — pip install into the conda environment (recommended)

```bash
cd ../install/API/python/cp3.13
pip install .
```

### Option B — add to PYTHONPATH (no install step)

```bash
export PYTHONPATH="$PWD/install/API/python/cp3.13:$PYTHONPATH"
```

---

## Quick-start Examples

### Pure water — H2O

```python
from xThermal import H2O

water = H2O.cIAPS84()

# Fluid constants
print(f"T range: {water.Tmin():.1f} – {water.Tmax():.1f} K")
print(f"p range: {water.pmin():.0f} – {water.pmax():.0f} Pa")
print(f"T_crit = {water.T_critical():.2f} K")
print(f"p_crit = {water.p_critical():.0f} Pa")
print(f"rho_crit = {water.rhomass_critical():.2f} kg/m³")

# Single-point density at 300 °C, 500 bar
water.UpdateState_TP(573.15, 500e5)
print(f"rho = {water.rhomass():.2f} kg/m³")
```

See full example: [`Library/H2O/test_H2O.py`](../Library/H2O/test_H2O.py)

---

### H₂O–NaCl mixture

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")   # or "IAPWS95"

# Single-point calculation: T=300 °C, p=300 bar, X=10 wt% NaCl
T = 573.15        # K
p = 300e5         # Pa
X = 0.1           # kg NaCl / kg solution

props = sw.UpdateState_TPX(T, p, X)

print(f"Phase:   {props.phase}")
print(f"Density: {props.Rho:.2f} kg/m³")
print(f"Enthalpy:{props.H/1e3:.2f} kJ/kg")
print(f"Cp:      {props.Cp:.2f} J/kg/K")
print(f"Mu:      {props.Mu:.4e} Pa·s")
```

See full example: [`Library/H2ONaCl/test_H2ONaCl.py`](../Library/H2ONaCl/test_H2ONaCl.py)

---
```

---

### Phase identification

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

# Phase region at a given T, p, X
T, p, X = 400.0, 250e5, 0.05
props = sw.UpdateState_TPX(T, p, X)
print(f"Phase region index: {props.phase}")

# Two-phase properties (VL region)
if props.phase == H2ONaCl.PhaseRegion_TwoPhase_VL:
    print(f"Liquid density:  {props.Rho_l:.2f} kg/m³")
    print(f"Vapour density:  {props.Rho_v:.2f} kg/m³")
    print(f"Liquid salinity: {props.X_l * 100:.2f} wt%")
    print(f"Vapour salinity: {props.X_v * 100:.2f} wt%")
```

---

### Phase boundaries

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

# Pressure of the V+L+H triple curve at T = 600 K
P_vlh = sw.P_VLH(600.0)
print(f"P_VLH(600 K) = {P_vlh/1e5:.1f} bar")

# Liquid and vapour salinity on the VL surface
T, p = 500.0, 200e5
X_l = sw.XL_VL(T, p)
X_v = sw.XV_VL(T, p)
print(f"X_liquid = {X_l*100:.2f} wt%,  X_vapour = {X_v*100:.4f} wt%")

# Critical point at a given pressure
T_crit = sw.T_Critical(p)
X_crit = sw.X_Critical(p)
print(f"T_crit = {T_crit:.1f} K,  X_crit = {X_crit*100:.2f} wt%")
```

See gallery example: [`doc/sphinx/source/gallery/H2ONaCl/pT/plot_VL.py`](sphinx/source/gallery/H2ONaCl/pT/plot_VL.py)

---

### Enthalpy-based flash (H–p–X)

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

H = 1.5e6    # J/kg
p = 300e5    # Pa
X = 0.1      # kg/kg

props = sw.UpdateState_HPX(H, p, X)
print(f"Temperature: {props.T - 273.15:.1f} °C")
print(f"Phase:       {props.phase}")
```

See gallery example: [`doc/sphinx/source/gallery/H2ONaCl/HPX/plot_PhaseDiagram_HPX.py`](sphinx/source/gallery/H2ONaCl/HPX/plot_PhaseDiagram_HPX.py)

---

### Pure NaCl

```python
from xThermal import NaCl

halite = NaCl.cNaCl("IAPS84")

T = 1100.0   # K
p = 1000e5   # Pa

props = halite.UpdateState_TPX(T, p)
print(f"Density:  {props.Rho:.2f} kg/m³")
print(f"Enthalpy: {props.H/1e3:.2f} kJ/kg")
print(f"Cp:       {props.Cp:.2f} J/kg/K")

# Phase boundaries
print(f"Melting T at {p/1e5:.0f} bar: {halite.Melting_T(p) - 273.15:.1f} °C")
print(f"Boiling p at {T:.0f} K:       {halite.Boiling_p(T)/1e5:.1f} bar")
```

See gallery example: [`doc/sphinx/source/gallery/NaCl/pT/plot_PhaseDiagram_NaCl.py`](sphinx/source/gallery/NaCl/pT/plot_PhaseDiagram_NaCl.py)

---

## ThermodynamicProperties fields

Every `UpdateState_*` call returns a `ThermodynamicProperties` object:

| Field | Unit | Description |
|---|---|---|
| `T`, `p`, `X` | K, Pa, kg/kg | Input state |
| `phase` | — | Phase region index |
| `Rho` | kg/m³ | Bulk density |
| `H` | J/kg | Bulk specific enthalpy |
| `Cp` | J/kg/K | Bulk specific heat |
| `Mu` | Pa·s | Bulk dynamic viscosity |
| `Rho_l`, `Rho_v`, `Rho_h` | kg/m³ | Phase densities (multi-phase) |
| `H_l`, `H_v`, `H_h` | J/kg | Phase enthalpies (multi-phase) |
| `X_l`, `X_v` | kg/kg | Phase salinities (multi-phase) |
| `S_l`, `S_v`, `S_h` | — | Phase saturations (multi-phase) |
| `IsobaricExpansivity` | 1/K | Thermal expansion coefficient |
| `IsothermalCompressibility` | 1/Pa | Isothermal compressibility |

---

## Troubleshooting

**`ImportError: No module named 'xThermal'`**
Make sure the `py313_xthermo` conda environment is active and you either ran `pip install .` or set `PYTHONPATH` to the `cp3.XX` directory.

**SWIG version mismatch after `brew upgrade`**
Delete the CMake cache and reconfigure:
```bash
cd build
rm CMakeCache.txt
cmake .. -DBuild_API_Python=ON -DUSE_NUMPY=ON
make -j$(nproc) && make install
```

**`Invalid version` error during `pip install`**
The generated `setup.py` has a stale empty version. Re-run `cmake .. && make install` to regenerate it, then `pip install .` again.
