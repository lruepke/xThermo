# xThermo

C++ implementation of the Driesner & Heinrich (2007) equation of state for the H₂O–NaCl system, with Python bindings via SWIG.

Created by Zhikui Guo (GEOMAR).

---

## What it does

xThermo calculates thermodynamic properties (density, enthalpy, heat capacity, viscosity, phase state) of:

- **Pure water** — IAPWS-84 and IAPWS-95 formulations
- **Pure NaCl** — solid, liquid, and vapour phases
- **H₂O–NaCl mixtures** — full phase diagram including liquid, vapour, halite, and three-phase regions

The valid range covers 0–1000 °C, 1–5000 bar, and 0–100 wt% NaCl.

---

## Building from source

### 1. Dependencies

| Dependency | Required | Notes |
|---|---|---|
| gcc ≥ 11 | yes | C++11 standard |
| CMake ≥ 3.18 | yes | |
| GSL | yes | Build from `ThirdParties/` — see [ThirdParties/readme.md](ThirdParties/readme.md) |
| SWIG ≥ 4 | for Python bindings | macOS: `brew install swig`; Linux: `apt install swig` |
| conda / miniconda | for Python bindings | manages the Python environment |
| CoolProp | no | optional comparison library |
| freesteam / PROST | no | optional comparison libraries |

### 2. Build GSL first

GSL must be built before xThermo. A build script is provided:

```bash
cd ThirdParties
bash gsl.sh
```

See [ThirdParties/readme.md](ThirdParties/readme.md) for details and manual build instructions.

### 3. Build xThermo

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc)
make install
```

Key CMake options (all default to OFF unless noted):

| Option | Default | Description |
|---|---|---|
| `Build_API_Python` | ON | Build Python bindings via SWIG |
| `USE_NUMPY` | ON | Enable NumPy support in Python bindings |
| `USE_GSL` | ON | Link against GSL (required for full functionality) |
| `USE_COOLPROP` | OFF | Link against CoolProp |
| `Build_IAPWS_Others` | OFF | Build freesteam and PROST comparison libraries |
| `USE_OMP` | OFF | Enable OpenMP parallelisation |

### 4. Set the install path

```bash
export xThermo_DIR=<path_to_repo>/install
```

---

## Python bindings

The Python bindings expose the full H₂O–NaCl EOS to Python as the `xThermal` package. For detailed installation and usage instructions, see **[doc/python_install.md](doc/python_install.md)**.

### Quick setup

```bash
# Create and activate the conda environment (Python 3.13 + poetry)
conda env create -f environment.yaml
conda activate py313_xthermo

# Install Python package dependencies
poetry install

# Build with Python bindings enabled
mkdir build && cd build
cmake .. -DBuild_API_Python=ON -DUSE_NUMPY=ON -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc) && make install

# Install the xThermal package
cd ../install/API/python/cp3.13
pip install .
```

### Usage

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
props = sw.UpdateState_TPX(573.15, 300e5, 0.1)  # T [K], p [Pa], X [kg/kg]
print(props.Rho, props.H, props.phase)
```

---

## Documentation

Build the HTML documentation (requires the conda environment and `poetry install`):

```bash
cd doc/sphinx
make html
```

---

## Julia bindings

Julia bindings are under development. To enable them:

```bash
cmake .. -DBuild_API_Julia=ON -DJlCxx_DIR=<path_to_JlCxx>/lib/cmake/JlCxx
```

The Julia `CxxWrap` package is required (`] add CxxWrap` in the Julia REPL).
