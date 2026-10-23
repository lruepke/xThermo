# Installation

## Dependencies

| Dependency | Required | Notes |
|---|---|---|
| gcc ≥ 11 | yes | C++11 standard |
| CMake ≥ 3.18 | yes | |
| GSL | yes | Build from `ThirdParties/` — see below |
| SWIG ≥ 4 | for Python bindings | macOS: `brew install swig`; Linux: `apt install swig` |
| conda / miniconda | for Python bindings | manages the Python environment |

## 1. Build GSL

GSL must be built before xThermo. A build script is provided in `ThirdParties/`:

```bash
cd ThirdParties
bash gsl.sh
```

This downloads the GSL source, builds it, and installs it to `ThirdParties/install/gsl_<OS>_<arch>`.
See [ThirdParties/readme.md](https://github.com/lruepke/xThermo/blob/master/ThirdParties/readme.md) for manual build instructions.

## 2. Build xThermo

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc)
make install
```

Key CMake options:

| Option | Default | Description |
|---|---|---|
| `Build_API_Python` | ON | Build Python bindings via SWIG |
| `USE_NUMPY` | ON | Enable NumPy array support in Python bindings |
| `USE_GSL` | ON | Link against GSL (required for full functionality) |
| `USE_COOLPROP` | OFF | Link against CoolProp (optional comparison library) |
| `USE_OMP` | OFF | Enable OpenMP parallelisation |

## 3. Python bindings

The Python bindings require the conda + poetry environment:

```bash
# Create and activate the conda environment
conda env create -f environment.yaml
conda activate py313_xthermo

# Install Python package dependencies
poetry install

# Build with Python bindings enabled
mkdir build && cd build
cmake .. -DBuild_API_Python=ON -DUSE_NUMPY=ON -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc) && make install

# Install the xThermal Python package
cd ../install/API/python/cp3.13
pip install .
```

!!! note
    The Python package name is `xThermal`. Import it as `from xThermal import H2ONaCl`.

## 4. Set the install path (C++ use)

After building, export the install directory so downstream CMake projects can find xThermo:

```bash
export xThermo_DIR=<path_to_repo>/install
```
