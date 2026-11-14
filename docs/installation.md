# Installation

xThermo is a C++ library with a Python interface. Python users build and install the package `xthermo` with pip;
C++ users build and install the library with CMake. Both compile the library from source, which takes a few
minutes.

## Python

### From source (current)

You need Python ≥ 3.10, a C/C++ compiler (gcc ≥ 11 or Apple clang) and `make`. pip takes CMake and SWIG from PyPI if
they are not installed, and builds GSL from the official, checksum-verified release.

```bash
git clone https://github.com/lruepke/xThermo.git
cd xThermo
pip install .                 # or: pip install ".[plotting]" to add matplotlib
```

Use any environment you like (venv, conda, …). NumPy is installed as a dependency; Matplotlib is only needed for the
plotting helpers.

Check the installation:

```python
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
s = xthermo.tpx(sw, 573.15, 300e5, 0.1)      # T [K], p [Pa], X [kg/kg]
print(xthermo.__version__, s.phase_name, s.Rho)
```

The former package name still works: `from xThermal import H2ONaCl` gives the same modules.

### From PyPI (planned)

!!! info "Not yet published"
    A package on PyPI is planned: prebuilt wheels for Linux and macOS, installed with

    ```text
    pip install xthermo
    ```

    It is **not yet published**, so this command does not work yet; install from source as above. The licensing of
    binary packages, which contain GSL (GPL) and Triangle (non-commercial licence), is being settled first; see
    [About](about.md#license).

## C++ library

### Dependencies

| Dependency | Required | Notes |
|---|---|---|
| C++ compiler (gcc ≥ 11 or Apple clang) | yes | the code is compiled as C++11 |
| CMake ≥ 3.18 | yes | |
| GSL 2.x | yes | from the system, or built from source by CMake (see below) |
| SWIG ≥ 4 and Python ≥ 3.10 with headers | for the Python modules | skip with `-DXTHERMO_PYTHON=OFF` |
| PROST 4.1 | bundled | the IAPS-84 water backend; compiled automatically |

GSL from your package manager:

```bash
brew install gsl              # macOS
sudo apt install libgsl-dev   # Debian/Ubuntu
conda install -c conda-forge gsl
```

For a GSL in a non-standard location add `-DGSL_ROOT_DIR=<prefix>`. Without a system GSL, add
`-DXTHERMO_BUNDLED_GSL=ON`: CMake then downloads the official GSL 2.8 release from ftp.gnu.org, checks its SHA-256
checksum and builds it. For offline builds, point `-DXTHERMO_GSL_ARCHIVE=<path to gsl-2.8.tar.gz>` to a local copy
of the same archive.

### Build and install

From the repository root:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j
cmake --install build
```

The build type is Release unless `-DCMAKE_BUILD_TYPE` is given. Options:

| Option | Default | Description |
|---|---|---|
| `XTHERMO_PYTHON` | ON | Build the Python modules (SWIG) |
| `XTHERMO_BUNDLED_GSL` | OFF | Build the pinned GSL release from source instead of using the system GSL |
| `XTHERMO_OPENMP` | OFF | Use OpenMP for batch and lookup-table evaluation |
| `BUILD_SHARED_LIBS` | OFF | Build the C++ library as a shared library |
| `XTHERMO_BUILD_TESTS` | OFF | Build the developer test programs (`ctest`) |

`cmake --install` creates:

```
install/
├── bin/xThermal                 # command-line tool
├── include/xThermo/             # C++ headers
├── lib/libxThermal.a            # C++ library (lib64/ on some Linux distributions)
├── lib/cmake/xThermo/           # CMake package: find_package(xThermo)
└── API/python/cp3.13/           # Python modules; the name follows the Python version used by CMake
    ├── xthermo/                 # __init__.py, H2O.py, NaCl.py, H2ONaCl.py, the compiled modules, ...
    └── xThermal/                # compatibility name
```

To use these Python modules without pip, set `PYTHONPATH=/path/to/xThermo/install/API/python/cp3.13`.

### Use it from your CMake project

```cmake
find_package(xThermo REQUIRED)
target_link_libraries(my_app PRIVATE xThermo::xThermal)
```

Point CMake to the installation with `-DCMAKE_PREFIX_PATH=/path/to/xThermo/install` (or `export xThermo_DIR=/path/to/xThermo/install`).
More in the [C++ API overview](cpp/index.md).

### Other language bindings

The Julia wrapper (`Library/H2O/H2O.cxx`, CxxWrap) and the Matlab MEX functions (`Library/H2ONaCl/Matlab/`) are kept
in the source tree but are not part of the current CMake build.

## Developer environment (optional)

The tools used for development (documentation, verification) come from a conda + poetry environment:
- **conda** (`environment.yaml`) creates the environment `py313_xthermo` with Python 3.13 and poetry;
- **poetry** (`pyproject.toml`, `poetry.lock`) installs the Python packages at pinned versions.

```bash
conda env create -f environment.yaml   # once
conda activate py313_xthermo
poetry install --no-root
mkdocs serve                           # this documentation, locally
```

## Troubleshooting

- **`ModuleNotFoundError: No module named 'xthermo'`:** run `pip install .` in the environment you use, or set
  `PYTHONPATH` to the `cp3.XX` directory that matches your Python version.
- **SWIG or Python version changed (e.g. after `brew upgrade`):** delete `build/CMakeCache.txt`, then configure and
  build again.
- **`GSL was not found`:** install GSL, pass `-DGSL_ROOT_DIR=<prefix>`, or use `-DXTHERMO_BUNDLED_GSL=ON`.
