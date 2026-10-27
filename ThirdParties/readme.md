# Third-party libraries

## GSL — required

[GNU Scientific Library](https://www.gnu.org/software/gsl/) is used for root finding (Brent's method and the hybrid
multidimensional solver, e.g. for saturation pressures and phase boundaries). Any GSL 2.x works; xThermo uses only
long-stable parts of it.

CMake uses the system GSL by default (`brew install gsl`, `apt install libgsl-dev`, `conda install gsl`; a
non-standard location with `-DGSL_ROOT_DIR=<prefix>`). With `-DXTHERMO_BUNDLED_GSL=ON` it instead builds the official
GSL 2.8 release from source (`https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz`, SHA-256 pinned in
`cmake/xThermoGSL.cmake`) as a static library inside the build tree; `-DXTHERMO_GSL_ARCHIVE=<path>` uses a local copy
of that archive. No GSL is built or installed in `ThirdParties/`.

---

## CoolProp — optional

*Not part of the current CMake build: the CoolProp comparison class (`Library/H2O/IAPWS95_CoolProp.*`) is compiled
without CoolProp. The notes below describe how CoolProp was built for the former `USE_COOLPROP` option.*

[CoolProp](http://www.coolprop.org) is only used for validation and performance comparisons. It is disabled by default (`USE_COOLPROP=OFF`). Skip this section unless you specifically need it.

### Build with the script (recommended)

`CoolProp.sh` applies the required patch, builds both static and shared variants, installs to the correct platform-specific path, and creates the architecture symlink automatically:

```bash
cd ThirdParties
bash CoolProp.sh
```

The script installs to `ThirdParties/install/CoolProp_<OS>_<arch>`.

### What the patch does

The CoolProp source requires two small modifications before it can be compiled as part of xThermo:

1. **DLL export symbols** — `COOLPROP_VAR` must be defined in `DataStructures.h` so that low-level functions are accessible from shared libraries.
2. **C++11 detection** — `HAS_MOVE_SEMANTICS=1` must be added to `CMakeLists.txt`.

Both changes are bundled in `ThirdParties/CoolProp_patches/381c8535.patch` and applied automatically by the script.

### Build manually

```bash
cd ThirdParties/coolprop
git apply ../CoolProp_patches/381c8535.patch

INSTALL_DEST=$(cd .. && pwd)/install/CoolProp_$(uname -s)_$(uname -m)

# Static library
mkdir build-static && cd build-static
cmake .. -DCOOLPROP_INSTALL_PREFIX="$INSTALL_DEST" -DCOOLPROP_STATIC_LIBRARY=ON
make -j$(nproc) && make install
cd ..

# Shared library
mkdir build-shared && cd build-shared
cmake .. -DCOOLPROP_INSTALL_PREFIX="$INSTALL_DEST" -DCOOLPROP_SHARED_LIBRARY=ON
make -j$(nproc) && make install
cd ..

git reset --hard   # restore patched files
```

---

## freesteam — optional

*Not part of the current CMake build.*

[freesteam](http://freesteam.sourceforge.net) is only used for validation comparisons.

```bash
cd ThirdParties/freesteam-2.1
mkdir build && cd build
cmake ..
make -j$(nproc) && make install
```

---

## PROST — bundled, compiled automatically

PROST 4.1 (O. Bauer and O. Engel, TU Hamburg-Harburg, 1998) is an implementation of the IAPS-84 equation of state
for water. Its C sources in `ThirdParties/PROST/` are compiled directly into libxThermal and provide xThermo's
`IAPS84` water backend; no separate build is needed. PROST is licensed under the GNU LGPL v2.0
(`ThirdParties/PROST/COPYING.LIB-2.0`), not under xThermo's MIT licence; the changes made for xThermo are listed in
`ThirdParties/PROST/XTHERMO_MODIFICATIONS.txt`. The original distribution (http://www.en.tu-harburg.de/PROST) is no
longer available. `ThirdParties/PROST/CMakeLists.txt` (a standalone PROST library build) is not used by xThermo.
