# Third-party libraries

All third-party libraries are compiled and installed into `ThirdParties/install/`.
Build scripts are provided for GSL and CoolProp that automate the full process.
Run all scripts from the `ThirdParties/` directory.

---

## GSL — required

[GNU Scientific Library](https://www.gnu.org/software/gsl/) is used for multidimensional root-finding (e.g. saturation pressure, liquid/vapour density along the boiling curve of H₂O).

### Prerequisites

`automake` is required to generate the GSL build system:

```bash
# macOS
brew install automake

# Linux (Debian/Ubuntu)
sudo apt install automake
```

### Build with the script (recommended)

`gsl.sh` downloads the GSL source from a GitHub mirror, builds it, installs it to the correct platform-specific path, and creates the architecture symlink automatically:

```bash
cd ThirdParties
bash gsl.sh
```

The script installs to `ThirdParties/install/gsl_<OS>_<arch>` (e.g. `gsl_Darwin_arm64`) and also creates the alternate-name symlink needed for cross-platform compatibility.

### Build manually

```bash
cd ThirdParties/gsl
./autogen.sh
mkdir build && cd build
../configure --prefix=${PWD}/../../install/gsl_$(uname -s)_$(uname -m)
make -j$(nproc)
make install
```

Then create the architecture symlink manually (check your platform with `uname -sm`):

```bash
cd ThirdParties/install

# macOS Apple Silicon
ln -s gsl_Darwin_arm64 gsl_Darwin_aarch64

# Linux x86_64 (no symlink needed, arm64/aarch64 alias only applies to Apple)
```

---

## CoolProp — optional

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
for water. Its C sources in `ThirdParties/PROST/` are compiled directly into libxThermal (`cmake/PROST.cmake`,
`USE_PROST_SRC=ON` by default) and provide xThermo's `IAPS84` water backend; no separate build is needed. PROST is
licensed under the GNU LGPL v2.0 (`ThirdParties/PROST/COPYING.LIB-2.0`), not under xThermo's MIT licence; the changes
made for xThermo are listed in `ThirdParties/PROST/XTHERMO_MODIFICATIONS.txt`. The original distribution
(http://www.en.tu-harburg.de/PROST) is no longer available.

Only with `-DUSE_PROST_SRC=OFF` does xThermo link a separately built PROST library from `ThirdParties/install/PROST`:

```bash
cd ThirdParties/PROST
mkdir build && cd build
cmake ..
make -j$(nproc) && make install
```
