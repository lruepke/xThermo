#!/usr/bin/env bash
# Build the xThermo WebAssembly module (Library/API/wasm) reproducibly:
#   1. GSL (the release and SHA-256 pinned in cmake/xThermoGSL.cmake), static, with emconfigure/emmake;
#   2. the module with the standalone CMake project in this directory;
#   3. checks: no Triangle symbol in the module; prints sizes and SHA-256 of the files.
# Requires Emscripten on PATH (source <emsdk>/emsdk_env.sh), version XTHERMO_EMSCRIPTEN (default 4.0.10, the pinned
# version; set XTHERMO_EMSCRIPTEN=any to accept another), curl or a local archive (XTHERMO_GSL_ARCHIVE), make.
#
#   bash Library/API/wasm/make_wasm.sh <out-dir>      # -> <out-dir>/dist/{xthermo_wasm.js,xthermo_wasm.wasm,xthermo.js}
#
# GSL is built once per <out-dir> (kept in <out-dir>/gsl). XTHERMO_WASM_COMMIT (default: git HEAD) is recorded in
# xt_build_info(). Two builds from the same commit with the same tool versions give byte-identical files.
set -euo pipefail
OUT="${1:?usage: make_wasm.sh <out-dir>}"
mkdir -p "$OUT" && OUT="$(cd "$OUT" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
JOBS="${JOBS:-4}"

want="${XTHERMO_EMSCRIPTEN:-4.0.10}"
have="$(emcc --version | head -1 | sed -E 's/.* ([0-9]+\.[0-9]+\.[0-9]+).*/\1/')"
if [ "$want" != any ] && [ "$have" != "$want" ]; then
  echo "make_wasm.sh: Emscripten $have found, $want required (XTHERMO_EMSCRIPTEN=any to override)" >&2; exit 1
fi

# GSL release and checksum from the native build's pin
GSL_RELEASE="$(sed -nE 's/^set\(XTHERMO_GSL_RELEASE ([0-9.]+)\)/\1/p' "$ROOT/cmake/xThermoGSL.cmake")"
GSL_SHA256="$(sed -nE 's/^set\(XTHERMO_GSL_SHA256 ([0-9a-f]{64})\)/\1/p' "$ROOT/cmake/xThermoGSL.cmake")"
GSL="$OUT/gsl"; GSL_PREFIX="$GSL/install"
if [ ! -f "$GSL_PREFIX/lib/libgsl.a" ]; then
  mkdir -p "$GSL"
  archive="${XTHERMO_GSL_ARCHIVE:-}"
  if [ -z "$archive" ]; then
    archive="$GSL/gsl-$GSL_RELEASE.tar.gz"
    [ -f "$archive" ] || curl -fsSL -o "$archive" "https://ftp.gnu.org/gnu/gsl/gsl-$GSL_RELEASE.tar.gz"
  fi
  sum="$( (sha256sum "$archive" 2>/dev/null || shasum -a 256 "$archive") | cut -d' ' -f1)"
  [ "$sum" = "$GSL_SHA256" ] || { echo "make_wasm.sh: GSL archive checksum $sum, expected $GSL_SHA256" >&2; exit 1; }
  rm -rf "$GSL/src" && mkdir -p "$GSL/src" && tar xzf "$archive" -C "$GSL/src" --strip-components 1
  ( cd "$GSL/src"
    emconfigure ./configure --prefix="$GSL_PREFIX" --disable-shared --enable-static CFLAGS=-O2 > "$GSL/configure.log" 2>&1
    emmake make -j "$JOBS" > "$GSL/make.log" 2>&1
    emmake make install > "$GSL/install.log" 2>&1 ) || { echo "make_wasm.sh: GSL build failed, logs in $GSL" >&2; exit 1; }
fi

COMMIT="${XTHERMO_WASM_COMMIT:-$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)}"
emcmake cmake -S "$HERE" -B "$OUT/build" -DGSL_PREFIX="$GSL_PREFIX" -DXTHERMO_WASM_COMMIT="$COMMIT" \
  > "$OUT/configure.log" 2>&1 || { tail -30 "$OUT/configure.log"; exit 1; }
cmake --build "$OUT/build" -j "$JOBS" > "$OUT/compile.log" 2>&1 || { grep -E "error|undefined" "$OUT/compile.log" | head -30; exit 1; }

if grep -qiE "triangulate|triangleinit|incircle|delaunay" "$OUT/build/xthermo_wasm.js.symbols"; then
  echo "make_wasm.sh: Triangle symbols in the module" >&2; exit 1
fi
rm -rf "$OUT/dist" && mkdir -p "$OUT/dist/licenses"
cp "$OUT/build/xthermo_wasm.js" "$OUT/build/xthermo_wasm.wasm" "$OUT/build/xthermo.js" "$OUT/dist/"

# licence texts of everything in the module (GPL-3.0 section 4/6: the licence and notices travel with the module);
# taken from the sources that were built
EM="$(dirname "$(command -v emcc)")"
cp "$ROOT/LICENSE" "$OUT/dist/licenses/xThermo-MIT.txt"
cp "$GSL/src/COPYING" "$OUT/dist/licenses/GSL-GPL-3.0.txt"
cp "$ROOT/ThirdParties/PROST/COPYING.LIB-2.0" "$OUT/dist/licenses/PROST-LGPL-2.0.txt"
cp "$ROOT/ThirdParties/PROST/XTHERMO_MODIFICATIONS.txt" "$OUT/dist/licenses/PROST-modifications.txt"
cp "$EM/LICENSE" "$OUT/dist/licenses/Emscripten-MIT-or-NCSA.txt"
cp "$EM/system/lib/libc/musl/COPYRIGHT" "$OUT/dist/licenses/musl-MIT.txt"
cp "$EM/system/lib/libcxx/LICENSE.TXT" "$OUT/dist/licenses/LLVM-Apache-2.0-with-LLVM-exception.txt"
cat > "$OUT/dist/NOTICE.txt" <<EOF
xThermo WebAssembly module (xthermo_wasm.wasm, xthermo_wasm.js) and its JavaScript interface (xthermo.js)

xThermo $(sed -nE 's/^project\(xThermoWasm VERSION ([0-9.]+).*/\1/p' "$HERE/CMakeLists.txt"), commit $COMMIT
Built with Emscripten $have (emsdk) and GSL $GSL_RELEASE by Library/API/wasm/make_wasm.sh.

The module is a combined work. Because it contains the GNU Scientific Library (GPL-3.0-or-later), it is distributed
under the terms of the GNU General Public License, version 3 or (at your option) any later version
(licenses/GSL-GPL-3.0.txt). It contains:
  - xThermo (Zhikui Guo, Lars Rüpke, GEOMAR): MIT licence (licenses/xThermo-MIT.txt)
  - GNU Scientific Library $GSL_RELEASE (the GSL authors): GPL-3.0-or-later (licenses/GSL-GPL-3.0.txt)
  - PROST 4.1 (O. Bauer, O. Engel, TU Hamburg-Harburg), modified for xThermo: GNU Library General Public License 2.0
    (licenses/PROST-LGPL-2.0.txt; modifications: licenses/PROST-modifications.txt)
  - Emscripten runtime and JavaScript glue: MIT or University of Illinois/NCSA (licenses/Emscripten-MIT-or-NCSA.txt);
    musl libc: MIT (licenses/musl-MIT.txt); libc++, libc++abi, compiler-rt: Apache-2.0 with LLVM exceptions
    (licenses/LLVM-Apache-2.0-with-LLVM-exception.txt); dlmalloc: public domain
It does not contain Triangle, CoolProp or freesteam.

Corresponding source:
  - xThermo, including PROST and the build scripts: https://github.com/lruepke/xThermo, commit $COMMIT
    (Library/API/wasm/make_wasm.sh, Library/API/wasm/CMakeLists.txt)
  - GSL $GSL_RELEASE: https://ftp.gnu.org/gnu/gsl/gsl-$GSL_RELEASE.tar.gz, SHA-256 $GSL_SHA256
    (where this module is published, a copy of that archive is provided next to it as source/gsl-$GSL_RELEASE.tar.gz)
  - Emscripten $have: https://github.com/emscripten-core/emsdk (system libraries and JavaScript glue)
EOF
cd "$OUT/dist"
echo "xThermo WebAssembly module: commit $COMMIT, Emscripten $have, GSL $GSL_RELEASE"
for f in xthermo_wasm.wasm xthermo_wasm.js xthermo.js; do
  printf "%-18s %8d bytes  gzip %7d  sha256 %s\n" "$f" "$(wc -c < "$f")" "$(gzip -9cn "$f" | wc -c)" \
    "$( (sha256sum "$f" 2>/dev/null || shasum -a 256 "$f") | cut -d' ' -f1)"
done
