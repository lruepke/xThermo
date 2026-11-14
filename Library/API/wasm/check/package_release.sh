#!/bin/bash
# Explorer release package (W2c): the module files that passed ci_wasm.sh, unchanged, plus the corresponding GSL
# source, laid out as they are served next to the page (docs/app/):
#   xthermo_wasm.js, xthermo_wasm.wasm, xthermo.js, NOTICE.txt, licenses/   (ci_wasm.sh's dist, build a)
#   source/gsl-<release>.tar.gz                                              (the checksum-verified archive the
#                                                                             module was built from)
#   SHA256SUMS                                                               (every file above)
# The module is not rebuilt here. Usage: package_release.sh <ci_wasm OUT dir> <package dir>
set -euo pipefail
CI="${1:?usage: package_release.sh <ci_wasm OUT dir> <package dir>}"; PKG="${2:?package dir}"
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
REL="$(sed -nE 's/^set\(XTHERMO_GSL_RELEASE ([0-9.]+)\)/\1/p' "$ROOT/cmake/xThermoGSL.cmake")"
[ -f "$CI/dist/xthermo_wasm.wasm" ] || { echo "package_release.sh: no verified module in $CI/dist" >&2; exit 1; }
cmp -s "$CI/dist/xthermo_wasm.wasm" "$CI/wasm-a/dist/xthermo_wasm.wasm" \
  || { echo "package_release.sh: $CI/dist is not the verified build a" >&2; exit 1; }
rm -rf "$PKG" && mkdir -p "$PKG/source"
cp -R "$CI/dist/." "$PKG/"
cp "$CI/wasm-a/gsl/gsl-$REL.tar.gz" "$PKG/source/"
( cd "$PKG" && find . -type f ! -name SHA256SUMS | sed 's#^\./##' | LC_ALL=C sort \
    | while read -r f; do (sha256sum "$f" 2>/dev/null || shasum -a 256 "$f"); done > SHA256SUMS )
echo "release package $PKG:"; cat "$PKG/SHA256SUMS"
