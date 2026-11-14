#!/bin/bash
# Explorer release checks (W2c): the files served next to the page match the module that was verified and its
# corresponding source.
#   check_release.sh <package dir | deployed app URL> <commit> [--online]
# <package dir> is the output of package_release.sh; with a URL (e.g. https://lruepke.github.io/xThermo/app) the
# files are downloaded first, so the deployed binary itself is checked. Checks:
#   1. SHA256SUMS lists the module, interface, NOTICE, every licence text and the GSL archive, and all match
#   2. the module's build information (xt_build_info) and NOTICE.txt name <commit>; GSL and Emscripten as pinned
#   3. the GSL archive has the SHA-256 pinned in cmake/xThermoGSL.cmake and the release named by the module
#   4. every licence, notice and source link of the page (docs/app/index.html, explorer.js) has its file
#   5. --online: the source links of the page resolve on GitHub for <commit> (tree and make_wasm.sh)
# Prints PASS/FAIL per check; exit status 1 if any fails. Needs node, curl (URL or --online), sha256sum or shasum.
set -uo pipefail
SRC="${1:?usage: check_release.sh <package dir | URL> <commit> [--online]}"; COMMIT="${2:?commit}"; ONLINE="${3:-}"
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
FAILS=0
check() { if [ "$1" = 0 ]; then echo "PASS $2"; else echo "FAIL $2"; FAILS=$((FAILS + 1)); fi; }
sha() { (sha256sum "$1" 2>/dev/null || shasum -a 256 "$1") | cut -d' ' -f1; }
REL="$(sed -nE 's/^set\(XTHERMO_GSL_RELEASE ([0-9.]+)\)/\1/p' "$ROOT/cmake/xThermoGSL.cmake")"
GSL_SHA="$(sed -nE 's/^set\(XTHERMO_GSL_SHA256 ([0-9a-f]{64})\)/\1/p' "$ROOT/cmake/xThermoGSL.cmake")"

if [[ "$SRC" =~ ^https?:// ]]; then
  DIR="$(mktemp -d)"; trap 'rm -rf "$DIR"' EXIT
  curl -fsSL -o "$DIR/SHA256SUMS" "$SRC/SHA256SUMS" || { echo "FAIL download $SRC/SHA256SUMS"; exit 1; }
  while read -r _ f; do mkdir -p "$DIR/$(dirname "$f")"; curl -fsSL -o "$DIR/$f" "$SRC/$f" || echo "  download failed: $f"; done < "$DIR/SHA256SUMS"
  echo "downloaded $(wc -l < "$DIR/SHA256SUMS") files from $SRC"
else
  DIR="$SRC"
fi

# 1. checksums and completeness
need="xthermo_wasm.wasm xthermo_wasm.js xthermo.js NOTICE.txt licenses/xThermo-MIT.txt licenses/GSL-GPL-3.0.txt
      licenses/PROST-LGPL-2.0.txt licenses/PROST-modifications.txt licenses/Emscripten-MIT-or-NCSA.txt
      licenses/musl-MIT.txt licenses/LLVM-Apache-2.0-with-LLVM-exception.txt source/gsl-$REL.tar.gz"
bad=0
for f in $need; do grep -q "  $f\$" "$DIR/SHA256SUMS" || { echo "  not in SHA256SUMS: $f"; bad=1; }; done
while read -r h f; do [ -f "$DIR/$f" ] && [ "$(sha "$DIR/$f")" = "$h" ] || { echo "  mismatch: $f"; bad=1; }; done < "$DIR/SHA256SUMS"
check $bad "files complete and identical to SHA256SUMS"
echo "  xthermo_wasm.wasm sha256 $(sha "$DIR/xthermo_wasm.wasm")"

# 2. build information and notice name the commit
info="$(cd "$DIR" && node -e "
  const X = require('./xthermo.js');
  X.load(require('./xthermo_wasm.js')).then(xt => console.log(JSON.stringify(xt.buildInfo)), e => { console.error(e); process.exit(1); });")"
echo "  build info: $info"
node -e "const b = JSON.parse(process.argv[1]); process.exit(b.commit === process.argv[2] && b.gsl === process.argv[3] && b.emscripten === '4.0.10' ? 0 : 1)" \
  "$info" "$COMMIT" "$REL"
check $? "module build information: commit $COMMIT, GSL $REL, Emscripten 4.0.10"
grep -q "commit $COMMIT\$" "$DIR/NOTICE.txt" && grep -q "https://github.com/lruepke/xThermo, commit $COMMIT\$" "$DIR/NOTICE.txt" \
  && grep -q "SHA-256 $GSL_SHA" "$DIR/NOTICE.txt" && grep -q "source/gsl-$REL.tar.gz" "$DIR/NOTICE.txt"
check $? "NOTICE.txt: commit, public source URL, GSL checksum and the shipped source archive"

# 3. GSL source archive
[ "$(sha "$DIR/source/gsl-$REL.tar.gz")" = "$GSL_SHA" ]
check $? "source/gsl-$REL.tar.gz has the pinned SHA-256"

# 4. links of the page
miss=0
for f in $(grep -oE 'href="(licenses/[A-Za-z0-9._-]+\.txt|NOTICE\.txt)"' "$ROOT/docs/app/index.html" | sed -E 's/href="|"//g' | sort -u); do
  [ -s "$DIR/$f" ] || { echo "  page links missing file: $f"; miss=1; }
done
grep -q 'source/gsl-\${encodeURIComponent(b.gsl)}.tar.gz' "$ROOT/docs/app/explorer.js" && [ -s "$DIR/source/gsl-$REL.tar.gz" ] || miss=1
check $miss "every licence, notice and source link of the page has its file"

# 5. public source of the commit
if [ "$ONLINE" = --online ]; then
  ok=0
  for u in "https://github.com/lruepke/xThermo/tree/$COMMIT" "https://github.com/lruepke/xThermo/blob/$COMMIT/Library/API/wasm/make_wasm.sh" \
           "https://github.com/lruepke/xThermo/blob/$COMMIT/ThirdParties/PROST/COPYING.LIB-2.0"; do
    code="$(curl -s -o /dev/null -w '%{http_code}' -L "$u")"; echo "  $code $u"; [ "$code" = 200 ] || ok=1
  done
  check $ok "public source of commit $COMMIT on GitHub (tree, build script, PROST)"
fi
echo "release checks failed: $FAILS"
exit $((FAILS > 0))
