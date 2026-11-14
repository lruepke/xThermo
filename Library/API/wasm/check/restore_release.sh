#!/bin/bash
# Keep the last verified Explorer release when a new one could not be built or verified.
#   restore_release.sh <previous app dir> <out dir>
# <previous app dir> is the app/ directory of the currently deployed site (from the gh-pages branch). The release set
# is restored only as a whole and only if:
#   1. its SHA256SUMS lists the complete set (module, glue, interface, NOTICE, licence texts, GSL source) and every
#      listed file matches (sha256sum -c --strict);
#   2. check_release.sh passes on it for the commit named in its NOTICE.txt (build information, NOTICE, GSL archive,
#      page links), so binary, notices and corresponding source belong to the same verified build;
#   3. its JavaScript interface xthermo.js is identical to Library/API/wasm/xthermo.js of this revision, so the page
#      of this revision drives the module through the interface it was tested with.
# Only the files listed in SHA256SUMS (and SHA256SUMS itself) are copied to <out dir>; nothing else from the old
# site. Exit status 0: restored; 1: nothing restored (no previous release, or it failed a check), <out dir> left empty.
set -uo pipefail
PREV="${1:?usage: restore_release.sh <previous app dir> <out dir>}"; OUT="${2:?out dir}"
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
refuse() { echo "restore_release.sh: previous Explorer release not restored: $*" >&2; rm -rf "$OUT"; exit 1; }

[ -f "$PREV/SHA256SUMS" ] || refuse "no SHA256SUMS in $PREV (no previous release deployed)"
(cd "$PREV" && sha256sum -c --strict --quiet SHA256SUMS) || refuse "checksums do not match"
commit="$(sed -nE 's/^.*https:\/\/github\.com\/lruepke\/xThermo, commit ([0-9a-f]{40})$/\1/p' "$PREV/NOTICE.txt" 2>/dev/null | head -1)"
[ -n "$commit" ] || refuse "NOTICE.txt names no source commit"
cmp -s "$PREV/xthermo.js" "$ROOT/Library/API/wasm/xthermo.js" \
  || refuse "its interface xthermo.js differs from Library/API/wasm/xthermo.js of this revision"

rm -rf "$OUT" && mkdir -p "$OUT"
cp "$PREV/SHA256SUMS" "$OUT/"
while read -r _ f; do mkdir -p "$OUT/$(dirname "$f")" && cp "$PREV/$f" "$OUT/$f" || refuse "copy failed: $f"; done < "$PREV/SHA256SUMS"
bash "$ROOT/Library/API/wasm/check/check_release.sh" "$OUT" "$commit" || refuse "release checks failed for commit $commit"
(cd "$OUT" && sha256sum -c --strict --quiet SHA256SUMS) || refuse "checksums do not match after copying"
echo "previous verified Explorer release restored: commit $commit, $(grep ' xthermo_wasm.wasm$' "$OUT/SHA256SUMS")"
