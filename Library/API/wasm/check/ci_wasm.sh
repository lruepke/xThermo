#!/bin/bash
# W2a CI job (explorer_wasm; review/W2A_WASM_MODULE.md): WebAssembly module build, reproducibility and
# native-versus-WebAssembly verification. Runs from the repository root in the Emscripten image (emsdk 4.0.10, with
# g++, make, cmake, Python 3 + NumPy and SWIG for the native module). Every check prints PASS or FAIL; the exit
# status is 0 only if all pass. Release readiness of the Explorer module = this job passing on the deployed commit.
#   1. native library and Python module, GSL from the same pinned source as the module (XTHERMO_BUNDLED_GSL)
#   2. module build twice in separate directories (GSL included): byte-identical files (reproducibility)
#   3. module contents: size limits, no Triangle symbol
#   4. native reference (native_reference.py; glibc libm FMA variants off) and comparison (verify.mjs)
#   5. Explorer page interface tests (test_explorer.mjs; W2b)
# Shared by the GitLab job explorer_wasm and the GitHub documentation workflow (job explorer, which deploys the
# module of step 2 that passed these checks). XTHERMO_WASM_COMMIT (environment) sets the commit recorded in the module.
# Overrides: OUT (default ci-wasm), JOBS (default 4), PYTHON (default python3).
set -uo pipefail
ROOT="$(pwd)"
OUT="${OUT:-$ROOT/ci-wasm}"
JOBS="${JOBS:-4}"
PYTHON="${PYTHON:-python3}"
mkdir -p "$OUT"
FAILS=0
check() { if [ "$1" = 0 ]; then echo "PASS $2"; else echo "FAIL $2"; FAILS=$((FAILS + 1)); fi; }
T_LAST=$(date +%s)
timing() { local now; now=$(date +%s); echo "TIMING $1 $((now - T_LAST))" | tee -a "$OUT/timings.txt"; T_LAST=$now; }

# 1. native. The bundled GSL (an ExternalProject running make -j8) is built first on its own: inside a parallel
#    outer make its install step fails with GNU make 4.3 ("write jobserver: Bad file descriptor"; W2a finding)
PYEXE="$($PYTHON -c 'import sys; print(sys.executable)')"
cmake -S "$ROOT" -B "$OUT/native" -DXTHERMO_BUNDLED_GSL=ON -DPython3_EXECUTABLE="$PYEXE" \
      -DCMAKE_INSTALL_PREFIX="$OUT/native-install" > "$OUT/native-configure.log" 2>&1 \
  && cmake --build "$OUT/native" --target gsl_bundled > "$OUT/native-gsl.log" 2>&1 \
  && cmake --build "$OUT/native" -j "$JOBS" > "$OUT/native-compile.log" 2>&1 \
  && cmake --install "$OUT/native" > "$OUT/native-install.log" 2>&1
rc=$?; check $rc "native build"
[ $rc = 0 ] || { tail -n 40 "$OUT"/native-*.log; exit 1; }
timing native

# 2. module, twice
for d in a b; do
  JOBS=$JOBS bash "$ROOT/Library/API/wasm/make_wasm.sh" "$OUT/wasm-$d" > "$OUT/wasm-$d.log" 2>&1
  rc=$?; check $rc "module build $d"; tail -n 4 "$OUT/wasm-$d.log"
  [ $rc = 0 ] || exit 1
done
timing module
diff -r "$OUT/wasm-a/dist" "$OUT/wasm-b/dist" > "$OUT/reproducibility.diff" 2>&1
check $? "reproducible: the two builds' files (module, glue, interface, NOTICE, licences) byte-identical"
missing=0
for f in NOTICE.txt licenses/GSL-GPL-3.0.txt licenses/PROST-LGPL-2.0.txt licenses/PROST-modifications.txt \
         licenses/xThermo-MIT.txt licenses/Emscripten-MIT-or-NCSA.txt licenses/musl-MIT.txt \
         licenses/LLVM-Apache-2.0-with-LLVM-exception.txt; do
  [ -s "$OUT/wasm-a/dist/$f" ] || { echo "  missing $f"; missing=1; }
done
grep -q "GNU GENERAL PUBLIC LICENSE" "$OUT/wasm-a/dist/licenses/GSL-GPL-3.0.txt" && grep -q "Version 3" "$OUT/wasm-a/dist/licenses/GSL-GPL-3.0.txt" || missing=1
check $missing "licence texts and NOTICE in the module directory"
cmp -s "$OUT/wasm-a/gsl/install/lib/libgsl.a" "$OUT/wasm-b/gsl/install/lib/libgsl.a"; check $? "reproducible: GSL archives byte-identical"

# 3. contents
size=$(wc -c < "$OUT/wasm-a/dist/xthermo_wasm.wasm"); gz=$(gzip -9cn "$OUT/wasm-a/dist/xthermo_wasm.wasm" | wc -c)
echo "  xthermo_wasm.wasm $size bytes, gzip $gz bytes"
[ "$size" -le 1048576 ] && [ "$gz" -le 393216 ]; check $? "size: wasm <= 1 MiB, gzip <= 384 KiB"
! grep -qiE "triangulate|triangleinit|incircle|delaunay" "$OUT/wasm-a/build/xthermo_wasm.js.symbols"; check $? "no Triangle symbol in the module"

# 4. verification
#    The native reference runs with glibc's FMA/AVX2 variants of exp, log, pow, sin, ... switched off
#    (GLIBC_TUNABLES). glibc selects them at run time on CPUs with FMA, so a native reference taken on such a runner
#    differs from one taken without FMA by up to 3e-7 relative in H-P-X (pipeline 1721672; W2a report section 4).
#    The module is the same on every CPU; with the generic code paths the reference is too (the compiled code is
#    baseline x86-64 either way). The FMA bit (CPUID.1:ECX bit 12) must be inactive under the tunable.
PYTAG="cp$($PYTHON -c 'import sys; print(f"{sys.version_info[0]}.{sys.version_info[1]}")')"
LDSO=/lib64/ld-linux-x86-64.so.2
NOFMA="glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4"
ecx() { env ${1:+GLIBC_TUNABLES=$1} "$LDSO" --list-diagnostics | sed -nE 's/^x86\.cpu_features\.features\[0x0\]\.active\[0x2\]=//p'; }
e0=$(ecx ""); e1=$(ecx "$NOFMA")
echo "  CPU: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-); FMA active: default $(( (e0 >> 12) & 1 )), reference $(( (e1 >> 12) & 1 ))"
[ -n "$e1" ] && [ $(( (e1 >> 12) & 1 )) = 0 ]; check $? "FMA variants of glibc libm off for the native reference"
XTHERMAL_PYTHONPATH="$OUT/native-install/API/python/$PYTAG" GLIBC_TUNABLES="$NOFMA" \
  $PYTHON "$ROOT/Library/API/wasm/check/native_reference.py" "$OUT/reference.json.gz" > "$OUT/reference.log" 2>&1
rc=$?; check $rc "native reference"; grep -v Warning "$OUT/reference.log" | tail -16
timing reference
if [ $rc = 0 ]; then
  node "$ROOT/Library/API/wasm/check/verify.mjs" "$OUT/wasm-a/dist" "$OUT/reference.json.gz" "$OUT/verification-report.txt" \
    > "$OUT/verify.log" 2> "$OUT/verify-stderr.log"
  rc=$?; sed -n '/^points/,$p' "$OUT/verification-report.txt" | head -80
  check $rc "native-versus-WebAssembly verification (report: verification-report.txt)"
fi
timing verification
# 5. Explorer page interface (docs/app, W2b) with this module
node "$ROOT/Library/API/wasm/check/test_explorer.mjs" "$OUT/wasm-a/dist" > "$OUT/explorer-tests.txt" 2>&1
rc=$?; tail -n 1 "$OUT/explorer-tests.txt"; grep '^FAIL' "$OUT/explorer-tests.txt"
check $rc "Explorer page interface tests (explorer-tests.txt)"
# context, not a check: on a runner with FMA, the same comparison against the native default (FMA) code paths
if [ $(( (e0 >> 12) & 1 )) = 1 ]; then
  XTHERMAL_PYTHONPATH="$OUT/native-install/API/python/$PYTAG" $PYTHON "$ROOT/Library/API/wasm/check/native_reference.py" \
    "$OUT/reference-fma.json.gz" > "$OUT/reference-fma.log" 2>&1 \
    && node "$ROOT/Library/API/wasm/check/verify.mjs" "$OUT/wasm-a/dist" "$OUT/reference-fma.json.gz" \
         "$OUT/context-native-fma.txt" > /dev/null 2>&1
  echo "  context (not a check): against the native FMA code paths: $(grep -E '^(points|RESULT|FAIL )' "$OUT/context-native-fma.txt" 2>/dev/null | tr '\n' ' ')"
  timing context
fi
rm -rf "$OUT/dist" && cp -R "$OUT/wasm-a/dist" "$OUT/dist"
echo "checks failed: $FAILS"
exit $((FAILS > 0))
