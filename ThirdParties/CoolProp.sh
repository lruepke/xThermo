cd coolprop
git apply ../CoolProp_patches/381c8535.patch

# figure out system name
#system_name=`uname -s`

# build only for the current host architecture
#ARCH="$(uname -m)"
#case "$ARCH" in
#  aarch64) ARCH=arm64 ;;
#  arm64) ARCH=arm64 ;;
#  x86_64|amd64) ARCH=x86_64 ;;
#  i386|i686) ARCH=x86 ;;
#  *) ARCH="$ARCH" ;;
#esac

system_name="$(uname -s)"     # e.g. Linux
ARCH="$(uname -m)"            # e.g. aarch64
INSTALL_DEST="$(cd .. && pwd)/install/CoolProp_${system_name}_${ARCH}"

build_one () {
  local which="$1"   # static | shared
  rm -rf "build-${which}"
  mkdir -p "build-${which}"
  pushd "build-${which}" >/dev/null

  cmake .. \
    -DCOOLPROP_INSTALL_PREFIX="${INSTALL_DEST}" \
    -DCOOLPROP_STATIC_LIBRARY=$([ "$which" = static ] && echo ON || echo OFF) \
    -DCOOLPROP_SHARED_LIBRARY=$([ "$which" = shared ] && echo ON || echo OFF) \
    -DCOOLPROP_OBJECT_LIBRARY=OFF \
    -DBITNESS=NATIVE

  # Strip any rogue -m64 the project might inject (harmless if not present)
  grep -R --null -l -- '-m64' . | while IFS= read -r -d '' f; do
    sed -i 's/-m64//g' "$f"
  done

  make -j"$(nproc)"
  make install
  popd >/dev/null
}

build_one static
build_one shared

#popd >/dev/null

# create alternate-name symlink for cross-platform cases with M chips (linux vs macOS)
cd ../install || { echo "ERROR: install dir not found from $(pwd)" >&2; exit 1; }
INSTALL_DIR="$(pwd)"

case "$ARCH" in
  arm64)
    src="CoolProp_${system_name}_arm64"
    dst="CoolProp_${system_name}_aarch64"
    echo "DEBUG: ARCH=arm64; install_dir=${INSTALL_DIR}; will create symlink ${dst} -> ${src}" >&2
    if [ -e "${INSTALL_DIR}/${dst}" ]; then
      echo "DEBUG: ${dst} already exists; skipping" >&2
    else
      if [ -e "${INSTALL_DIR}/${src}" ]; then
        ln -sfn "${src}" "${dst}" && echo "Created symlink ${dst} -> ${src}" >&2 || echo "ERROR: ln failed" >&2
      else
        echo "WARN: source ${src} does not exist; creating dangling symlink ${dst} -> ${src}" >&2
        ln -sfn "${src}" "${dst}" || echo "ERROR: ln failed" >&2
      fi
    fi
    ;;
  aarch64)
    src="CoolProp_${system_name}_aarch64"
    dst="CoolProp_${system_name}_arm64"
    echo "DEBUG: ARCH=aarch64; install_dir=${INSTALL_DIR}; will create symlink ${dst} -> ${src}" >&2
    if [ -e "${INSTALL_DIR}/${dst}" ]; then
      echo "DEBUG: ${dst} already exists; skipping" >&2
    else
      if [ -e "${INSTALL_DIR}/${src}" ]; then
        ln -sfn "${src}" "${dst}" && echo "Created symlink ${dst} -> ${src}" >&2 || echo "ERROR: ln failed" >&2
      else
        echo "WARN: source ${src} does not exist; creating dangling symlink ${dst} -> ${src}" >&2
        ln -sfn "${src}" "${dst}" || echo "ERROR: ln failed" >&2
      fi
    fi
    ;;
  *)
    echo "DEBUG: ARCH=${ARCH}; no alternate-name symlink rule" >&2
    ;;
esac

cd - >/dev/null || true