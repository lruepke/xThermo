cd coolprop
git apply ../CoolProp_patches/381c8535.patch

# figure out system name
system_name=`uname -s`

# build only for the current host architecture
ARCH="$(uname -m)"
case "$ARCH" in
  aarch64) ARCH=arm64 ;;
  arm64) ARCH=arm64 ;;
  x86_64|amd64) ARCH=x86_64 ;;
  i386|i686) ARCH=x86 ;;
  *) ARCH="$ARCH" ;;
esac

BUILD_DIR="build"
mkdir -p "${BUILD_DIR}"
pushd "${BUILD_DIR}" >/dev/null

# configure + build static
cmake  \
      -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp_${system_name}_${ARCH} \
      -DCOOLPROP_STATIC_LIBRARY=ON ..
make -j "$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)"
make install

# clean and build shared
rm -rf ./*
cmake  \
      -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp_${system_name}_${ARCH} \
      -DCOOLPROP_SHARED_LIBRARY=ON ..
make -j "$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)"
make install

popd >/dev/null

# create alternate-name symlink for cross-platform cases with M chips (linux vs macOS)
cd ../../install || exit 0
case "$ARCH" in
  arm64)
    [ -d "CoolProp_${system_name}_aarch64" ] || ln -s "CoolProp_${system_name}_arm64" "CoolProp_${system_name}_aarch64" || true
    ;;
  aarch64)
    [ -d "CoolProp_${system_name}_arm64" ] || ln -s "CoolProp_${system_name}_aarch64" "CoolProp_${system_name}_arm64" || true
    ;;
esac

cd - >/dev/null || true