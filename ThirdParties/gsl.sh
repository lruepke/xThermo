# download and install GSL
GSL_VER="${GSL_VER:-2-7}"
MIRROR_URL="https://github.com/lruepke/xThermo/releases/download/gsl-release-2-7/gsl-release-${GSL_VER}.tar.gz"
UPSTREAM_URL="https://cgit.git.savannah.gnu.org/cgit/gsl.git/snapshot/gsl-release${GSL_VER}.tar.gz"

download() {
  url="$1"
  echo "Trying $url"
  curl -L --progress-bar --fail --retry 5 --retry-delay 3 -o "gsl-${GSL_VER}.tar.gz" "$url"
}

download "$MIRROR_URL" || download "$UPSTREAM_URL"

# extract and rename to gsl
tar -xzf gsl-${GSL_VER}.tar.gz
mv gsl-release-${GSL_VER} gsl

#figure out architecture
ARCH="$(uname -m)"
#case "$ARCH" in
#  aarch64) ARCH=arm64 ;;
#  arm64) ARCH=arm64 ;;
#  x86_64|amd64) ARCH=x86_64 ;;
#  i386|i686) ARCH=x86 ;;
#  *) ARCH="$ARCH" ;;
#esac

PREFIX="${PWD}/../../install/gsl_$(uname -s)_${ARCH}"
echo "Installing GSL to $PREFIX"

# build and install
cd gsl
./autogen.sh
mkdir build
cd build
../configure --prefix=${PWD}/../../install/gsl_`uname -s`_`arch`
# add -fPIC to cblas compile flag
sed -i '/mode=compile/a\        -fPIC \\' cblas/Makefile
# make -j 8
make
make install

# create alternate-name symlink for cross-platform cases with M chips (linux vs macOS)
cd ../../install || { echo "ERROR: install dir not found from $(pwd)" >&2; exit 1; }
INSTALL_DIR="$(pwd)"
TARGET_OS=$(uname -s)

case "$ARCH" in
  arm64)
    src="gsl_${TARGET_OS}_arm64"
    dst="gsl_${TARGET_OS}_aarch64"
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
    src="gsl_${TARGET_OS}_aarch64"
    dst="gsl_${TARGET_OS}_arm64"
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