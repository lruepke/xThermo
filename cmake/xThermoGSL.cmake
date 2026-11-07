# GSL for xThermo. xThermo uses only the 1-D root finder (Brent), the multi-dimensional root finder (hybrids), vectors
# and gsl_strerror, which are available in every GSL release from 1.x on; no particular version is required.
#
# Default: a system GSL found by CMake's FindGSL (hint: -DGSL_ROOT_DIR=<prefix>).
# XTHERMO_BUNDLED_GSL=ON: build the official GSL release below from source, as a static, position-independent library.
# Requires a C compiler and make. XTHERMO_GSL_ARCHIVE may point to a local copy of the same archive (offline builds);
# the checksum is verified in either case.
# Defines the targets GSL::gsl and GSL::gslcblas and the variables XTHERMO_GSL_VERSION and XTHERMO_GSL_ORIGIN.

set(XTHERMO_GSL_RELEASE 2.8)
set(XTHERMO_GSL_SHA256 6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190)
set(XTHERMO_GSL_ARCHIVE "https://ftp.gnu.org/gnu/gsl/gsl-${XTHERMO_GSL_RELEASE}.tar.gz"
    CACHE STRING "URL or local path of the GSL source archive used with XTHERMO_BUNDLED_GSL")

if(NOT XTHERMO_BUNDLED_GSL)
  find_package(GSL)
  if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL was not found. Install it (e.g. 'brew install gsl', 'apt install libgsl-dev', "
                        "'conda install gsl'), point CMake to it with -DGSL_ROOT_DIR=<prefix>, "
                        "or build it from source with -DXTHERMO_BUNDLED_GSL=ON.")
  endif()
  set(XTHERMO_GSL_VERSION ${GSL_VERSION})
  set(XTHERMO_GSL_ORIGIN system)
  return()
endif()

include(ExternalProject)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)   # extracted files get the extraction time (CMake >= 3.24)
endif()
find_program(XTHERMO_MAKE NAMES gmake make REQUIRED)
set(_gsl_prefix ${CMAKE_BINARY_DIR}/gsl)
set(_gsl_cflags "-O2")
if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET)
  string(APPEND _gsl_cflags " -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
  foreach(_arch IN LISTS CMAKE_OSX_ARCHITECTURES)
    string(APPEND _gsl_cflags " -arch ${_arch}")
  endforeach()
endif()
ExternalProject_Add(gsl_bundled
  URL ${XTHERMO_GSL_ARCHIVE}
  URL_HASH SHA256=${XTHERMO_GSL_SHA256}
  PREFIX ${_gsl_prefix}
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-shared --enable-static --with-pic
                    CC=${CMAKE_C_COMPILER} CFLAGS=${_gsl_cflags}
  BUILD_COMMAND ${XTHERMO_MAKE} -j8
  INSTALL_COMMAND ${XTHERMO_MAKE} install
  BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libgsl.a <INSTALL_DIR>/lib/libgslcblas.a
  LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON)
ExternalProject_Get_Property(gsl_bundled INSTALL_DIR)
file(MAKE_DIRECTORY ${INSTALL_DIR}/include)   # must exist at configure time for the imported targets

add_library(GSL::gslcblas STATIC IMPORTED GLOBAL)
set_target_properties(GSL::gslcblas PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/lib/libgslcblas.a
                                               INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
add_library(GSL::gsl STATIC IMPORTED GLOBAL)
set_target_properties(GSL::gsl PROPERTIES IMPORTED_LOCATION ${INSTALL_DIR}/lib/libgsl.a
                                          INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include
                                          INTERFACE_LINK_LIBRARIES "GSL::gslcblas;$<$<NOT:$<PLATFORM_ID:Darwin>>:m>")
set(XTHERMO_GSL_BUNDLED_PREFIX ${INSTALL_DIR})
set(XTHERMO_GSL_VERSION ${XTHERMO_GSL_RELEASE})
set(XTHERMO_GSL_ORIGIN "bundled, built from source")
