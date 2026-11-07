# xThermo: H2O-NaCl equation of state (Driesner & Heinrich 2007; Driesner 2007), C++ library, command-line tool and
# Python bindings.
#
#   cmake -S . -B build                      # Release by default; GSL from the system
#   cmake --build build -j
#   cmake --install build --prefix <prefix>
#
# Options (all optional):
#   XTHERMO_PYTHON=ON|OFF        Python modules (SWIG); default ON
#   XTHERMO_BUNDLED_GSL=ON|OFF   build a pinned GSL release from source instead of using the system GSL; default OFF
#   XTHERMO_OPENMP=ON|OFF        OpenMP for batch and lookup-table evaluation; default OFF
#   XTHERMO_BUILD_TESTS=ON|OFF   developer test programs (Library/test, Library/*/test.cxx); default OFF
#   BUILD_SHARED_LIBS=ON|OFF     shared instead of static C++ library; default OFF
cmake_minimum_required(VERSION 3.18)
project(xThermo VERSION 1.0.0 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(XTHERMO_PYTHON "Build the Python modules (SWIG)" ON)
option(XTHERMO_BUNDLED_GSL "Build a pinned GSL release from source instead of using the system GSL" OFF)
option(XTHERMO_OPENMP "Use OpenMP for batch and lookup-table evaluation" OFF)
option(XTHERMO_BUILD_TESTS "Build the developer test programs" OFF)
option(BUILD_SHARED_LIBS "Build the C++ library as a shared library" OFF)

include(GNUInstallDirs)
include(xThermoGSL)

if(XTHERMO_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# installed programs and modules keep the run-time paths of dependencies outside the system directories (e.g. a GSL
# given with GSL_ROOT_DIR); a shared xThermo library is found relative to the executable
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
if(BUILD_SHARED_LIBS)
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()
endif()

add_subdirectory(Library)
add_subdirectory(Application/commandline)
if(XTHERMO_PYTHON)
  add_subdirectory(Library/API/python)
endif()
if(XTHERMO_BUILD_TESTS)
  enable_testing()
  add_subdirectory(Library/test)
endif()

message(STATUS "xThermo ${PROJECT_VERSION}: ${CMAKE_BUILD_TYPE}, ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}, "
               "GSL ${XTHERMO_GSL_VERSION} (${XTHERMO_GSL_ORIGIN}), Python ${XTHERMO_PYTHON}, OpenMP ${XTHERMO_OPENMP}, "
               "shared ${BUILD_SHARED_LIBS}")
