# Python package "xthermo" (SWIG modules H2O, NaCl, H2ONaCl) and the compatibility package "xThermal".
# Build tree: importable from ${PROJECT_BINARY_DIR}/python (PYTHONPATH=<build>/python).
# Install:    <prefix>/${XTHERMO_PYTHON_INSTALL_DIR}/{xthermo,xThermal}, install component "python". A wheel build
#             (scikit-build-core, pyproject.toml) installs only this component, at the root of the wheel.
# Each module links the library objects and GSL, so the extension modules have no xThermo or GSL runtime dependency.

find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(SWIG 4.0 REQUIRED COMPONENTS python)
include(UseSWIG)

if(SKBUILD)
  set(_python_install_default .)
else()
  set(_python_install_default API/python/cp${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR})
endif()
set(XTHERMO_PYTHON_INSTALL_DIR ${_python_install_default} CACHE STRING
    "Install directory (relative to the prefix) of the directory that contains the Python packages")
set(XTHERMO_PYTHON_DIR ${PROJECT_BINARY_DIR}/python)
set(XTHERMO_PACKAGE_DIR ${XTHERMO_PYTHON_DIR}/xthermo)

# files SWIG reads through %include, so that changing them regenerates the wrappers
file(GLOB XTHERMO_SWIG_DEPENDS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/Library/thermo/*.i ${PROJECT_SOURCE_DIR}/Library/*/*.h
     ${PROJECT_SOURCE_DIR}/Library/H2ONaCl/numpy.i)

foreach(Module H2O NaCl H2ONaCl)
  set(_i ${PROJECT_SOURCE_DIR}/Library/${Module}/${Module}.i)
  set(SWIG_MODULE_${Module}_python_EXTRA_DEPS ${XTHERMO_SWIG_DEPENDS})
  set_property(SOURCE ${_i} PROPERTY CPLUSPLUS ON)
  set_property(SOURCE ${_i} PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)
  swig_add_library(${Module}_python TYPE MODULE LANGUAGE python
                   OUTPUT_DIR ${XTHERMO_PACKAGE_DIR}
                   OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${Module}
                   SOURCES ${_i})
  # -doxygen turns the C++ documentation comments into Python docstrings
  set_property(TARGET ${Module}_python PROPERTY SWIG_COMPILE_OPTIONS -doxygen)
  set_target_properties(${Module}_python PROPERTIES OUTPUT_NAME ${Module} LIBRARY_OUTPUT_DIRECTORY ${XTHERMO_PACKAGE_DIR})
  target_link_libraries(${Module}_python PRIVATE xThermal_objects Python3::Module)
  target_compile_options(${Module}_python PRIVATE
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-writable-strings -Wno-deprecated-declarations>)
  install(TARGETS ${Module}_python LIBRARY DESTINATION ${XTHERMO_PYTHON_INSTALL_DIR}/xthermo COMPONENT python)
  install(FILES ${XTHERMO_PACKAGE_DIR}/${Module}.py DESTINATION ${XTHERMO_PYTHON_INSTALL_DIR}/xthermo COMPONENT python)
endforeach()

# pure-Python files of both packages (copied into the build tree, so that it is importable as is) and the version
foreach(_f xthermo/__init__.py xthermo/plotting.py xthermo/states.py xThermal/__init__.py)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_f} ${XTHERMO_PYTHON_DIR}/${_f} COPYONLY)
endforeach()
file(WRITE ${XTHERMO_PACKAGE_DIR}/_version.py "__version__ = \"${PROJECT_VERSION}\"\n")
install(FILES ${XTHERMO_PACKAGE_DIR}/__init__.py ${XTHERMO_PACKAGE_DIR}/plotting.py ${XTHERMO_PACKAGE_DIR}/states.py ${XTHERMO_PACKAGE_DIR}/_version.py
        DESTINATION ${XTHERMO_PYTHON_INSTALL_DIR}/xthermo COMPONENT python)
install(FILES ${XTHERMO_PYTHON_DIR}/xThermal/__init__.py DESTINATION ${XTHERMO_PYTHON_INSTALL_DIR}/xThermal COMPONENT python)
