# xThermo WebAssembly module (xThermo Explorer; review/W1_WASM_FEASIBILITY.md). A standalone project, configured
# only with Emscripten and independent of the top-level build; normally run through make_wasm.sh:
#
#   emcmake cmake -S Library/API/wasm -B <dir> -DGSL_PREFIX=<GSL built with Emscripten> [-DXTHERMO_WASM_COMMIT=<sha>]
#   cmake --build <dir>        # -> <dir>/xthermo_wasm.js, xthermo_wasm.wasm, xthermo.js
#
# Sources: the library's own source list (Library/CMakeLists.txt: thermo, H2O, NaCl, H2ONaCl, PROST) without
# Triangle (thermo/triangle.cpp, used only by the lookup-table code; not part of this module). The link fails on any
# reachable undefined symbol, so a module that needed Triangle could not be built. Same language standard and
# Release optimisation as the native build; C++ exceptions as WebAssembly exceptions (typed errors are caught).
cmake_minimum_required(VERSION 3.18)
project(xThermoWasm VERSION 1.0.0 LANGUAGES C CXX)

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "Library/API/wasm builds only with Emscripten (emcmake cmake ...); see make_wasm.sh")
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GSL_PREFIX "" CACHE PATH "Prefix of a static GSL built with Emscripten (include/, lib/libgsl.a, lib/libgslcblas.a)")
set(XTHERMO_WASM_COMMIT "unknown" CACHE STRING "xThermo commit recorded in xt_build_info()")
if(NOT EXISTS "${GSL_PREFIX}/lib/libgsl.a" OR NOT EXISTS "${GSL_PREFIX}/include/gsl/gsl_version.h")
  message(FATAL_ERROR "GSL_PREFIX='${GSL_PREFIX}' does not contain a GSL built with Emscripten (see make_wasm.sh)")
endif()

get_filename_component(XTHERMO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(L ${XTHERMO_ROOT}/Library)
file(GLOB XTHERMO_SOURCES CONFIGURE_DEPENDS
  ${L}/thermo/*.cpp ${L}/H2O/*.cpp ${L}/NaCl/*.cpp ${L}/H2ONaCl/*.cpp ${XTHERMO_ROOT}/ThirdParties/PROST/*.c)
list(FILTER XTHERMO_SOURCES EXCLUDE REGEX "/thermo/triangle\\.cpp$")

add_executable(xthermo_wasm ${XTHERMO_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/xthermo_wasm.cpp)
set_target_properties(xthermo_wasm PROPERTIES SUFFIX ".js")
target_include_directories(xthermo_wasm PRIVATE ${L} ${L}/thermo ${L}/AMR_LUT ${L}/H2O ${L}/NaCl ${L}/H2ONaCl
                                                ${XTHERMO_ROOT}/ThirdParties/PROST ${GSL_PREFIX}/include)
target_compile_definitions(xthermo_wasm PRIVATE XTHERMO_WASM_VERSION="${PROJECT_VERSION}"
                                                XTHERMO_WASM_COMMIT="${XTHERMO_WASM_COMMIT}")
target_compile_options(xthermo_wasm PRIVATE -fwasm-exceptions $<$<COMPILE_LANGUAGE:CXX>:-Wno-writable-strings>
                                            -ffile-prefix-map=${XTHERMO_ROOT}=xthermo)
target_link_libraries(xthermo_wasm PRIVATE ${GSL_PREFIX}/lib/libgsl.a ${GSL_PREFIX}/lib/libgslcblas.a)
set(XTHERMO_WASM_EXPORTS
  _xt_abi_version _xt_build_info _xt_new _xt_free _xt_last_error _xt_nfields _xt_field_name _xt_eval _xt_message
  _xt_phase_name _xt_status_name _malloc _free)
string(REPLACE ";" "," XTHERMO_WASM_EXPORTS "${XTHERMO_WASM_EXPORTS}")
target_link_options(xthermo_wasm PRIVATE
  -fwasm-exceptions
  -sMODULARIZE=1 -sEXPORT_NAME=createXThermoModule   # a factory; classic script, usable with importScripts in a Worker
  -sENVIRONMENT=web,worker,node
  -sALLOW_MEMORY_GROWTH=1
  -sFILESYSTEM=0
  -sERROR_ON_UNDEFINED_SYMBOLS=1
  -sEXPORTED_FUNCTIONS=${XTHERMO_WASM_EXPORTS}
  -sEXPORTED_RUNTIME_METHODS=UTF8ToString,HEAPF64,HEAP32
  --emit-symbol-map)

# the JavaScript interface next to the module
add_custom_command(TARGET xthermo_wasm POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/xthermo.js $<TARGET_FILE_DIR:xthermo_wasm>/xthermo.js)
