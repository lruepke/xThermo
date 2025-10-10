

# set(CMAKE_MACOSX_RPATH 1)
# set(CMAKE_INSTALL_RPATH 1)
# set(CMAKE_SKIP_RPATH TRUE) #do not handle rpath when installing
# set build types
if(WIN32)
  # MSVC has some issues for linking gsl, coolprop dll, so only provide static library.
  if(BUILD_SHARED)
    list(APPEND buildTypes SHARED)
  else()
    list(APPEND buildTypes STATIC)
  endif()
else()
  list(APPEND buildTypes STATIC SHARED)
endif()
list(APPEND Modules thermo H2O NaCl H2ONaCl)
# =====================================================================================
# global configuration
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

# =====================================================================================
#             GSL and openmp and some other configurations
# -------------------------------------------------------------------------------------
#include(gsl)
##include(openmp)
#include(apple)
#include(PROST)
#include(CoolProp)
##include(xThermal)
# =====================================================================================
#                            API generation
# -------------------------------------------------------------------------------------
include(swig)
include(python)
include(julia)
include(tcl)
# =====================================================================================
if(WIN32)
  add_compile_definitions("NO_TIMER=1")
  add_compile_definitions("FMT_USE_WINDOWS_H=0")
  add_compile_definitions("xTHERMO_DLL=1")
  add_compile_definitions("xThermal_DLL_EXPORT=1")
endif()

#2. Source code dir
aux_source_directory(. SRC_MAIN)
include_directories(.)
include_directories(AMR_LUT)
INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AMR_LUT/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")

# add submodules
list(APPEND SC_FILES_ALLMODULES "${PROST_SRC_FILES}")
foreach(Module IN LISTS Modules)
  include_directories(${Module})
  add_subdirectory(${Module})
  FILE(GLOB SC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${Module}/*.cpp")
  list(APPEND SC_FILES_ALLMODULES ${SC_FILE})
  # head files for installation
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${Module}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include
          FILES_MATCHING PATTERN "*.h")
endforeach()

# package all source code to xThermal lib
foreach(buildType IN LISTS buildTypes)
  set(CMAKE_LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/${buildType}" )
  set(LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${buildType}")
  add_library ("${LIB_NAME}_${buildType}" ${buildType} ${SC_FILES_ALLMODULES} )
  set_target_properties("${LIB_NAME}_${buildType}" PROPERTIES
          OUTPUT_NAME ${LIB_NAME}
          LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildType}
          ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildType}
          )
  # set proper lib type of GSL according to the build type
  if(${buildType} STREQUAL "SHARED")
    target_link_libraries("${LIB_NAME}_${buildType}" ${GSL_LIBRARIES_SHARED} ${CoolProp_LIBRARIES_SHARED} ${LIBomp})
  elseif(${buildType} STREQUAL "STATIC")
    target_link_libraries("${LIB_NAME}_${buildType}" ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp})
  else()
    message(FATAL_ERROR "buildType=${buildType} is neither STATIC nor SHARED, please set correct buildType: list(APPEND buildTypes SHARED)")
  endif()
  install (TARGETS "${LIB_NAME}_${buildType}" DESTINATION ${CMAKE_INSTALL_PREFIX}/${buildType})
  
endforeach()
# Install head files and third parties
FILE(GLOB HeadFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${HeadFiles} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
FILE(GLOB HeadFiles_PROST "${PROST_DIR}/*.h")
INSTALL(FILES ${HeadFiles_PROST} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/PROST)

# testing
# add_subdirectory(benchmark_test)
 add_subdirectory(test)

# add matlab mex
# matlab_add_mex(NAME lookup SHARED SRC ./API/Matlab/lookup.cpp LINK_TO ${GSL_LIBRARIES_STATIC} ${LIB_NAME}_STATIC ${Matlab_MEX_BLAS})
# add_test(NAME C_MEX_BLAS COMMAND ${Matlab_MAIN_PROGRAM} -batch "addpath('$<TARGET_FILE_DIR:matrixMultiply>'), r=runtests('${CMAKE_CURRENT_SOURCE_DIR}', Name='TestMex/test_blas'); assertSuccess(r)")

# ===================================================================================================
message(STATUS "")
message(STATUS "xThermal ${xThermal_VERSION} has been configured for ${CMAKE_SYSTEM_NAME}")
message(STATUS "")
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "C compiler: " ${CMAKE_C_COMPILER})
message(STATUS "C++ compiler: " ${CMAKE_CXX_COMPILER})
if(USE_OMP)
  message(STATUS "Use OpenMP: " ${LIBomp})
endif(USE_OMP)

message(STATUS "Build API: " ${Build_API_MultiLanguage})
message(STATUS "Build Python API: " ${Build_API_Python})
message(STATUS "Build Julia API: " ${Build_API_Julia})
message(STATUS "Build TCL API: " ${Build_API_tcl})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
# ===================================================================================================


