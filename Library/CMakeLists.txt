# The xThermal C++ library. All sources are compiled once, into the object library xThermal_objects; the installed
# library (static by default, shared with BUILD_SHARED_LIBS), the command-line tool and the Python modules use those
# objects.

set(XTHERMO_SOURCE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/thermo
  ${CMAKE_CURRENT_SOURCE_DIR}/AMR_LUT
  ${CMAKE_CURRENT_SOURCE_DIR}/H2O
  ${CMAKE_CURRENT_SOURCE_DIR}/NaCl
  ${CMAKE_CURRENT_SOURCE_DIR}/H2ONaCl)
set(XTHERMO_PROST_DIR ${PROJECT_SOURCE_DIR}/ThirdParties/PROST)

file(GLOB XTHERMO_SOURCES CONFIGURE_DEPENDS
  thermo/*.cpp H2O/*.cpp NaCl/*.cpp H2ONaCl/*.cpp ${XTHERMO_PROST_DIR}/*.c)

# version header for the command-line tool (written to the build tree)
set(xThermal_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(xThermal_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(xThermal_PATCH_VERSION ${PROJECT_VERSION_PATCH})
set(xThermal_VERSION ${PROJECT_VERSION})
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.git)
  execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE _git_hash OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_git_hash)
    set(xThermal_EXTRA_VERSION "-git-${_git_hash}")
    string(APPEND xThermal_VERSION ${xThermal_EXTRA_VERSION})
  endif()
endif()
string(TIMESTAMP xThermal_DATE %Y-%m-%d)
set(xThermal_GSL_VERSION ${XTHERMO_GSL_VERSION})
configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/generated/config.h)

# common compile settings of the library objects and their users
add_library(xThermal_settings INTERFACE)
foreach(_dir ${PROJECT_BINARY_DIR}/generated ${XTHERMO_SOURCE_DIRS})
  target_include_directories(xThermal_settings INTERFACE $<BUILD_INTERFACE:${_dir}>)
endforeach()
target_include_directories(xThermal_settings INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/xThermo>)
target_link_libraries(xThermal_settings INTERFACE GSL::gsl)
if(XTHERMO_OPENMP)
  target_compile_definitions(xThermal_settings INTERFACE USE_OMP=1)
  target_link_libraries(xThermal_settings INTERFACE OpenMP::OpenMP_CXX)
endif()

add_library(xThermal_objects OBJECT ${XTHERMO_SOURCES})
set_target_properties(xThermal_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(xThermal_objects PRIVATE ${XTHERMO_PROST_DIR})
target_link_libraries(xThermal_objects PUBLIC xThermal_settings)
# the code converts string literals to char* in many places; clang warns about every one
target_compile_options(xThermal_objects PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-writable-strings>)
if(TARGET gsl_bundled)
  add_dependencies(xThermal_objects gsl_bundled)
endif()

add_library(xThermal $<TARGET_OBJECTS:xThermal_objects>)
add_library(xThermo::xThermal ALIAS xThermal)
target_link_libraries(xThermal PUBLIC xThermal_settings)

# ---------------------------------------------------------------------------------------------------------------------
# installation and the CMake package for downstream projects: find_package(xThermo) -> xThermo::xThermal
include(CMakePackageConfigHelpers)
set(XTHERMO_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/xThermo)

install(TARGETS xThermal xThermal_settings EXPORT xThermoTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(EXPORT xThermoTargets NAMESPACE xThermo:: DESTINATION ${XTHERMO_CMAKE_DIR})
file(GLOB XTHERMO_PUBLIC_HEADERS
  stdfunc.h thermo/*.h thermo/*.H thermo/*.hpp AMR_LUT/*.h AMR_LUT/*.H H2O/*.h NaCl/*.h H2ONaCl/*.h)
install(FILES ${XTHERMO_PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xThermo)
# PROST headers (IAPS-84 backend) for programs that call PROST directly, e.g. the review verification harnesses
install(DIRECTORY ${XTHERMO_PROST_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xThermo/PROST FILES_MATCHING PATTERN "*.h")

# a bundled GSL is installed with the library, so that downstream projects link the same static GSL
if(XTHERMO_GSL_BUNDLED_PREFIX)
  install(DIRECTORY ${XTHERMO_GSL_BUNDLED_PREFIX}/include/gsl DESTINATION ${CMAKE_INSTALL_LIBDIR}/xThermo/gsl/include)
  install(FILES ${XTHERMO_GSL_BUNDLED_PREFIX}/lib/libgsl.a ${XTHERMO_GSL_BUNDLED_PREFIX}/lib/libgslcblas.a
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/xThermo/gsl/lib)
endif()

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/xThermoConfig.cmake.in
  ${PROJECT_BINARY_DIR}/xThermoConfig.cmake INSTALL_DESTINATION ${XTHERMO_CMAKE_DIR})
write_basic_package_version_file(${PROJECT_BINARY_DIR}/xThermoConfigVersion.cmake COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/xThermoConfig.cmake ${PROJECT_BINARY_DIR}/xThermoConfigVersion.cmake
        DESTINATION ${XTHERMO_CMAKE_DIR})
