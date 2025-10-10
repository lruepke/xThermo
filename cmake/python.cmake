if(NOT Build_API_Python)
    return()
endif(NOT Build_API_Python)

# ensure SWIG executable is known in the cache BEFORE CMake probes the SWIG language
if(NOT CMAKE_SWIG_COMPILER)
  find_program(_SWIG_EXECUTABLE swig HINTS /opt/homebrew/bin /usr/local/bin /usr/bin)
  if(_SWIG_EXECUTABLE)
    set(CMAKE_SWIG_COMPILER "${_SWIG_EXECUTABLE}" CACHE FILEPATH "SWIG executable")
    set(CMAKE_SWIG_COMPILER_ENV_VAR "SWIG_EXECUTABLE" CACHE STRING "Environment variable name for SWIG executable")
    message(STATUS "Using SWIG executable: ${CMAKE_SWIG_COMPILER}")
  endif()
endif()

#make sure SWIG is available (now the cache variables exist)
find_package(SWIG REQUIRED)
include(UseSWIG)
#enable_language(SWIG) # LHR: this may be a problem


# ---------------------------- Find Python 3 ---------------------------------
# Development.Module is useful, while CMake v3.16 doesn't have this feature, but CMake v3.22 works. Doesn't test on other versions.
set(numpy )
if (USE_NUMPY)
#    include_directories(/Users/zguo/miniconda3/lib/python3.9/site-packages/numpy/core/include)
    set(numpy NumPy)
endif (USE_NUMPY)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development Development.Module ${numpy})
#Check if numpy is found
if (Python3_NumPy_FOUND)
    message(STATUS "Numpy found: " ${Python3_NumPy_VERSION})
    message(STATUS "Numpy include dir: " ${Python3_NumPy_INCLUDE_DIRS})
    add_compile_definitions("HAVE_numpy=1")
endif (Python3_NumPy_FOUND)
list(APPEND CMAKE_SWIG_FLAGS "-py3" "-DPY3")
message(STATUS "Python include: " ${Python3_INCLUDE_DIRS})
message(STATUS "Python version: " ${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}.${Python3_VERSION_PATCH})
set(pythonTag cp${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR})
set(platTAG "none")
if(LINUX)
  set(platTAG "manylinux2010_x86_64")
elseif(APPLE)
  execute_process(COMMAND uname -m
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} ERROR_QUIET
    OUTPUT_VARIABLE arch_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(COMMAND system_profiler SPSoftwareDataType
                  COMMAND grep macOS
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} ERROR_QUIET
                  OUTPUT_VARIABLE macOS_version
                  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  string(REGEX MATCH "[0-9][0-9]\.[0-9]" macOS_version "${macOS_version}")
  string(REGEX REPLACE "([0-9]+)\\..*" "\\1" macOS_version_major "${macOS_version}")
  string(REGEX REPLACE "^.[0-9]+\\.([0-9]+).*" "\\1" macOS_version_minor "${macOS_version}")
  set(macOS_version ${macOS_version_major}_${macOS_version_minor})
  set(platTAG "macosx_${macOS_version}_${arch_name}")
  message(STATUS "platTAG: " ${platTAG})
elseif(WIN32)
  set(platTAG "win_amd64")
endif()

# ----------------------------------------------------------------------------
set(INSTALL_PATH_API_PYTHON ${CMAKE_INSTALL_PREFIX}/API/python/${pythonTag}/${PROJECT_NAME})
set(PYTHON_PROJECT ${PROJECT_NAME})
set(PYTHON_PROJECT_DIR ${PROJECT_BINARY_DIR}/python/${PROJECT_NAME})
function(Build_Python_API Module source_depend libs_depend)
    message(STATUS "Build Python module ${Module}")
    set_property(SOURCE ${Module}.i PROPERTY CPLUSPLUS ON)
    set_property(SOURCE ${Module}.i PROPERTY SWIG_MODULE_NAME ${Module}) # The name is the final generated .py file name
    # add dynamic library path
    set(Target_name "${Module}_Python")
    swig_add_library(${Target_name}
      TYPE SHARED
      LANGUAGE python
      OUTPUT_DIR ${PYTHON_PROJECT_DIR}/${Module}
      SOURCES ${Module}.i ${source_depend}  #Option1: only set .i file here and then link dependent libraries; Option2: set .i and all dependent source files, which are passed through the Build_Python_API function (cmake function, called in other place)
    )
    if(WIN32)
    else()
        target_compile_options(${Target_name} PRIVATE "-Wno-deprecated-declarations") #cancle warnings
    endif()
    target_link_libraries(${Target_name} PRIVATE ${libs_depend})
    target_include_directories(${Target_name} PRIVATE ${Python3_INCLUDE_DIRS} )
    if(APPLE)
      set_target_properties(${Target_name} PROPERTIES
        SUFFIX ".so"
        OUTPUT_NAME "${Module}"
        # INSTALL_RPATH "@loader_path;@loader_path/../../${PYTHON_PROJECT}/.libs"
        )
      # This configuration is important for Mac system.
      set_property(TARGET ${Target_name} APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
    elseif(UNIX)
        #   set_target_properties(${Target_name} PROPERTIES
        #     INSTALL_RPATH "$ORIGIN:$ORIGIN/../../${PYTHON_PROJECT}/.libs"
        #     )
    endif()
    if(MSVC)
      target_link_libraries(${Target_name} PRIVATE ${Python3_LIBRARIES})
    endif()
    # install
    install(TARGETS ${Target_name} DESTINATION ${INSTALL_PATH_API_PYTHON})
    install(FILES ${PYTHON_PROJECT_DIR}/${Module}/${Module}.py DESTINATION ${INSTALL_PATH_API_PYTHON})
    INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${Module}/ DESTINATION ${INSTALL_PATH_API_PYTHON}
      FILES_MATCHING PATTERN "*.py"
      PATTERN "IAPWS_Others" EXCLUDE
      )
    # copy benchmark test data from Driesner(2007)
    INSTALL(DIRECTORY 
    ${PROJECT_SOURCE_DIR}/Library/benchmark_test/Driesner2007a
      ${PROJECT_SOURCE_DIR}/Library/benchmark_test/Driesner2007b 
      DESTINATION ${INSTALL_PATH_API_PYTHON}/data
    )
endfunction(Build_Python_API )

#######################
## Python Packaging  ##
#######################
#file(MAKE_DIRECTORY python/${PYTHON_PROJECT})
#file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/API/python/__init__.py CONTENT "__version__ = \"${xThermal_MAJOR_VERSION}.${xThermal_MINOR_VERSION}.${xThermal_PATCH_VERSION}\"\n")
# foreach(Module IN LISTS Modules)
#   file(GENERATE OUTPUT ${PYTHON_PROJECT_DIR}/${Module}/__init__.py CONTENT "")
# endforeach()
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/API/python/__init__.py.in
    ${CMAKE_BINARY_DIR}/API/python/__init__.py
    @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/API/python/setup.py.in
  ${CMAKE_BINARY_DIR}/API/python/setup.py
  @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/API/python/setup.cfg.in
  ${CMAKE_BINARY_DIR}/API/python/setup.cfg
  @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/API/python/MANIFEST.in
  ${CMAKE_BINARY_DIR}/API/python/MANIFEST.in
  @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/API/python/publish.sh.in
  ${CMAKE_BINARY_DIR}/API/python/publish.sh
  @ONLY)

# install
install(FILES 
      ${CMAKE_BINARY_DIR}/API/python/setup.cfg
      ${CMAKE_BINARY_DIR}/API/python/setup.py
      ${CMAKE_CURRENT_SOURCE_DIR}/API/python/readme.rst
      ${CMAKE_BINARY_DIR}/API/python/MANIFEST.in
      ${CMAKE_BINARY_DIR}/API/python/publish.sh
      DESTINATION ${INSTALL_PATH_API_PYTHON}/../
      )
install(FILES ${CMAKE_BINARY_DIR}/API/python/__init__.py DESTINATION ${INSTALL_PATH_API_PYTHON})
# run command to publish python package
# python setup.py sdist
# twine upload dist/* --verbose
