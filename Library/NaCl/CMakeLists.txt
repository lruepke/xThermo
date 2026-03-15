if(WIN32)
  if(BUILD_SHARED)
    return()
  endif()
endif()

set(Module NaCl)

FILE(GLOB SC_FILES "*.cpp" "../H2O/*.cpp" "../thermo/*.cpp" ${PROST_SRC_FILES})

# ================== Add test ==================
if(UNIX AND NOT APPLE)
  # fix static link issue for linux system
  set(LINK_FLAGS "-Wl,--no-as-needed -ldl")
endif()
add_executable(test_${Module} test.cxx  ${SC_FILES})
target_link_libraries(test_${Module} ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp})
set_target_properties(test_${Module} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} LINK_FLAGS "${LINK_FLAGS}")
if(USE_COOLPROP)
  CopyCoolPropLib("test_${Module}")
endif()

# ================== build wrapper for other language, e.g. python, tcl =====================
file(GLOB API_SRC "*.cpp" "../thermo/*.cpp" "../H2O/*.cpp")
list(APPEND API_SRC ${PROST_SRC_FILES})

list(APPEND Source_depends_module ${API_SRC})
list(APPEND libs_depend ${PYTHON_LIBRARIES} ${GSL_LIBRARIES_STATIC}  ${CoolProp_LIBRARIES_STATIC} ${LIBomp})

if(NOT BUILD_MATLAB_SEPARATELY)
  # ------ Python ------------------
  if(Build_API_Python)
    Build_Python_API(${Module} "${Source_depends_module}" "${libs_depend}")
  endif()
  # ------ Julia ------------------
  # if(Build_API_Julia)
  #   Build_Julia_API(${Module} "${SC_FILES}" "${libs_depend}")
  # endif(Build_API_Julia)
endif()
