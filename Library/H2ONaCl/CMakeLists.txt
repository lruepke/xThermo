if(WIN32)
  if(BUILD_SHARED)
    return()
  endif()
endif()

set(Module H2ONaCl)

FILE(GLOB SC_FILES "*.cpp" "../thermo/*.cpp" "../H2O/*.cpp" "../NaCl/*.cpp" ${PROST_SRC_FILES})

# ================== Add test ==================
if(UNIX AND NOT APPLE)
  # fix static link issue for linux system
  set(LINK_FLAGS "-Wl,--no-as-needed -ldl")
endif()
add_executable(test_${Module} test.cxx  ${SC_FILES})
target_link_libraries(test_${Module} ${GSL_LIBRARIES_SHARED} ${CoolProp_LIBRARIES_SHARED} ${LIBomp})
set_target_properties(test_${Module} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} LINK_FLAGS "${LINK_FLAGS}")
if(USE_COOLPROP)
  CopyCoolPropLib("test_${Module}")
endif()

# ================== build wrapper for other language, e.g. python, tcl =====================
file(GLOB API_SRC "*.cpp" "../thermo/*.cpp" "../H2O/*.cpp" "../NaCl/*.cpp")
list(APPEND API_SRC ${PROST_SRC_FILES})

list(APPEND Source_depends_module ${API_SRC})
list(APPEND libs_depend ${PYTHON_LIBRARIES} ${GSL_LIBRARIES_STATIC}  ${CoolProp_LIBRARIES_STATIC} ${LIBomp})

# ----- Matlab -------
if(Build_API_Matlab)
  list(APPEND Functions_Matlab prop_TPX prop_HPX prop_water_TP P_X_Critical X_VL lookup)
  # 同时构建静态和动态库两种
  foreach(buildType IN LISTS buildTypes)
    foreach(Func_Matlab IN LISTS Functions_Matlab)
      matlab_add_mex(NAME matlab_${Func_Matlab}_${buildType} SHARED SRC  Matlab/${Func_Matlab}.cpp LINK_TO "xThermal_${buildType}")
      set_target_properties("matlab_${Func_Matlab}_${buildType}" PROPERTIES
          OUTPUT_NAME ${Func_Matlab}
          LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildType}
          ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildType}
          )
      INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Matlab/ DESTINATION ${CMAKE_INSTALL_PREFIX}/API/Matlab/${buildType}
            FILES_MATCHING PATTERN "*.m")
      install (TARGETS matlab_${Func_Matlab}_${buildType} DESTINATION ${CMAKE_INSTALL_PREFIX}/API/Matlab/${buildType})
    endforeach()
  endforeach()
endif(Build_API_Matlab)

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

