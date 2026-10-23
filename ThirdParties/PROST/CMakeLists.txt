# 0. CMake Minimum version
cmake_minimum_required(VERSION 2.8...3.12 FATAL_ERROR)


project(PROST)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../install" CACHE PATH "..." FORCE)
endif()

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/*.c" "${PROJECT_SOURCE_DIR}/*.h")

list(APPEND buildTypes STATIC SHARED) # STATIC SHARED
foreach(buildType IN LISTS buildTypes)
  add_library ("${PROJECT_NAME}_${buildType}" ${buildType} ${SC_FILES} )
  set_target_properties("${PROJECT_NAME}_${buildType}" PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  # set proper lib type of GSL according to the build type
  if(${buildType} STREQUAL "SHARED")
    target_link_libraries("${PROJECT_NAME}_${buildType}" ${GSL_LIBRARIES_SHARED} ${LIBS_REQUIRED} ${LIBomp})
  elseif(${buildType} STREQUAL "STATIC")
    target_link_libraries("${PROJECT_NAME}_${buildType}" ${GSL_LIBRARIES_STATIC} ${LIBS_REQUIRED} ${LIBomp})
  else()
    message(FATAL_ERROR "buildType=${buildType} is neither STATIC nor SHARED, please set correct buildType: list(APPEND buildTypes SHARED)")
  endif()
  install (TARGETS "${PROJECT_NAME}_${buildType}" DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/lib)
endforeach()
# add_library (${PROJECT_NAME} STATIC ${SC_FILES})
# set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".so")

# install (TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/lib)
# INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/include
#           FILES_MATCHING PATTERN "*.h")
FILE(GLOB HEAD_FILES "${PROJECT_SOURCE_DIR}/*.h" )
install(FILES ${HEAD_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/include)