# 0. CMake Minimum version
cmake_minimum_required(VERSION 2.8...3.12 FATAL_ERROR)


project(freesteam)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../install" CACHE PATH "..." FORCE)
endif()

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/src/*.c" "${PROJECT_SOURCE_DIR}/include/*.h")

# 1. GSL package
set(suffix "${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR}")
if(NOT ${CMAKE_OSX_ARCHITECTURES})
    message(STATUS "CMAKE_OSX_ARCHITECTURES: ${CMAKE_OSX_ARCHITECTURES}")
    set(suffix "${CMAKE_SYSTEM_NAME}_${CMAKE_OSX_ARCHITECTURES}")
endif()
set(GSL_DIR "${CMAKE_SOURCE_DIR}/../install/gsl_${suffix}" CACHE FILEPATH "GSL path contains include, lib")
if(EXISTS ${GSL_DIR})
  set(GSL_INCLUDE_DIR "${GSL_DIR}/include" CACHE FILEPATH "include path of GSL")
  set(GSL_LIBRARY_DIR "${GSL_DIR}/lib" CACHE FILEPATH "lib path of GSL")
  include_directories(${GSL_INCLUDE_DIR})
  if(APPLE)
    set( GSL_LIBRARIES_SHARED ${GSL_LIBRARY_DIR}/libgsl.dylib ${GSL_LIBRARY_DIR}/libgslcblas.dylib )
  else()
    set( GSL_LIBRARIES_SHARED ${GSL_LIBRARY_DIR}/libgsl.so ${GSL_LIBRARY_DIR}/libgslcblas.so )
  endif()
  set( GSL_LIBRARIES_STATIC ${GSL_LIBRARY_DIR}/libgsl.a ${GSL_LIBRARY_DIR}/libgslcblas.a )
  message(STATUS "GSL include path found: " ${GSL_INCLUDE_DIR} )  
  message(STATUS "GSL library path found: " ${GSL_DIR}/lib )     
else()
  message(FATAL_ERROR "GSL_DIR = ${GSL_DIR} doesn't exist.")
endif()
# include_directories(${GSL_INCLUDE_DIR})

add_library (${PROJECT_NAME} STATIC ${SC_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".so")
# specify gsl libraries path for freesteam
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
  install (TARGETS "${PROJECT_NAME}_${buildType}" DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}_${suffix}/lib)
endforeach()
# set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH ${GSL_ROOT_DIR}/lib)
# target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES_STATIC})

# install (TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/lib)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/example/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}_${suffix}/include
          FILES_MATCHING PATTERN "*.h")