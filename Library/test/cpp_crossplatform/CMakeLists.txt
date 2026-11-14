cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(AppleDeveloper)

PROJECT(example)
set(APP_NAME example_ios)
# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install" CACHE PATH "..." FORCE)
endif()


# OPENFOAM 
add_definitions( -DWM_LABEL_SIZE=32 -DWM_DP )

# 1. generate framework
# add_subdirectory(cpp_framework)
# 2. add ios app
add_subdirectory(app_ios)


# ===================================================================================================
message(STATUS "")
message(STATUS "${PROJECT_NAME} has been configured for ${CMAKE_SYSTEM_NAME}")
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
# ===================================================================================================