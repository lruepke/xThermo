cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

PROJECT(example_generic)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../install" CACHE PATH "..." FORCE)
endif()


# ==================================================================
#                generate cpp generic library
# ------------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED example_generic.cpp example_generic.h)

install (TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/generic)
install(FILES example_generic.h DESTINATION ${CMAKE_INSTALL_PREFIX}/generic)
# ==================================================================
