# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(templateClass LANGUAGES CXX C)
set(LIB_NAME "${PROJECT_NAME}")
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

add_executable(test main.cpp)