# 0. CMake Minimum version
cmake_minimum_required(VERSION 2.8...3.12 FATAL_ERROR)


set(PROGRAM_NAME waterEOS)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
    project(${PROGRAM_NAME})
else()
    project(${PROGRAM_NAME} LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

include_directories(..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)


#2. Source code dir
aux_source_directory(. SRC_MAIN)
# add Executable target
add_executable(${PROGRAM_NAME}  ${Srcs_bash}  example.cpp)

target_link_libraries(${PROGRAM_NAME} freesteam)


