# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
set(AUTHOR "Zhikui Guo")
set(AUTHOR_EMAIL "zguo@geomar.de")

project(demo_swig LANGUAGES CXX C)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

#ignore some warnings
if(WIN32)
else()
    add_definitions("-Wno-writable-strings") # This is invalid for MSVC
    add_compile_options(-fPIC) # this option is important of Linux system.
endif()

# ///////////////// configure swig /////////////////////////
set(CMAKE_SWIG_FLAGS)
find_package(SWIG REQUIRED)
include(UseSWIG)
# //////////////////////////////////////////////////////////

# ///////////////// configure python /////////////////////////
# ---------------------------- Find Python 3 ---------------------------------
#set(PYTHON_EXECUTABLE "path/to/specific/python3")
find_package(Python3 REQUIRED COMPONENTS Interpreter Development Development.Module)
list(APPEND CMAKE_SWIG_FLAGS "-py3" "-DPY3")
message(STATUS "Python include: " ${Python3_INCLUDE_DIRS})
message(STATUS "Python library: " ${Python3_LIBRARIES}) #Python3_LIBRARIES works only if "Development" keyword in the find_package
message(STATUS "Python version: " ${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}.${Python3_VERSION_PATCH})

include_directories(${PROJECT_SOURCE_DIR})
set(PYTHON_PROJECT_DIR ${PROJECT_BINARY_DIR}/python/${PROJECT_NAME})
set(Module demo)
message(STATUS "Build Python module ${Module}")
set_property(SOURCE ${Module}.i PROPERTY CPLUSPLUS ON)
set_property(SOURCE ${Module}.i PROPERTY SWIG_MODULE_NAME ${Module}) # The name is the final generated .py file name
# add dynamic library path
set(Target_name "${Module}_Python")
swig_add_library(${Target_name}
    TYPE SHARED
    LANGUAGE python
    OUTPUT_DIR ${PYTHON_PROJECT_DIR}/${Module}
    SOURCES ${Module}.i demo.cpp  #Option1: only set .i file here and then link dependent libraries; Option2: set .i and all dependent source files, which are passed through the Build_Python_API function (cmake function, called in other place)
)
if(WIN32)
else()
    target_compile_options(${Target_name} PRIVATE "-Wno-deprecated-declarations") #cancle warnings
endif()
target_link_libraries(${Target_name} PRIVATE ${Python3_LIBRARIES})
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
set(INSTALL_PATH_API_PYTHON ${CMAKE_BINARY_DIR}/install/${PROJECT_NAME})
install(TARGETS ${Target_name} DESTINATION ${INSTALL_PATH_API_PYTHON})
install(FILES ${PYTHON_PROJECT_DIR}/${Module}/${Module}.py DESTINATION ${INSTALL_PATH_API_PYTHON})
install(FILES ${CMAKE_SOURCE_DIR}/test.py DESTINATION ${INSTALL_PATH_API_PYTHON})
# ////////////////////////////////////////////////////////////

# add c++ test run
add_executable(demo main.cpp demo.cpp demo.h)