# -------非常重要的设置，必须设置最小版本3.15，否则会出现非常奇怪的问题!!!: ImportError: dlopen(******.so, 0x0002): symbol not found in flat namespace '_PyCMethod_New'
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)


project(test_swig LANGUAGES CXX)

# add_executable(${PROJECT_NAME} main.cpp test.cpp)
add_library(${PROJECT_NAME}_main test.cpp test.H)
set_target_properties(${PROJECT_NAME}_main PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

set(PYTHON_MODULE_NAME ${PROJECT_NAME})
option(Build_API_Python "whether build API for Python" ON)
# -------------- SWIG ---------------------------
set(CMAKE_SWIG_FLAGS)
find_package(SWIG REQUIRED)
include(UseSWIG)
if(${SWIG_VERSION} VERSION_GREATER_EQUAL 4)
  list(APPEND CMAKE_SWIG_FLAGS "-doxygen") #翻译代码注释到最终的.py文件，显示为python注释的格式
endif()
if(UNIX AND NOT APPLE)
  list(APPEND CMAKE_SWIG_FLAGS "-DSWIGWORDSIZE64")
endif()
# -------------- Find Python 3 ---------------------------
# Development.Module这个很关键，但是Cmake版本16没有这个功能，版本22测试通过，其他没试过
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
list(APPEND CMAKE_SWIG_FLAGS "-py3" "-DPY3")
# --------------------------------------------------------
# 设置模块名字，需要保证动态库、.i文件里面的module关键词、最后的.py文件名一致；但是动态库.so文件名会自动加一个下划线的前缀
set(PYTHON_MODULE_NAME test_swig)
# -- 设置源代码为C++
set_property(SOURCE test.i PROPERTY CPLUSPLUS ON)
set_property(SOURCE test.i PROPERTY SWIG_MODULE_NAME ${PYTHON_MODULE_NAME}) #这里的名字是最终的.py文件的名字
# 这里的名字只影响最后的.so文件的名字
swig_add_library(${PYTHON_MODULE_NAME}
  TYPE SHARED
  LANGUAGE python
  OUTPUT_DIR ${CMAKE_BINARY_DIR}
  SOURCES test.i #这里的source可以只写.i文件，然后在链接需要的动态库；也可以在这里把所有需要的源文件写在这里直接编译
  ) 
# 添加一些编译选项使其不显示警告信息。在编译过程中有一些警告信息会提示，后面有一个中括号显示[-Wxxxx]，如果不想显示这类警告则添加编译选项 "-Wno-xxxx" 即可
target_compile_options(${PYTHON_MODULE_NAME} PRIVATE "-Wno-deprecated-declarations")
target_link_libraries(${PYTHON_MODULE_NAME} PRIVATE ${PROJECT_NAME}_main)
# 头文件路径
target_include_directories(${PYTHON_MODULE_NAME}
  PRIVATE
  .     #可以写多个路径
  ${Python3_INCLUDE_DIRS}
  )
# set_property(TARGET test_swig PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON) #可选

if(APPLE)
  set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES
    SUFFIX ".so"
    # INSTALL_RPATH "@loader_path;@loader_path/../../${PYTHON_PROJECT}/.libs"
    )
  # 这一个对于mac系统非常关键
  set_property(TARGET ${PYTHON_MODULE_NAME} APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
elseif(UNIX)
#   set_target_properties(${PYTHON_MODULE_NAME} PROPERTIES
#     INSTALL_RPATH "$ORIGIN:$ORIGIN/../../${PYTHON_PROJECT}/.libs"
#     )
endif()

if(MSVC)
  target_link_libraries(${PYTHON_MODULE_NAME} PRIVATE ${Python3_LIBRARIES})
endif()