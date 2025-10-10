# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0") #指定将编译到哪个版本的macOS，也可以消除类似"ld: warning: object file (xxx) was built for newer macOS version (13.0) than being linked (12.0)，"这样的警告信息
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(AUTHOR "Zhikui Guo")
set(AUTHOR_EMAIL "zguo@geomar.de")

project(xThermal LANGUAGES CXX C)
set(LIB_NAME "${PROJECT_NAME}")
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
enable_testing()
#ignore some warnings
if(WIN32)
else()
    add_definitions("-Wno-writable-strings") # This is invalid for MSVC
    add_compile_options(-fPIC) # this option is important of Linux system.
endif()



# =====================================================================================
#                            Build Options
# -------------------------------------------------------------------------------------
#  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install" CACHE PATH "..." FORCE)
#  endif()
#  option(USE_OMP "Enable OpenMP to parallel computing" ON)
option(USE_GSL "Enable GSL" ON)
option(USE_COOLPROP "Enable CoolProp" ON)
option(USE_NUMPY "Enable Numpy" OFF)
option(Build_API_MultiLanguage "whether build API for other langulages" OFF)
option(Build_API_Python "whether build API for Python" OFF)
option(Build_API_Matlab "whether build API for Matlab" OFF)
option(Build_API_Julia "whether build API for Julia" OFF) #cxx is exclusively used as Julia warpper extension
option(Build_API_tcl "whether build API for tcl" OFF)
option(Build_IAPWS_Others "whether build other libs of IAPWS code" OFF)
option(BUILD_SHARED "Build shared library on windows system" OFF)
option(USE_OMP "Enable OpenMP to parallel computing" OFF)
set(TARGET_SUFFIX "${CMAKE_SYSTEM_NAME}_${CMAKE_SYSTEM_PROCESSOR}")
set(BUILD_MATLAB_SEPARATELY OFF)

# --------------------------------------------------------------------------------------

# =====================================================================================
#                            Options
# -------------------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install" CACHE PATH "..." FORCE)
endif()
# --------------------------------------------------------------------------------------

# =====================================================================================
#             configurations
# -------------------------------------------------------------------------------------
include(Matlab)
include(gsl)
include(openmp)
include(apple)
include(PROST)
include(CoolProp)
include(xThermal)
# =====================================================================================

# add library
add_subdirectory(Library)

#add applications
# Do not builld anyother application when build Matlab API
if(NOT BUILD_MATLAB_SEPARATELY)
    add_subdirectory(Application/commandline)
    #add_subdirectory(Application/desktop)
endif()

