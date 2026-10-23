if(NOT Build_API_Matlab)
    return()
endif(NOT Build_API_Matlab)

message(STATUS "Build Matlab API: TRUE")
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM ENG_LIBRARY MAT_LIBRARY MAIN_PROGRAM MEX_COMPILER MCC_COMPILER )
if(Matlab_FOUND)
    matlab_get_release_name_from_version(${Matlab_VERSION_STRING} release)
    message(STATUS "Matlab release: ${release}")
    message(STATUS "Matlab include: ${Matlab_INCLUDE_DIRS}" )
    message(STATUS "Matlab main program: ${Matlab_MAIN_PROGRAM}" )
    message(STATUS "Matlab Matlab_MEX_LIBRARY: ${Matlab_MEX_LIBRARY}" )
    message(STATUS "Matlab Matlab_MX_LIBRARY: ${Matlab_MX_LIBRARY}" )
    message(STATUS "Matlab Matlab_MEX_COMPILER: ${Matlab_MEX_COMPILER}" )
    message(STATUS "Matlab Matlab_LIBRARIES: ${Matlab_LIBRARIES}" )
endif(Matlab_FOUND)

find_library(Matlab_MEX_BLAS
NAMES libmwblas mwblas
NO_DEFAULT_PATH
PATHS ${Matlab_EXTERN_LIBRARY_DIR} ${Matlab_BINARIES_DIR}
PATH_SUFFIXES ${_matlab_libdir_suffix}
REQUIRED
)
message(STATUS "Matlab BLAS library: ${Matlab_MEX_BLAS}")

# In OSX system, the matlab only has x86_64 version, therefore there are some problem when build matlab api for OSX system with arm64 architecture. 
# So we have to specify target architecture to x86_64 only for matlab api building. That means when building matlab API, we need x86_64 version of gsl, CoolProp, ... library, and do not build any application, such as commandline, desktop
if(APPLE)
    if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        message(STATUS "This is a x86_64 ${CMAKE_SYSTEM_NAME}")
    else()
        message(STATUS "This is not a x86_64 ${CMAKE_SYSTEM_NAME}, the CMAKE_OSX_ARCHITECTURES will be set to x86_64")
        set(CMAKE_OSX_ARCHITECTURES x86_64)
        set(TARGET_SUFFIX "${CMAKE_SYSTEM_NAME}_${CMAKE_OSX_ARCHITECTURES}")
        set(BUILD_MATLAB_SEPARATELY ON)
        message(STATUS "Matlab API will be built separately.")
        set(USE_OMP OFF)
        message(STATUS "In this case, the omp is also not available because the system omp lib is ${CMAKE_SYSTEM_PROCESSOR} version")
    endif()
endif()
message(STATUS "TARGET_SUFFIX: ${TARGET_SUFFIX}")

