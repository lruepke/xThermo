All the third party libs will be compiled and install in the path `install` in current directory.

## [GSL](https://www.gnu.org/software/gsl/): [requried]
GSL is used for multidimensional root-finding, e.g., find saturation pressure, liquid and vapor density by given specific T for boling curve of H2O.

* First step: fetch gsl submodule. 

```
git submodule update --init --recursive
```

Do this in the xthermo directory. I will checkout all submodules.


* Need `automake` installation.
* install path: ThirdParties/install/gsl

```
cd ThirdParties/gsl
./autogen.sh
git reset --hard
mkdir build
cd build
../configure --prefix=${PWD}/../../install/gsl
make -j 8
make install
```

This installation will results in a folder `gsl` in the `ThirdParties/install` directory. We need to distinguish between different architectures at the moment. You therefore need to create a soft link that includes the target architecture in the name. 

On Mac with M chips, the extension will be `_Darwin_arm64` and the command is something like this:

```
cd ThirdParties/install
ln -s gsl gsl_Darwin_arm64
```

On Linux systems, it will be something `_Linux_x86_64` or alike and the command is:

```
cd ThirdParties/install
ln -s gsl gsl_Linux_x86_64
```

You can get information on your system with this command: `uname -sm`

If xThermo compilation fails (or the compilation of applications like hydrothermalFoam that link against xthermo), please check the error message for missing gsl library, maybe your system name differs to the examples above.


## [CoolProp](http://www.coolprop.org): [optional]

Only used to compare results and numerical performance. We don't use Windows very much and have done only limited testing. There were some issues ralted to the dyanmic library (DLL) generation. A simple workaround is to simply not link against coolprop; deactivate in the CMakeList.txt file in the root directory.

We need to make a few modifications to the coolprop library, which are handled by the patch below. 

The CoolProp library only provides limited access to the functions declared in `CoolPropLib.h`. If you want access to the low-level functions, you will need to add DLL export statements before the required functions.

For example:

```cpp
COOLPROP_VAR double Tmin(void);
```

The `COOLPROP_VAR` is defined in `DataStructures.h` using the patch below:

```cpp
#ifndef COOLPROP_VAR

#ifdef WIN32
#  ifdef COOLPROP_DLL
#    ifdef COOLPROP_DLL_EXPORT
#      define COOLPROP_VAR __declspec(dllexport)
#    else
#      define COOLPROP_VAR __declspec(dllimport)
#    endif
#  else
#    define COOLPROP_VAR
#  endif
#else
#  define COOLPROP_VAR
#endif

#endif
```


Furthermore, the c++11 detection in the [CPstrings.h](CoolProp/include/CPstrings.h) is performed by a macro named `HAS_MOVE_SEMANTICS`, so need to add a definition in the end of [CMakeLists.txt](CoolProp/CMakeLists.txt): 

```
add_definitions("-DHAS_MOVE_SEMANTICS=1")
if(WIN32)
  if(COOLPROP_SHARED_LIBRARY)
    INSTALL(DIRECTORY ${CMAKE_BINARY_DIR}/Release/ DESTINATION ${COOLPROP_INSTALL_PREFIX}/shared )
  else()
    INSTALL(TARGETS ${LIB_NAME} DESTINATION ${COOLPROP_INSTALL_PREFIX}/static )
  endif()
else()
  INSTALL(TARGETS ${LIB_NAME} DESTINATION ${COOLPROP_INSTALL_PREFIX}/lib )
endif(WIN32)

INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${COOLPROP_INSTALL_PREFIX}/include
      FILES_MATCHING PATTERN "*.h"
      )
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/externals/fmtlib/fmt/ DESTINATION ${COOLPROP_INSTALL_PREFIX}/include/fmt
      FILES_MATCHING PATTERN "*.*"
      )
add_compile_options(-fPIC) # this option is important of Linux system.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
```

**Apply patch**: 

Go into `ThirdParties\coolProp`and apply the patch:

```
git apply ../CoolProp_patches/381c8535.patch
```

Now you should be able to complile the coolprop library:

```
git submodule update --init --recursive
mkdir build
cd build
cmake -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp -DCOOLPROP_STATIC_LIBRARY=ON ..
make
make install
rm -rf *
cmake -DCOOLPROP_INSTALL_PREFIX=${PWD}/../../install/CoolProp -DCOOLPROP_SHARED_LIBRARY=ON ..
make
make install
```

Remember to run command `git reset --hard` to reset the repository after compiling the code.

This installation will results in a folder `coolProp` in the `ThirdParties/install` directory. We again need to distinguish between different architectures at the moment. You therefore need to create a soft link that includes the target architecture in the name. 

On Mac with M chips, the extension will be `_Darwin_arm64` and the command is something like this:

```
cd ThirdParties/install
ln -s CoolProp CoolProp_Darwin_arm64
```

On Linux systems, it will be something `_Linux_x86_64` or alike and the command is:

```
cd ThirdParties/install
ln -s CoolProp CoolProp_Linux_x86_64
```


## [freesteam](http://freesteam.sourceforge.net): [optional]
Only used to compare results and numerical performance.

```
cd freesteam-2.1
mkdir build
cd build
cmake ..
make
make install
```

## [PROST](http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas/PROST%20Properties%20of%20Water%20and%20Steam.htm): [optional]
Only used to compare results and numerical performance.

```
cd PROST
mkdir build
cd build
cmake ..
make
make install
```

If you install freesteam and PROST, you will also have to create the respective platform dependent links using variations of the commands above.