
## xThermo

Implementation of the various water EOS and the H2O-NaCl EOS of Driesner and Heinrich, 2007. 

Created by Zhikui Guo.


## Requirements for building from source

1. gcc >=11 [required]
2. GSL [required]. Build it in [ThirdParites](./ThirdParties), see [readme](./ThirdParties/readme.md).
3. Freesteam and PROST [Optional]. cmake variable `Build_IAPWS_Others`. Use `cmake -DBuild_IAPWS_Others=OFF` to disable it.
4. Python installation [Optional]. cmake variable `Build_API_Python`. Use `cmake -DBuild_API_Python=OFF` to disable it.
5. TCL installation [Optinoal]. cmake variable `Build_API_tcl`. Use `cmake -DBuild_API_tcl=OFF` to disable it.
6. OpenMP installation [Optional]. cmake variable `USE_OMP`. Use `cmake -DUSE_OMP=OFF` to disable it.


Finally, build the xthermo package:

```
mkdir build
cd build
cmake ..
make
make install
````

Add an environment variable:

```
export xThermo_DIR=<path_to_installation>/xThermo/install
```

## Building the python bindings

We use swig for linking against multiple languages. To compile the python (and all other) binding, you need to have a swig [required] compiler available. Under Mac, simple install it with brew.

To activate the python binding compilation, change the respective lines in CMakeLists.txt

```
option(Build_API_Python "whether build API for Python" ON)
````

Then compile normally following the workflow above. 

If everything went well, you will have a directory matching your active python version in `install/API/python/`. You can pip install it, or you make it available by updating the python path:

For example.

```
export PYTHONPATH="$PWD/install/API/python/cp3.13:$PYTHONPATH"
```

Afterwards, you should be able to import the modules using statements like this:

```
from xThermal import H2O
```

To make the documentation, the easiest way is to create a matching python environment. Assuming you are having conda or miniconda installed:

```
conda env create -f environment.yaml
```

Creates a python 3.13 virtual environment called py313_xtermo.

Aferwards, you can install the required packages using poetry:

```
poetry install
```

Now you can create the docs:

```
cd doc/sphinx
make html
```



## Julia 

Julia installation is under development. 

Julia installation [Optional]. cmake variable `Build_API_Julia`. Use `cmake -DBuild_API_Julia=OFF` to disable it.
   1. Julia package `CxxWrap` installation [Required]. Environment variable `JlCxx_DIR`. Julia install CxxWrap: `add CxxWrap`
```bash
  export JlCxx_DIR=/Users/zguo/.julia/artifacts/0eb6ffcbe0cde10007a2678013611da80643b728/lib/cmake/JlCxx
  # search path of local package of julia
  export JULIA_LOAD_PATH=/Users/zguo/MyData/Research/3-CodeProject/Hydrothermal-OpenFOAM/xthermo/Library/build:$JULIA_LOAD_PATH
```
