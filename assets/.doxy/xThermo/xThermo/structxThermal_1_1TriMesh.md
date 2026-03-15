

# Struct xThermal::TriMesh



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**TriMesh**](structxThermal_1_1TriMesh.md)



[More...](#detailed-description)

* `#include <DataStructures.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; std::vector&lt; [**int**](classxThermal_1_1xThermalError.md) &gt; &gt; | [**connection**](#variable-connection)  <br> |
|  std::vector&lt; std::string &gt; | [**names**](#variable-names)   = `{"x","y","z"}`<br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**x**](#variable-x)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**y**](#variable-y)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**z**](#variable-z)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**void**](classxThermal_1_1xThermalError.md) | [**append**](#function-append) ([**const**](classxThermal_1_1xThermalError.md) [**TriMesh**](structxThermal_1_1TriMesh.md) & trimesh) <br> |




























## Detailed Description


2D triangle mesh structure. Relate to matplotlib ax.tri\_plot() in 2D or ax.plot\_trisurf() in 3D 


    
## Public Attributes Documentation




### variable connection 

```C++
std::vector<std::vector<int> > xThermal::TriMesh::connection;
```




<hr>



### variable names 

```C++
std::vector<std::string> xThermal::TriMesh::names;
```




<hr>



### variable x 

```C++
std::vector<double> xThermal::TriMesh::x;
```




<hr>



### variable y 

```C++
std::vector<double> xThermal::TriMesh::y;
```




<hr>



### variable z 

```C++
std::vector<double> xThermal::TriMesh::z;
```




<hr>
## Public Functions Documentation




### function append 

```C++
inline void xThermal::TriMesh::append (
    const  TriMesh & trimesh
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

