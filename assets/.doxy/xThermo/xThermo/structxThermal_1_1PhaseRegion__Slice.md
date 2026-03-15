

# Struct xThermal::PhaseRegion\_Slice



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md)



[More...](#detailed-description)

* `#include <DataStructures.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::map&lt; std::string, std::vector&lt; [**Line\_slice**](structxThermal_1_1Line__slice.md) &gt; &gt; | [**lines**](#variable-lines)  <br> |
|  std::map&lt; std::string, std::vector&lt; [**Point\_slice**](structxThermal_1_1Point__slice.md) &gt; &gt; | [**points**](#variable-points)  <br> |
|  std::map&lt; std::string, std::vector&lt; [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) &gt; &gt; | [**regions**](#variable-regions)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::string | [**info**](#function-info) ([**const**](classxThermal_1_1xThermalError.md) [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) & slice) <br> |




























## Detailed Description


Struct of phase region of a slice at constant T, P, H, or X. Contains phase region polygon, and some visualization properties. 


    
## Public Attributes Documentation




### variable lines 

```C++
std::map<std::string, std::vector<Line_slice> > xThermal::PhaseRegion_Slice::lines;
```




<hr>



### variable points 

```C++
std::map<std::string, std::vector<Point_slice> > xThermal::PhaseRegion_Slice::points;
```




<hr>



### variable regions 

```C++
std::map<std::string, std::vector<Polygon_slice> > xThermal::PhaseRegion_Slice::regions;
```




<hr>
## Public Functions Documentation




### function info 

```C++
inline std::string xThermal::PhaseRegion_Slice::info (
    const  PhaseRegion_Slice & slice
) 
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::PhaseRegion_Slice::operator<< (
    std::ostream & os,
    const  PhaseRegion_Slice & slice
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

