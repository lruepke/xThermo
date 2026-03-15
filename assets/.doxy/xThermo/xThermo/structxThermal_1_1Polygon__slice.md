

# Struct xThermal::Polygon\_slice



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**H**](#variable-h)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**P**](#variable-p)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**T**](#variable-t)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**X**](#variable-x)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**ec**](#variable-ec)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**fc**](#variable-fc)  <br> |
|  std::string | [**name**](#variable-name)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Polygon\_slice**](#function-polygon_slice-12) () <br> |
|   | [**Polygon\_slice**](#function-polygon_slice-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; H0, std::string name0, [**COLOR**](structxThermal_1_1COLOR.md) fc0, [**COLOR**](structxThermal_1_1COLOR.md) ec0) <br> |
|  std::string | [**info**](#function-info) ([**const**](classxThermal_1_1xThermalError.md) [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) & self) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator&lt;**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) & other) const<br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator==**](#function-operator_1) ([**const**](classxThermal_1_1xThermalError.md) [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) & other) const<br> |




























## Public Attributes Documentation




### variable H 

```C++
std::vector<double> xThermal::Polygon_slice::H;
```




<hr>



### variable P 

```C++
std::vector<double> xThermal::Polygon_slice::P;
```




<hr>



### variable T 

```C++
std::vector<double> xThermal::Polygon_slice::T;
```




<hr>



### variable X 

```C++
std::vector<double> xThermal::Polygon_slice::X;
```




<hr>



### variable ec 

```C++
std::vector<double> xThermal::Polygon_slice::ec;
```




<hr>



### variable fc 

```C++
std::vector<double> xThermal::Polygon_slice::fc;
```




<hr>



### variable name 

```C++
std::string xThermal::Polygon_slice::name;
```




<hr>
## Public Functions Documentation




### function Polygon\_slice [1/2]

```C++
inline xThermal::Polygon_slice::Polygon_slice () 
```




<hr>



### function Polygon\_slice [2/2]

```C++
inline xThermal::Polygon_slice::Polygon_slice (
    std::vector< double > T0,
    std::vector< double > P0,
    std::vector< double > X0,
    std::vector< double > H0,
    std::string name0,
    COLOR fc0,
    COLOR ec0
) 
```




<hr>



### function info 

```C++
inline std::string xThermal::Polygon_slice::info (
    const  Polygon_slice & self
) 
```




<hr>



### function operator&lt; 

```C++
inline bool xThermal::Polygon_slice::operator< (
    const  Polygon_slice & other
) const
```




<hr>



### function operator== 

```C++
inline bool xThermal::Polygon_slice::operator== (
    const  Polygon_slice & other
) const
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::Polygon_slice::operator<< (
    std::ostream & os,
    const  Polygon_slice & self
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

