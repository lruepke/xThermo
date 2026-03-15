

# Struct xThermal::Point\_slice



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**Point\_slice**](structxThermal_1_1Point__slice.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**H**](#variable-h)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P**](#variable-p)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T**](#variable-t)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**X**](#variable-x)  <br> |
|  std::string | [**marker**](#variable-marker)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**mec**](#variable-mec)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**mfc**](#variable-mfc)  <br> |
|  std::string | [**name**](#variable-name)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Point\_slice**](#function-point_slice-12) () <br> |
|   | [**Point\_slice**](#function-point_slice-22) ([**double**](classxThermal_1_1xThermalError.md) T, [**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) X, [**double**](classxThermal_1_1xThermalError.md) H, std::string name, std::string marker, [**COLOR**](structxThermal_1_1COLOR.md) mfc0, [**COLOR**](structxThermal_1_1COLOR.md) mec0) <br> |
|  std::string | [**info**](#function-info) ([**const**](classxThermal_1_1xThermalError.md) [**Point\_slice**](structxThermal_1_1Point__slice.md) & self) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator&lt;**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**Point\_slice**](structxThermal_1_1Point__slice.md) & other) const<br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator==**](#function-operator_1) ([**const**](classxThermal_1_1xThermalError.md) [**Point\_slice**](structxThermal_1_1Point__slice.md) & other) const<br> |




























## Public Attributes Documentation




### variable H 

```C++
double xThermal::Point_slice::H;
```




<hr>



### variable P 

```C++
double xThermal::Point_slice::P;
```




<hr>



### variable T 

```C++
double xThermal::Point_slice::T;
```




<hr>



### variable X 

```C++
double xThermal::Point_slice::X;
```




<hr>



### variable marker 

```C++
std::string xThermal::Point_slice::marker;
```




<hr>



### variable mec 

```C++
std::vector<double> xThermal::Point_slice::mec;
```




<hr>



### variable mfc 

```C++
std::vector<double> xThermal::Point_slice::mfc;
```




<hr>



### variable name 

```C++
std::string xThermal::Point_slice::name;
```




<hr>
## Public Functions Documentation




### function Point\_slice [1/2]

```C++
inline xThermal::Point_slice::Point_slice () 
```




<hr>



### function Point\_slice [2/2]

```C++
inline xThermal::Point_slice::Point_slice (
    double T,
    double P,
    double X,
    double H,
    std::string name,
    std::string marker,
    COLOR mfc0,
    COLOR mec0
) 
```




<hr>



### function info 

```C++
inline std::string xThermal::Point_slice::info (
    const  Point_slice & self
) 
```




<hr>



### function operator&lt; 

```C++
inline bool xThermal::Point_slice::operator< (
    const  Point_slice & other
) const
```




<hr>



### function operator== 

```C++
inline bool xThermal::Point_slice::operator== (
    const  Point_slice & other
) const
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::Point_slice::operator<< (
    std::ostream & os,
    const  Point_slice & self
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

