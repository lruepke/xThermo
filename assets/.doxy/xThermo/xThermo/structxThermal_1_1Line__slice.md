

# Struct xThermal::Line\_slice



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**Line\_slice**](structxThermal_1_1Line__slice.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**H**](#variable-h)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**P**](#variable-p)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**T**](#variable-t)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**X**](#variable-x)  <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**color**](#variable-color)  <br> |
|  std::string | [**ls**](#variable-ls)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**lw**](#variable-lw)  <br> |
|  std::string | [**name**](#variable-name)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Line\_slice**](#function-line_slice-12) () <br> |
|   | [**Line\_slice**](#function-line_slice-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X0, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; H0, std::string name, [**COLOR**](structxThermal_1_1COLOR.md) c, std::string ls="solid", [**double**](classxThermal_1_1xThermalError.md) lw=1) <br> |
|  std::string | [**info**](#function-info) ([**const**](classxThermal_1_1xThermalError.md) [**Line\_slice**](structxThermal_1_1Line__slice.md) & self) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator&lt;**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**Line\_slice**](structxThermal_1_1Line__slice.md) & other) const<br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator==**](#function-operator_1) ([**const**](classxThermal_1_1xThermalError.md) [**Line\_slice**](structxThermal_1_1Line__slice.md) & other) const<br> |




























## Public Attributes Documentation




### variable H 

```C++
std::vector<double> xThermal::Line_slice::H;
```




<hr>



### variable P 

```C++
std::vector<double> xThermal::Line_slice::P;
```




<hr>



### variable T 

```C++
std::vector<double> xThermal::Line_slice::T;
```




<hr>



### variable X 

```C++
std::vector<double> xThermal::Line_slice::X;
```




<hr>



### variable color 

```C++
std::vector<double> xThermal::Line_slice::color;
```




<hr>



### variable ls 

```C++
std::string xThermal::Line_slice::ls;
```




<hr>



### variable lw 

```C++
double xThermal::Line_slice::lw;
```




<hr>



### variable name 

```C++
std::string xThermal::Line_slice::name;
```




<hr>
## Public Functions Documentation




### function Line\_slice [1/2]

```C++
inline xThermal::Line_slice::Line_slice () 
```




<hr>



### function Line\_slice [2/2]

```C++
inline xThermal::Line_slice::Line_slice (
    std::vector< double > T0,
    std::vector< double > P0,
    std::vector< double > X0,
    std::vector< double > H0,
    std::string name,
    COLOR c,
    std::string ls="solid",
    double lw=1
) 
```




<hr>



### function info 

```C++
inline std::string xThermal::Line_slice::info (
    const  Line_slice & self
) 
```




<hr>



### function operator&lt; 

```C++
inline bool xThermal::Line_slice::operator< (
    const  Line_slice & other
) const
```




<hr>



### function operator== 

```C++
inline bool xThermal::Line_slice::operator== (
    const  Line_slice & other
) const
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::Line_slice::operator<< (
    std::ostream & os,
    const  Line_slice & self
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

