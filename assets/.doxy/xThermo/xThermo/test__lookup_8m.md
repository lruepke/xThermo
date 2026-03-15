

# File test\_lookup.m



[**FileList**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**Matlab**](dir_a83d64805f95e67d3b37610e6eb407bf.md) **>** [**test\_lookup.m**](test__lookup_8m.md)

[Go to the source code of this file](test__lookup_8m_source.md)
























## Public Attributes

| Type | Name |
| ---: | :--- |
|  single point lookup | [**H**](#variable-h)   = `H0`<br> |
|   | [**H0**](#variable-h0)   = `1E6`<br> |
|   | [**LineColor**](#variable-linecolor)  <br> |
|   | [**P**](#variable-p)   = `P0`<br> |
|  J kg | [**P0**](#variable-p0)   = `100E5`<br> |
|   | [**PP**](#variable-pp)  <br> |
|  vector lookup | [**T**](#variable-t)   = `linspace(1, 500, 100) + 273.15`<br> |
|   | [**T0**](#variable-t0)   = `100 + 273.15`<br> |
|   | [**X**](#variable-x)   = `X0`<br> |
|  Pa | [**X0**](#variable-x0)   = `0.1`<br> |
|   | [**XX**](#variable-xx)   = `HH\*0 + X0`<br> |
|   | [**clc**](#variable-clc)  <br> |
|   | [**clear**](#variable-clear)  <br> |
|   | [**clf**](#variable-clf)  <br> |
|   | [**file\_LUT\_HPX**](#variable-file_lut_hpx)   = `'lut\_constX\_HP\_10.bin'`<br> |
|   | [**file\_LUT\_HP\_water**](#variable-file_lut_hp_water)   = `'lut\_constX\_HP\_11.bin'`<br> |
|   | [**file\_LUT\_TPX**](#variable-file_lut_tpx)   = `'TPX/lut\_TPX\_9.bin'`<br> |
|   | [**function**](#variable-function)   = `/* multi line expression */`<br> |
|  for | [**i**](#variable-i)   = `/* multi line expression */`<br> |
|  shading | [**interp**](#variable-interp)  <br> |
|   | [**k**](#variable-k)  <br> |
|   | [**ncol**](#variable-ncol)   = `2`<br> |
|   | [**nrow**](#variable-nrow)   = `ceil(length(prop\_names)/ncol)`<br> |
|  end hold | [**off**](#variable-off)  <br> |
|  hold | [**on**](#variable-on)  <br> |
|   | [**prop\_names**](#variable-prop_names)   = `fieldnames(props)`<br> |
|   | [**props**](#variable-props)   = `lookup(file\_LUT, fluidName, H0,P0)`<br> |
|  single point lookup | [**props0**](#variable-props0)   = `lookup(file\_LUT, fluidName, T0,P0,X0)`<br> |
|  kg kg | [**show**](#variable-show)   = `true`<br> |
|   | [**water**](#variable-water)   = `test\_lookup\_HP(file\_LUT\_HP\_water)`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**if**](#function-if) (show) <br> |
|   | [**if**](#function-if) (~ strcmp) <br> |
|   | [**pcolor**](#function-pcolor) (HH/ 1E6, PP/ 1E5, props.) <br> |
|   | [**pcolor**](#function-pcolor) (TT-273. 15, PP/ 1E5, props.) <br> |
|   | [**plot**](#function-plot) (H, props. Rho) <br> |
|   | [**plot**](#function-plot) (T, props. Rho) <br> |
|   | [**props**](#function-props) (cell2mat(prop\_names(i))) <br> |
|   | [**title**](#function-title) (prop\_names(i)) <br> |
|   | [**xlabel**](#function-xlabel) ('Bulk specific enthalpy(MJ/kg)') <br> |
|   | [**xlabel**](#function-xlabel) ('Temperature(deg.C)') <br> |
|   | [**ylabel**](#function-ylabel) ('Pressure(bar)') <br> |




























## Public Attributes Documentation




### variable H 

```Objective-C
end matrix lookup H;
```




<hr>



### variable H0 

```Objective-C
H0;
```




<hr>



### variable LineColor 

```Objective-C
LineColor;
```




<hr>



### variable P 

```Objective-C
P;
```




<hr>



### variable P0 

```Objective-C
J kg P0;
```




<hr>



### variable PP 

```Objective-C
PP;
```




<hr>



### variable T 

```Objective-C
end matrix lookup T;
```




<hr>



### variable T0 

```Objective-C
T0;
```




<hr>



### variable X 

```Objective-C
end X;
```




<hr>



### variable X0 

```Objective-C
Pa X0;
```




<hr>



### variable XX 

```Objective-C
XX;
```




<hr>



### variable clc 

```Objective-C
clc;
```




<hr>



### variable clear 

```Objective-C
clear;
```




<hr>



### variable clf 

```Objective-C
clf;
```




<hr>



### variable file\_LUT\_HPX 

```Objective-C
file_LUT_HPX;
```




<hr>



### variable file\_LUT\_HP\_water 

```Objective-C
file_LUT_HP_water;
```




<hr>



### variable file\_LUT\_TPX 

```Objective-C
file_LUT_TPX;
```




<hr>



### variable function 

```Objective-C
end function;
```




<hr>



### variable i 

```Objective-C
for i;
```




<hr>



### variable interp 

```Objective-C
shading interp;
```




<hr>



### variable k 

```Objective-C
k;
```




<hr>



### variable ncol 

```Objective-C
ncol;
```




<hr>



### variable nrow 

```Objective-C
nrow;
```




<hr>



### variable off 

```Objective-C
hold off;
```




<hr>



### variable on 

```Objective-C
hold on;
```




<hr>



### variable prop\_names 

```Objective-C
prop_names;
```




<hr>



### variable props 

```Objective-C
tic props;
```




<hr>



### variable props0 

```Objective-C
single point lookup props0;
```




<hr>



### variable show 

```Objective-C
kg kg show;
```




<hr>



### variable water 

```Objective-C
water[H, P, props];
```




<hr>
## Public Functions Documentation




### function if 

```Objective-C
if (
    show
) 
```




<hr>



### function if 

```Objective-C
if (
    ~ strcmp
) 
```




<hr>



### function pcolor 

```Objective-C
pcolor (
    HH/ 1E6,
    PP/ 1E5,
    props.
) 
```




<hr>



### function pcolor 

```Objective-C
pcolor (
    TT-273. 15,
    PP/ 1E5,
    props.
) 
```




<hr>



### function plot 

```Objective-C
plot (
    H,
    props. Rho
) 
```




<hr>



### function plot 

```Objective-C
plot (
    T,
    props. Rho
) 
```




<hr>



### function props 

```Objective-C
props (
    cell2mat(prop_names(i))
) 
```




<hr>



### function title 

```Objective-C
title (
    prop_names(i)
) 
```




<hr>



### function xlabel 

```Objective-C
xlabel (
    'Bulk specific enthalpy(MJ/kg)'
) 
```




<hr>



### function xlabel 

```Objective-C
xlabel (
    'Temperature(deg.C)'
) 
```




<hr>



### function ylabel 

```Objective-C
ylabel (
    'Pressure(bar)'
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/Matlab/test_lookup.m`

