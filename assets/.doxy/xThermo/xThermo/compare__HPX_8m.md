

# File compare\_HPX.m



[**FileList**](files.md) **>** [**benchmark\_test**](dir_c3298c2f310225057389213c5b89b78a.md) **>** [**compare\_matlabcode**](dir_c0af4561eefcccb3b947513fd062ac1d.md) **>** [**compare\_HPX.m**](compare__HPX_8m.md)

[Go to the source code of this file](compare__HPX_8m_source.md)
























## Public Attributes

| Type | Name |
| ---: | :--- |
|   | [**H**](#variable-h)   = `linspace(Hmin, Hmax, 100)`<br> |
|   | [**Hmax**](#variable-hmax)   = `4.5E6`<br> |
|   | [**Hmin**](#variable-hmin)   = `0.1E6`<br> |
|   | [**P**](#variable-p)   = `linspace(Pmin, Pmax ,100)`<br> |
|   | [**Pmax**](#variable-pmax)   = `600E5`<br> |
|   | [**Pmin**](#variable-pmin)   = `1E5`<br> |
|   | [**fpout\_H**](#variable-fpout_h)   = `fopen('HH.txt','w')`<br> |
|   | [**fpout\_Rho**](#variable-fpout_rho)   = `fopen('RHO.txt','w')`<br> |
|   | [**fpout\_T**](#variable-fpout_t)   = `fopen('TT.txt','w')`<br> |
|   | [**fpout\_p**](#variable-fpout_p)   = `fopen('PP.txt','w')`<br> |
|  for | [**i**](#variable-i)   = `/* multi line expression */`<br> |
|  for | [**j**](#variable-j)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**constX**](#function-constx) (0. 2) <br> |
|  function | [**constX**](#function-constx) (x) <br> |
|  end | [**fclose**](#function-fclose) (fpout\_H) <br> |
|   | [**fclose**](#function-fclose) (fpout\_p) <br> |
|   | [**fclose**](#function-fclose) (fpout\_Rho) <br> |
|   | [**fclose**](#function-fclose) (fpout\_T) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_H, '%.8E ', H) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_p, '%.8E ', p) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Rho, '%.8E ', PROP. rho) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_T, '%.8E ', PROP. T) <br> |
|  end | [**fprintf**](#function-fprintf) (fpout\_H, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_p, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Rho, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_T, '\n') <br> |




























## Public Attributes Documentation




### variable H 

```Objective-C
end matrix lookup H;
```




<hr>



### variable Hmax 

```Objective-C
Hmax;
```




<hr>



### variable Hmin 

```Objective-C
Hmin;
```




<hr>



### variable P 

```Objective-C
P;
```




<hr>



### variable Pmax 

```Objective-C
Pmax;
```




<hr>



### variable Pmin 

```Objective-C
Pmin;
```




<hr>



### variable fpout\_H 

```Objective-C
fpout_H;
```




<hr>



### variable fpout\_Rho 

```Objective-C
fpout_Rho;
```




<hr>



### variable fpout\_T 

```Objective-C
fpout_T;
```




<hr>



### variable fpout\_p 

```Objective-C
fpout_p;
```




<hr>



### variable i 

```Objective-C
for i;
```




<hr>



### variable j 

```Objective-C
for j;
```




<hr>
## Public Functions Documentation




### function constX 

```Objective-C
constX (
    0. 2
) 
```




<hr>



### function constX 

```Objective-C
function constX (
    x
) 
```




<hr>



### function fclose 

```Objective-C
end fclose (
    fpout_H
) 
```




<hr>



### function fclose 

```Objective-C
fclose (
    fpout_p
) 
```




<hr>



### function fclose 

```Objective-C
fclose (
    fpout_Rho
) 
```




<hr>



### function fclose 

```Objective-C
fclose (
    fpout_T
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_H,
    '%.8E ',
    H
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_p,
    '%.8E ',
    p
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_Rho,
    '%.8E ',
    PROP. rho
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_T,
    '%.8E ',
    PROP. T
) 
```




<hr>



### function fprintf 

```Objective-C
end fprintf (
    fpout_H,
    '\n'
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_p,
    '\n'
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_Rho,
    '\n'
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_T,
    '\n'
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/benchmark_test/compare_matlabcode/compare_HPX.m`

