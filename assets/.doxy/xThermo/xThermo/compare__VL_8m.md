

# File compare\_VL.m



[**FileList**](files.md) **>** [**benchmark\_test**](dir_c3298c2f310225057389213c5b89b78a.md) **>** [**compare\_matlabcode**](dir_c0af4561eefcccb3b947513fd062ac1d.md) **>** [**compare\_VL.m**](compare__VL_8m.md)

[Go to the source code of this file](compare__VL_8m_source.md)
























## Public Attributes

| Type | Name |
| ---: | :--- |
|   | [**T**](#variable-t)   = `cat(2, T1, T2)`<br> |
|   | [**T1**](#variable-t1)   = `linspace(0.1,Tcrit\_h2o,200)`<br> |
|   | [**T2**](#variable-t2)   = `linspace(Tcrit\_h2o, 1000, 200)`<br> |
|  end | [**fpout\_T**](#variable-fpout_t)   = `fopen('TT\_VL.txt','w')`<br> |
|   | [**fpout\_Xl**](#variable-fpout_xl)   = `fopen('Xl\_VL.txt','w')`<br> |
|   | [**fpout\_Xv**](#variable-fpout_xv)   = `fopen('Xv\_VL.txt','w')`<br> |
|   | [**fpout\_p**](#variable-fpout_p)   = `fopen('PP\_VL.txt','w')`<br> |
|  for | [**i**](#variable-i)   = `/* multi line expression */`<br> |
|  for | [**iT**](#variable-it)   = `/* multi line expression */`<br> |
|  for | [**ip**](#variable-ip)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  fluidprop\_crit\_T return pressure in MPa | [**P\_crit**](#function-p_crit) (i) <br> |
|   | [**P\_vlh**](#function-p_vlh) (i) <br> |
| virtual end | [**P\_vlh**](#function-p_vlh) (P\_vlh&lt;= 0) = 0<br> |
|  function | [**PhaseBoundary\_VL**](#function-phaseboundary_vl) () <br> |
|  end | [**fclose**](#function-fclose) (fpout\_T) <br> |
|   | [**fclose**](#function-fclose) (fpout\_p) <br> |
|   | [**fclose**](#function-fclose) (fpout\_Xl) <br> |
|   | [**fclose**](#function-fclose) (fpout\_Xv) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_T, '%.8E ', T(iT)+273. 15) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_p, '%.8E ', P(ip)) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Xl, '%.8E ', Xw\_l) <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Xv, '%.8E ', Xw\_v) <br> |
|  end | [**fprintf**](#function-fprintf) (fpout\_T, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_p, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Xl, '\n') <br> |
|   | [**fprintf**](#function-fprintf) (fpout\_Xv, '\n') <br> |
|   | [**if**](#function-if) () <br> |




























## Public Attributes Documentation




### variable T 

```Objective-C
end matrix lookup T;
```




<hr>



### variable T1 

```Objective-C
T1;
```




<hr>



### variable T2 

```Objective-C
T2;
```




<hr>



### variable fpout\_T 

```Objective-C
end fpout_T;
```




<hr>



### variable fpout\_Xl 

```Objective-C
fpout_Xl;
```




<hr>



### variable fpout\_Xv 

```Objective-C
fpout_Xv;
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



### variable iT 

```Objective-C
for iT;
```




<hr>



### variable ip 

```Objective-C
for ip;
```




<hr>
## Public Functions Documentation




### function P\_crit 

```Objective-C
fluidprop_crit_T return pressure in MPa P_crit (
    i
) 
```




<hr>



### function P\_vlh 

```Objective-C
P_vlh (
    i
) 
```




<hr>



### function P\_vlh 

```Objective-C
virtual end P_vlh (
    P_vlh<= 0
) = 0
```




<hr>



### function PhaseBoundary\_VL 

```Objective-C
function PhaseBoundary_VL () 
```




<hr>



### function fclose 

```Objective-C
end fclose (
    fpout_T
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
    fpout_Xl
) 
```




<hr>



### function fclose 

```Objective-C
fclose (
    fpout_Xv
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_T,
    '%.8E ',
    T(iT)+273. 15
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_p,
    '%.8E ',
    P(ip)
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_Xl,
    '%.8E ',
    Xw_l
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_Xv,
    '%.8E ',
    Xw_v
) 
```




<hr>



### function fprintf 

```Objective-C
end fprintf (
    fpout_T,
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
    fpout_Xl,
    '\n'
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    fpout_Xv,
    '\n'
) 
```




<hr>



### function if 

```Objective-C
if () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/benchmark_test/compare_matlabcode/compare_VL.m`

