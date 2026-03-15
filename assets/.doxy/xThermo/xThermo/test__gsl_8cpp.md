

# File test\_gsl.cpp



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_gsl.cpp**](test__gsl_8cpp.md)

[Go to the source code of this file](test__gsl_8cpp_source.md)

_This test example comes from GSL release 2.7 chapter 38.8._ [More...](#detailed-description)

* `#include <stdlib.h>`
* `#include <stdio.h>`
* `#include <gsl/gsl_vector.h>`
* `#include <gsl/gsl_multiroots.h>`
* `#include <gsl/gsl_errno.h>`
* `#include <gsl/gsl_math.h>`
* `#include <gsl/gsl_roots.h>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**params**](structparams.md) <br>_See_ [https://nlopt.readthedocs.io/en/latest/NLopt\_Tutorial/](https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/) __\(f=(2x + 0)^3 - (-x + 1)^3 = 0\) _._ |
| struct | [**quadratic\_params**](structquadratic__params.md) <br> |
| struct | [**rparams**](structrparams.md) <br> |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**func**](#function-func) (const gsl\_vector \* x, void \* params, gsl\_vector \* f) <br> |
|  int | [**func\_df**](#function-func_df) (const gsl\_vector \* x, void \* params, gsl\_matrix \* J) <br> |
|  int | [**func\_fdf**](#function-func_fdf) (const gsl\_vector \* x, void \* params, gsl\_vector \* f, gsl\_matrix \* J) <br> |
|  int | [**main**](#function-main) () <br> |
|  void | [**print\_state**](#function-print_state) (size\_t iter, gsl\_multiroot\_fsolver \* s) <br> |
|  void | [**print\_state**](#function-print_state) (size\_t iter, gsl\_multiroot\_fdfsolver \* s) <br> |
|  void | [**print\_state2**](#function-print_state2) (size\_t iter, gsl\_multiroot\_fsolver \* s) <br> |
|  void | [**print\_state2**](#function-print_state2) (size\_t iter, gsl\_multiroot\_fdfsolver \* s) <br> |
|  double | [**quadratic**](#function-quadratic) (double x, void \* params) <br> |
|  double | [**quadratic\_deriv**](#function-quadratic_deriv) (double x, void \* params) <br> |
|  void | [**quadratic\_fdf**](#function-quadratic_fdf) (double x, void \* params, double \* y, double \* dy) <br> |
|  int | [**rosenbrock\_df**](#function-rosenbrock_df) (const gsl\_vector \* x, void \* params, gsl\_matrix \* J) <br> |
|  int | [**rosenbrock\_f**](#function-rosenbrock_f) (const gsl\_vector \* x, void \* params, gsl\_vector \* f) <br> |
|  int | [**rosenbrock\_fdf**](#function-rosenbrock_fdf) (const gsl\_vector \* x, void \* params, gsl\_vector \* f, gsl\_matrix \* J) <br> |
|  void | [**test\_NLopt\_tutorial**](#function-test_nlopt_tutorial) () <br> |
|  void | [**test\_NLopt\_tutorial\_fdf**](#function-test_nlopt_tutorial_fdf) () <br> |
|  void | [**test\_gsl\_release2\_7\_chap38\_8**](#function-test_gsl_release2_7_chap38_8) () <br> |
|  void | [**test\_gsl\_release2\_7\_chap38\_8\_fdf**](#function-test_gsl_release2_7_chap38_8_fdf) () <br> |
|  int | [**test\_onedimensional\_root\_finding**](#function-test_onedimensional_root_finding) () <br> |




























## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de))  \(\begin{equation}
\begin{matrix}
f_1(x,y) & =  & a(1-x) \\
f_2(x,y) & = & b(y-x^2)
\end{matrix}, (a=1, b=10)
\end{equation}\) 




**Version:**

0.1 




**Date:**

2021-12-27




**Copyright:**

Copyright (c) 2021 





    
## Public Functions Documentation




### function func 

```C++
int func (
    const gsl_vector * x,
    void * params,
    gsl_vector * f
) 
```




<hr>



### function func\_df 

```C++
int func_df (
    const gsl_vector * x,
    void * params,
    gsl_matrix * J
) 
```




<hr>



### function func\_fdf 

```C++
int func_fdf (
    const gsl_vector * x,
    void * params,
    gsl_vector * f,
    gsl_matrix * J
) 
```




<hr>



### function main 

```C++
int main () 
```




<hr>



### function print\_state 

```C++
void print_state (
    size_t iter,
    gsl_multiroot_fsolver * s
) 
```




<hr>



### function print\_state 

```C++
void print_state (
    size_t iter,
    gsl_multiroot_fdfsolver * s
) 
```




<hr>



### function print\_state2 

```C++
void print_state2 (
    size_t iter,
    gsl_multiroot_fsolver * s
) 
```




<hr>



### function print\_state2 

```C++
void print_state2 (
    size_t iter,
    gsl_multiroot_fdfsolver * s
) 
```




<hr>



### function quadratic 

```C++
double quadratic (
    double x,
    void * params
) 
```




<hr>



### function quadratic\_deriv 

```C++
double quadratic_deriv (
    double x,
    void * params
) 
```




<hr>



### function quadratic\_fdf 

```C++
void quadratic_fdf (
    double x,
    void * params,
    double * y,
    double * dy
) 
```




<hr>



### function rosenbrock\_df 

```C++
int rosenbrock_df (
    const gsl_vector * x,
    void * params,
    gsl_matrix * J
) 
```




<hr>



### function rosenbrock\_f 

```C++
int rosenbrock_f (
    const gsl_vector * x,
    void * params,
    gsl_vector * f
) 
```




<hr>



### function rosenbrock\_fdf 

```C++
int rosenbrock_fdf (
    const gsl_vector * x,
    void * params,
    gsl_vector * f,
    gsl_matrix * J
) 
```




<hr>



### function test\_NLopt\_tutorial 

```C++
void test_NLopt_tutorial () 
```




<hr>



### function test\_NLopt\_tutorial\_fdf 

```C++
void test_NLopt_tutorial_fdf () 
```




<hr>



### function test\_gsl\_release2\_7\_chap38\_8 

```C++
void test_gsl_release2_7_chap38_8 () 
```




<hr>



### function test\_gsl\_release2\_7\_chap38\_8\_fdf 

```C++
void test_gsl_release2_7_chap38_8_fdf () 
```




<hr>



### function test\_onedimensional\_root\_finding 

```C++
int test_onedimensional_root_finding () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/test/test_gsl.cpp`

