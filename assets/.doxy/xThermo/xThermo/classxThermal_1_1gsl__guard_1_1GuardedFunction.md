

# Class xThermal::gsl\_guard::GuardedFunction



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**gsl\_guard**](namespacexThermal_1_1gsl__guard.md) **>** [**GuardedFunction**](classxThermal_1_1gsl__guard_1_1GuardedFunction.md)



[More...](#detailed-description)

* `#include <gsl_guard.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GuardedFunction**](#function-guardedfunction) ([**const**](classxThermal_1_1xThermalError.md) [**gsl\_function**](classxThermal_1_1xThermalError.md) & f, [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* where) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**iterate**](#function-iterate) ([**gsl\_root\_fsolver**](classxThermal_1_1xThermalError.md) \* s) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**set**](#function-set) ([**gsl\_root\_fsolver**](classxThermal_1_1xThermalError.md) \* s, [**double**](classxThermal_1_1xThermalError.md) x\_lower, [**double**](classxThermal_1_1xThermalError.md) x\_upper) <br> |




























## Detailed Description


One-dimensional function for a Brent solve. Must outlive every use of the solver it is set on. 


    
## Public Functions Documentation




### function GuardedFunction 

```C++
inline xThermal::gsl_guard::GuardedFunction::GuardedFunction (
    const  gsl_function & f,
    const  char * where
) 
```




<hr>



### function iterate 

```C++
inline int xThermal::gsl_guard::GuardedFunction::iterate (
    gsl_root_fsolver * s
) 
```



One Brent iteration; throws if the function failed during it. Returns GSL's status (GSL\_SUCCESS). 


        

<hr>



### function set 

```C++
inline void xThermal::gsl_guard::GuardedFunction::set (
    gsl_root_fsolver * s,
    double x_lower,
    double x_upper
) 
```



Check the bracket as gsl\_root\_fsolver\_set / brent\_init would, then set the solver. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/gsl_guard.hpp`

