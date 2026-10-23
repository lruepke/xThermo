

# Class xThermal::gsl\_guard::GuardedMultirootFunction



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**gsl\_guard**](namespacexThermal_1_1gsl__guard.md) **>** [**GuardedMultirootFunction**](classxThermal_1_1gsl__guard_1_1GuardedMultirootFunction.md)



[More...](#detailed-description)

* `#include <gsl_guard.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GuardedMultirootFunction**](#function-guardedmultirootfunction) ([**const**](classxThermal_1_1xThermalError.md) [**gsl\_multiroot\_function**](classxThermal_1_1xThermalError.md) & f) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**iterate**](#function-iterate) ([**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md) \* s) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**set**](#function-set) ([**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md) \* s, [**const**](classxThermal_1_1xThermalError.md) [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* x) <br> |




























## Detailed Description


Multidimensional function for a hybrid solve: contains callback exceptions; values pass through unchanged. 


    
## Public Functions Documentation




### function GuardedMultirootFunction 

```C++
inline explicit xThermal::gsl_guard::GuardedMultirootFunction::GuardedMultirootFunction (
    const  gsl_multiroot_function & f
) 
```




<hr>



### function iterate 

```C++
inline int xThermal::gsl_guard::GuardedMultirootFunction::iterate (
    gsl_multiroot_fsolver * s
) 
```




<hr>



### function set 

```C++
inline int xThermal::gsl_guard::GuardedMultirootFunction::set (
    gsl_multiroot_fsolver * s,
    const  gsl_vector * x
) 
```



Returns GSL's status unchanged (the caller's convergence handling is not changed). 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/gsl_guard.hpp`

