

# Class xThermal::OuterEvaluation



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**OuterEvaluation**](classxThermal_1_1OuterEvaluation.md)



_Reset of a caller's result object at the outermost public evaluation entry (C7)._ [More...](#detailed-description)

* `#include <thermo.h>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**OuterEvaluation**](#function-outerevaluation-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|   | [**OuterEvaluation**](#function-outerevaluation-22) ([**const**](classxThermal_1_1xThermalError.md) [**OuterEvaluation**](classxThermal_1_1OuterEvaluation.md) &) = delete<br> |
|  [**OuterEvaluation**](classxThermal_1_1OuterEvaluation.md) & | [**operator=**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**OuterEvaluation**](classxThermal_1_1OuterEvaluation.md) &) = delete<br> |
|   | [**~OuterEvaluation**](#function-outerevaluation) () <br> |




























## Detailed Description


Public functions that fill a [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) object create one of these first. If the call comes from outside the library (evaluation depth 0 on this thread), the object is reset to a default-constructed state, so no field of a previous call survives. Calls the library makes internally (H-P-X -&gt; T-P-X, brine -&gt; water backend) see a depth &gt; 0 and keep their intermediate state, so the result of a fresh evaluation is unchanged. 


    
## Public Functions Documentation




### function OuterEvaluation [1/2]

```C++
explicit xThermal::OuterEvaluation::OuterEvaluation (
    ThermodynamicProperties & props
) 
```




<hr>



### function OuterEvaluation [2/2]

```C++
xThermal::OuterEvaluation::OuterEvaluation (
    const  OuterEvaluation &
) = delete
```




<hr>



### function operator= 

```C++
OuterEvaluation & xThermal::OuterEvaluation::operator= (
    const  OuterEvaluation &
) = delete
```




<hr>



### function ~OuterEvaluation 

```C++
xThermal::OuterEvaluation::~OuterEvaluation () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/thermo.h`

