

# Struct xThermal::ThermodynamicProperty



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md)



_Data struct of a thermodynamic property._ 

* `#include <DataStructures.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) | [**dimension**](#variable-dimension)  <br> |
|  [**char**](classxThermal_1_1xThermalError.md) | [**name**](#variable-name)  <br> |
|  [**char**](classxThermal_1_1xThermalError.md) | [**name\_math**](#variable-name_math)  <br> |
|  [**char**](classxThermal_1_1xThermalError.md) | [**unit\_name**](#variable-unit_name)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**value**](#variable-value)  <br> |












































## Public Attributes Documentation




### variable dimension 

```C++
PhysicalDimension xThermal::ThermodynamicProperty::dimension;
```



Dimension of a property in OpenFOAM's format: [kg, m, s, K, mol, A, cd] 


        

<hr>



### variable name 

```C++
char xThermal::ThermodynamicProperty::name[LENGTH_CHAR_NAME];
```



name of a property 


        

<hr>



### variable name\_math 

```C++
char xThermal::ThermodynamicProperty::name_math[LENGTH_CHAR_NAME];
```



Latex format of the name 


        

<hr>



### variable unit\_name 

```C++
char xThermal::ThermodynamicProperty::unit_name[LENGTH_CHAR_NAME];
```



Latex format of the unit 


        

<hr>



### variable value 

```C++
double xThermal::ThermodynamicProperty::value;
```



value of a property 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

