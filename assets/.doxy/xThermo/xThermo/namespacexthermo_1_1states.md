

# Namespace xthermo::states



[**Namespace List**](namespaces.md) **>** [**xthermo**](namespacexthermo.md) **>** [**states**](namespacexthermo_1_1states.md)



[More...](#detailed-description)
















## Classes

| Type | Name |
| ---: | :--- |
| class | [**EvaluationError**](classxthermo_1_1states_1_1EvaluationError.md) <br> |
| class | [**InvalidInputError**](classxthermo_1_1states_1_1InvalidInputError.md) <br> |
| class | [**NotAvailableError**](classxthermo_1_1states_1_1NotAvailableError.md) <br> |
| class | [**NotConvergedError**](classxthermo_1_1states_1_1NotConvergedError.md) <br> |
| class | [**OutOfRangeError**](classxthermo_1_1states_1_1OutOfRangeError.md) <br> |
| class | [**States**](classxthermo_1_1states_1_1States.md) <br> |
| class | [**Status**](classxthermo_1_1states_1_1Status.md) <br> |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  tuple | [**FIELDS**](#variable-fields)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**hpx**](#function-hpx) (model model, H H, p p, X X=0.0, \* errors="raise", absent absent="zero") <br> |
|   | [**tpx**](#function-tpx) (model model, T T, p p, X X=0.0, \* errors="raise", absent absent="zero") <br> |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  dict | [**\_ERRORS**](#variable-_errors)   = `/* multi line expression */`<br> |
|  int | [**\_NCODES**](#variable-_ncodes)   = `14`<br> |
|  dict | [**\_PHASE\_FIELDS**](#variable-_phase_fields)   = `{f: f[-1] for f in FIELDS if f[-2:] in ("\_l", "\_v", "\_h") and not f.startswith("S\_")}`<br> |
|  dict | [**\_PRESENT**](#variable-_present)   = `{0: "l", 1: "v", 2: "h", 3: "l", 4: "l", 5: "l", 6: "l", 7: "lv", 8: "lh", 9: "vh", 10: "lv", 11: "lvh"}`<br> |
|  dict | [**\_PRESENT\_LUT**](#variable-_present_lut)   = `{slot: np.array([slot in \_PRESENT.get(code, "") for code in range(\_NCODES)]) for slot in "lvh"}`<br> |
|  dict | [**\_SLOT**](#variable-_slot)   = `{"liquid": "l", "vapour": "v", "vapor": "v", "halite": "h"}`<br> |
















## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**\_evaluate**](#function-_evaluate) (model model, name1 name1, batch batch, a a, p p, X X, errors errors, absent absent) <br> |




## Detailed Description


 


    
## Public Attributes Documentation




### variable FIELDS 

```Python
tuple xthermo::states.FIELDS;
```




<hr>
## Public Functions Documentation




### function hpx 

```Python
xthermo::states::hpx (
    model model,
    H H,
    p p,
    X X=0.0,
    * errors="raise",
    absent absent="zero"
) 
```



 


        

<hr>



### function tpx 

```Python
xthermo::states::tpx (
    model model,
    T T,
    p p,
    X X=0.0,
    * errors="raise",
    absent absent="zero"
) 
```



 


        

<hr>
## Protected Attributes Documentation




### variable \_ERRORS 

```Python
dict xthermo::states._ERRORS;
```




<hr>



### variable \_NCODES 

```Python
int xthermo::states._NCODES;
```




<hr>



### variable \_PHASE\_FIELDS 

```Python
dict xthermo::states._PHASE_FIELDS;
```




<hr>



### variable \_PRESENT 

```Python
dict xthermo::states._PRESENT;
```




<hr>



### variable \_PRESENT\_LUT 

```Python
dict xthermo::states._PRESENT_LUT;
```




<hr>



### variable \_SLOT 

```Python
dict xthermo::states._SLOT;
```




<hr>
## Protected Functions Documentation




### function \_evaluate 

```Python
xthermo::states::_evaluate (
    model model,
    name1 name1,
    batch batch,
    a a,
    p p,
    X X,
    errors errors,
    absent absent
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/API/python/xthermo/states.py`

