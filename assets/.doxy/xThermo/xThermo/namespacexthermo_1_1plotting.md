

# Namespace xthermo::plotting



[**Namespace List**](namespaces.md) **>** [**xthermo**](namespacexthermo.md) **>** [**plotting**](namespacexthermo_1_1plotting.md)



[More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  dict | [**UNITS**](#variable-units)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**createWaterAxes**](#function-createwateraxes) (ax ax, fc fc='whitesmoke') <br> |
|   | [**makeLogLinearAxes**](#function-makeloglinearaxes) (ax ax, width\_logaxis width\_logaxis=0.4, xlim\_log xlim\_log=(2E-3, 1), xmax xmax=100, xloc\_major xloc\_major=20, xloc\_minor xloc\_minor=4, yloc yloc='left', color\_log color\_log=(127/255, 0, 1), xlabel xlabel='Bulk salinity(wt % NaCl)', bkcolor bkcolor='None', ylim ylim=None, yloc\_major yloc\_major=None, yloc\_minor yloc\_minor=None, ylabel ylabel=None) <br> |
|   | [**plot\_Pslice**](#function-plot_pslice) (P0 P0, x0\_log\_linear x0\_log\_linear=0.01, sw sw=None) <br> |
|   | [**plot\_Tslice**](#function-plot_tslice) (T0 T0, x0\_log\_linear x0\_log\_linear=0.01, sw sw=None, plog plog=True) <br> |
|   | [**plot\_phases**](#function-plot_phases) (ax ax, states states, x x="X", y y="T", \* legend=True, failed failed="0.75") <br> |
|   | [**plot\_property**](#function-plot_property) (ax ax, states states, field field, x x="X", y y="T", \* cmap="viridis", norm norm=None, vmin vmin=None, vmax vmax=None, contours contours=None, contour\_kw contour\_kw=None, colorbar colorbar=True, failed failed="0.75", label label=None) <br> |
|   | [**plot\_slice**](#function-plot_slice) (ax ax, sl sl, x x="X", y y=None, \* regions=True, lines lines=True, points points=True, legend legend=False) <br> |
|   | [**slice\_kind**](#function-slice_kind) (sl sl) <br> |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  dict | [**\_ASSEMBLAGE\_NAME**](#variable-_assemblage_name)   = `{"l": "L", "v": "V", "h": "H", "lv": "V+L", "lh": "L+H", "vh": "V+H", "lvh": "V+L+H"}`<br> |
|  dict | [**\_FIELD\_UNITS**](#variable-_field_units)   = `/* multi line expression */`<br> |
|  dict | [**\_PHASE\_COLOURS**](#variable-_phase_colours)   = `/* multi line expression */`<br> |
|  str | [**\_SINGLE\_PHASE\_LABEL**](#variable-_single_phase_label)   = `"Single-phase fluid"`<br> |
|  dict | [**\_STATE\_FIELD**](#variable-_state_field)   = `{"T": "T", "P": "p", "X": "X", "H": "H"}`<br> |
















## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**\_assemblage**](#function-_assemblage) (states states, grid\_t grid\_t) <br> |
|   | [**\_axes\_list**](#function-_axes_list) (ax ax) <br> |
|   | [**\_coords**](#function-_coords) (obj obj, k k) <br> |
|   | [**\_edges**](#function-_edges) (v v) <br> |
|   | [**\_grid**](#function-_grid) (states states, x x, y y) <br> |
|   | [**\_in\_fixed\_limits**](#function-_in_fixed_limits) (a a, px px, py py) <br> |
|   | [**\_name**](#function-_name) (spec spec) <br> |
|   | [**\_water\_strip**](#function-_water_strip) (ax\_water ax\_water, sl sl, y y, P0 P0, sw sw) <br> |




## Detailed Description


 


    
## Public Attributes Documentation




### variable UNITS 

```Python
dict xthermo::plotting.UNITS;
```




<hr>
## Public Functions Documentation




### function createWaterAxes 

```Python
xthermo::plotting::createWaterAxes (
    ax ax,
    fc fc='whitesmoke'
) 
```



 


        

<hr>



### function makeLogLinearAxes 

```Python
xthermo::plotting::makeLogLinearAxes (
    ax ax,
    width_logaxis width_logaxis=0.4,
    xlim_log xlim_log=(2E-3, 1),
    xmax xmax=100,
    xloc_major xloc_major=20,
    xloc_minor xloc_minor=4,
    yloc yloc='left',
    color_log color_log=(127/255, 0, 1),
    xlabel xlabel='Bulk salinity(wt % NaCl)',
    bkcolor bkcolor='None',
    ylim ylim=None,
    yloc_major yloc_major=None,
    yloc_minor yloc_minor=None,
    ylabel ylabel=None
) 
```



 


        

<hr>



### function plot\_Pslice 

```Python
xthermo::plotting::plot_Pslice (
    P0 P0,
    x0_log_linear x0_log_linear=0.01,
    sw sw=None
) 
```



 


        

<hr>



### function plot\_Tslice 

```Python
xthermo::plotting::plot_Tslice (
    T0 T0,
    x0_log_linear x0_log_linear=0.01,
    sw sw=None,
    plog plog=True
) 
```



 


        

<hr>



### function plot\_phases 

```Python
xthermo::plotting::plot_phases (
    ax ax,
    states states,
    x x="X",
    y y="T",
    * legend=True,
    failed failed="0.75"
) 
```



 


        

<hr>



### function plot\_property 

```Python
xthermo::plotting::plot_property (
    ax ax,
    states states,
    field field,
    x x="X",
    y y="T",
    * cmap="viridis",
    norm norm=None,
    vmin vmin=None,
    vmax vmax=None,
    contours contours=None,
    contour_kw contour_kw=None,
    colorbar colorbar=True,
    failed failed="0.75",
    label label=None
) 
```



 


        

<hr>



### function plot\_slice 

```Python
xthermo::plotting::plot_slice (
    ax ax,
    sl sl,
    x x="X",
    y y=None,
    * regions=True,
    lines lines=True,
    points points=True,
    legend legend=False
) 
```



 


        

<hr>



### function slice\_kind 

```Python
xthermo::plotting::slice_kind (
    sl sl
) 
```



 


        

<hr>
## Protected Attributes Documentation




### variable \_ASSEMBLAGE\_NAME 

```Python
dict xthermo::plotting._ASSEMBLAGE_NAME;
```




<hr>



### variable \_FIELD\_UNITS 

```Python
dict xthermo::plotting._FIELD_UNITS;
```




<hr>



### variable \_PHASE\_COLOURS 

```Python
dict xthermo::plotting._PHASE_COLOURS;
```




<hr>



### variable \_SINGLE\_PHASE\_LABEL 

```Python
str xthermo::plotting._SINGLE_PHASE_LABEL;
```




<hr>



### variable \_STATE\_FIELD 

```Python
dict xthermo::plotting._STATE_FIELD;
```




<hr>
## Protected Functions Documentation




### function \_assemblage 

```Python
xthermo::plotting::_assemblage (
    states states,
    grid_t grid_t
) 
```




<hr>



### function \_axes\_list 

```Python
xthermo::plotting::_axes_list (
    ax ax
) 
```




<hr>



### function \_coords 

```Python
xthermo::plotting::_coords (
    obj obj,
    k k
) 
```




<hr>



### function \_edges 

```Python
xthermo::plotting::_edges (
    v v
) 
```



 


        

<hr>



### function \_grid 

```Python
xthermo::plotting::_grid (
    states states,
    x x,
    y y
) 
```



 


        

<hr>



### function \_in\_fixed\_limits 

```Python
xthermo::plotting::_in_fixed_limits (
    a a,
    px px,
    py py
) 
```



 


        

<hr>



### function \_name 

```Python
xthermo::plotting::_name (
    spec spec
) 
```




<hr>



### function \_water\_strip 

```Python
xthermo::plotting::_water_strip (
    ax_water ax_water,
    sl sl,
    y y,
    P0 P0,
    sw sw
) 
```



 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/API/python/xthermo/plotting.py`

