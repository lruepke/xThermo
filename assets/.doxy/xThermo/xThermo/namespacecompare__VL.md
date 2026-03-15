

# Namespace compare\_VL



[**Namespace List**](namespaces.md) **>** [**compare\_VL**](namespacecompare__VL.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|   | [**category**](#variable-category)  <br> |
|  float | [**compare**](#variable-compare)   = `lambda a,b (str('%.6e'%(a)))-float(str('%.6e'%(b)))`<br> |
|  str | [**figpath**](#variable-figpath)   = `'.'`<br> |
|  list | [**fmt\_figs**](#variable-fmt_figs)   = `['pdf']`<br> |
|  str | [**phaseBoundary**](#variable-phaseboundary)   = `'VL'`<br> |
|  str | [**result\_path**](#variable-result_path)   = `'VL'`<br> |
|   | [**sw\_84**](#variable-sw_84)   = `H2ONaCl.cH2ONaCl("IAPS84")`<br> |
|   | [**sw\_95**](#variable-sw_95)   = `H2ONaCl.cH2ONaCl("IAPWS95")`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**benchmark\_VL**](#function-benchmark_vl) (sw sw, mmc4 mmc4='../Driesner2007b/1-s2.0-S0016703707002955-mmc4.txt') <br> |
|   | [**calculate**](#function-calculate) (sw sw=sw\_84) <br> |
|   | [**contourf\_phase**](#function-contourf_phase) (ax ax, TT TT, pp pp, phase phase, phase\_name phase\_name, ax\_cb ax\_cb=None) <br> |
|   | [**plot\_3d\_matlab**](#function-plot_3d_matlab) () <br> |
|   | [**plot\_diff**](#function-plot_diff) () <br> |
|   | [**plot\_err**](#function-plot_err) (ax ax, x x, y y, data data, label label='', cmap cmap='rainbow', scale\_data scale\_data='log', vmin vmin=1E-6, vmax vmax=1, s s=1) <br> |
|   | [**plot\_props\_VLH**](#function-plot_props_vlh) (sw sw, water water, mmc5 mmc5='../Driesner2007a/1-s2.0-S0016703707002943-mmc5.txt', mmc3 mmc3='../Driesner2007b/1-s2.0-S0016703707002955-mmc3.txt') <br> |
|   | [**savefig**](#function-savefig) (figname figname) <br> |




























## Public Attributes Documentation




### variable category 

```Python
compare_VL.category;
```




<hr>



### variable compare 

```Python
float compare_VL.compare;
```




<hr>



### variable figpath 

```Python
str compare_VL.figpath;
```




<hr>



### variable fmt\_figs 

```Python
list compare_VL.fmt_figs;
```




<hr>



### variable phaseBoundary 

```Python
str compare_VL.phaseBoundary;
```




<hr>



### variable result\_path 

```Python
str compare_VL.result_path;
```




<hr>



### variable sw\_84 

```Python
compare_VL.sw_84;
```




<hr>



### variable sw\_95 

```Python
compare_VL.sw_95;
```




<hr>
## Public Functions Documentation




### function benchmark\_VL 

```Python
compare_VL::benchmark_VL (
    sw sw,
    mmc4 mmc4='../Driesner2007b/1-s2.0-S0016703707002955-mmc4.txt'
) 
```




<hr>



### function calculate 

```Python
compare_VL::calculate (
    sw sw=sw_84
) 
```




<hr>



### function contourf\_phase 

```Python
compare_VL::contourf_phase (
    ax ax,
    TT TT,
    pp pp,
    phase phase,
    phase_name phase_name,
    ax_cb ax_cb=None
) 
```




<hr>



### function plot\_3d\_matlab 

```Python
compare_VL::plot_3d_matlab () 
```




<hr>



### function plot\_diff 

```Python
compare_VL::plot_diff () 
```




<hr>



### function plot\_err 

```Python
compare_VL::plot_err (
    ax ax,
    x x,
    y y,
    data data,
    label label='',
    cmap cmap='rainbow',
    scale_data scale_data='log',
    vmin vmin=1E-6,
    vmax vmax=1,
    s s=1
) 
```




<hr>



### function plot\_props\_VLH 

```Python
compare_VL::plot_props_VLH (
    sw sw,
    water water,
    mmc5 mmc5='../Driesner2007a/1-s2.0-S0016703707002943-mmc5.txt',
    mmc3 mmc3='../Driesner2007b/1-s2.0-S0016703707002955-mmc3.txt'
) 
```




<hr>



### function savefig 

```Python
compare_VL::savefig (
    figname figname
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/benchmark_test/compare_matlabcode/compare_VL.py`

