

# Namespace PhaseDiagramSlice



[**Namespace List**](namespaces.md) **>** [**PhaseDiagramSlice**](namespacePhaseDiagramSlice.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|   | [**sw**](#variable-sw)   = `H2ONaCl.cH2ONaCl("IAPS84")`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PhaseDiagramSlice\_constP**](#function-phasediagramslice_constp) (P0 P0=250E5, X0 X0=0.4, Tmin Tmin=200+273.15, Tmax Tmax=H2ONaCl.T\_MAX\_VLH+50, fmt fmt='jpeg', dpi dpi=900, figpath figpath='.', x0\_log\_linear x0\_log\_linear=0.01, lw lw=0.5, points points=30) <br> |
|   | [**makeLogLinearAxes**](#function-makeloglinearaxes) (ax ax, width\_logaxis width\_logaxis=0.4, xlim\_log xlim\_log=(2E-3, 1), xmax xmax=100, xloc\_major xloc\_major=20, xloc\_minor xloc\_minor=4, yloc yloc='left', color\_log color\_log=(127/255, 0, 1), xlabel xlabel='Bulk salinity(wt % NaCl)', bkcolor bkcolor='None', ylim ylim=None, yloc\_major yloc\_major=None, yloc\_minor yloc\_minor=None, ylabel ylabel=None) <br> |
|   | [**plot\_PhaseRegions**](#function-plot_phaseregions) (PhaseRegions PhaseRegions, ax\_linear\_TPX ax\_linear\_TPX, ax\_log\_TPX ax\_log\_TPX, ax\_linear\_HPX ax\_linear\_HPX, ax\_log\_HPX ax\_log\_HPX, lw lw) <br> |




























## Public Attributes Documentation




### variable sw 

```Python
PhaseDiagramSlice.sw;
```




<hr>
## Public Functions Documentation




### function PhaseDiagramSlice\_constP 

```Python
PhaseDiagramSlice::PhaseDiagramSlice_constP (
    P0 P0=250E5,
    X0 X0=0.4,
    Tmin Tmin=200+273.15,
    Tmax Tmax=H2ONaCl.T_MAX_VLH+50,
    fmt fmt='jpeg',
    dpi dpi=900,
    figpath figpath='.',
    x0_log_linear x0_log_linear=0.01,
    lw lw=0.5,
    points points=30
) 
```




<hr>



### function makeLogLinearAxes 

```Python
PhaseDiagramSlice::makeLogLinearAxes (
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



### function plot\_PhaseRegions 

```Python
PhaseDiagramSlice::plot_PhaseRegions (
    PhaseRegions PhaseRegions,
    ax_linear_TPX ax_linear_TPX,
    ax_log_TPX ax_log_TPX,
    ax_linear_HPX ax_linear_HPX,
    ax_log_HPX ax_log_HPX,
    lw lw
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/example_H2ONaCl/PhaseDiagramSlice.py`

