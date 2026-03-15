

# File test\_prop\_HPX.m



[**FileList**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**Matlab**](dir_a83d64805f95e67d3b37610e6eb407bf.md) **>** [**test\_prop\_HPX.m**](test__prop__HPX_8m.md)

[Go to the source code of this file](test__prop__HPX_8m_source.md)
























## Public Attributes

| Type | Name |
| ---: | :--- |
|  use xThermal | [**API**](#variable-api)   = `/* multi line expression */`<br> |
|   | [**H**](#variable-h)   = `linspace(0.1E6,4.5E6,100)`<br> |
|   | [**Mu**](#variable-mu)   = `reshape(Mu, size(TT))`<br> |
|  end | [**P**](#variable-p)   = `linspace(1,500,250)\*1E5`<br> |
|   | [**PP**](#variable-pp)   = `TT\*0 + P0`<br> |
|   | [**Rho**](#variable-rho)   = `reshape(Rho, size(TT))`<br> |
|   | [**T**](#variable-t)   = `linspace(0.1,900,100)`<br> |
|  P | [**X**](#variable-x)   = `10.^linspace(-10,0,100)`<br> |
|   | [**XX**](#variable-xx)   = `PP\*0 + X0`<br> |
|  P | [**\_\_pad0\_\_**](#variable-__pad0__)  <br> |
|  tic | [**calculate**](#variable-calculate)   = `fluidprop\_NaCl\_PTX( P,T, Xw, visc\_on )`<br> |
|   | [**clf**](#variable-clf)  <br> |
|   | [**colorbar**](#variable-colorbar)  <br> |
|   | [**cols**](#variable-cols)   = `5`<br> |
|  if use shared compiled mex please add libxThermal | [**dylib**](#variable-dylib)  <br> |
|  if use shared compiled mex | [**function**](#variable-function)   = `/* multi line expression */`<br> |
|   | [**nan\_Mu**](#variable-nan_mu)   = `Mu(isnan(Mu))`<br> |
|   | [**nan\_Mu\_PROST**](#variable-nan_mu_prost)   = `props.Mu(isnan(props.Mu))`<br> |
|  hold | [**off**](#variable-off)  <br> |
|  if use shared compiled mex please add libxThermal libgsl libCoolProp dylib to the dynalic library search | [**path**](#variable-path)   = `SinglePointTest(1, 100E5)`<br> |
|  xThermal | [**props**](#variable-props)   = `prop\_TPX(T + 273.15, P, X)`<br> |
|   | [**rho0**](#variable-rho0)   = `1`<br> |
|  toc plot | [**rows**](#variable-rows)   = `3`<br> |
|   | [**xThermal**](#variable-xthermal)   = `prop\_water\_TP(T+273.15, P, name\_EOS)`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  add xThermal path | [**addpath**](#function-addpath) ('/Users/zguo/MyData/Research/3\_CodeProject/Hydrothermal-OpenFOAM/xThermal/Library/API/Matlab') <br> |
|   | [**addpath**](#function-addpath) ('/Users/zguo/MyData/Research/3\_CodeProject/Hydrothermal-OpenFOAM/H2O\_NaCl\_EOS\_Driesner\_2007') <br> |
|   | [**contourf**](#function-contourf) (HH, PP, reshape(PROP.Reg, size(PP))) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, reshape(PROP.Rho, size(PP))) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, reshape(PROP.h, size(PP))) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, reshape(PROP.mu\_l, size(PP))) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, reshape(PROP.mu\_v, size(PP))) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, props. phase) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, (reshape(PROP.Reg, size(PP)) - props.phase)) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, (reshape(PROP.Rho, size(PP)) - props.Rho)) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, (reshape(PROP.h, size(PP)) - props.H)) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, (reshape(PROP.mu\_l, size(PP)) - props.Mu\_l)) <br> |
|   | [**contourf**](#function-contourf) (HH, PP, (reshape(PROP.mu\_v, size(PP)) - props.Mu\_v)) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, reshape(PROP.Reg, size(TT))) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, reshape(PROP.Rho, size(TT))) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, reshape(PROP.h, size(TT))) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, reshape(PROP.mu\_l, size(TT))) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, reshape(PROP.mu\_v, size(TT))) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, props. phase) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, (reshape(PROP.Reg, size(TT)) - props.phase)) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, (reshape(PROP.Rho, size(TT)) - props.Rho)) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, (reshape(PROP.h, size(TT)) - props.H)) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, (reshape(PROP.mu\_l, size(TT)) - props.Mu\_l)) <br> |
|   | [**contourf**](#function-contourf) ((XX), TT, (reshape(PROP.mu\_v, size(TT)) - props.Mu\_v)) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, Rho) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, H) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, log10(Mu)) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, props. Rho) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, log10(props.Mu)) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, (Rho - props.Rho)./ Rho) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, (H - props.H)./ H) <br> |
|   | [**contourf**](#function-contourf) (TT, PP/ 1E5, (Mu - props.Mu)./ Mu) <br> |
|   | [**figure**](#function-figure) (200) <br> |
|   | [**figure**](#function-figure) (100) <br> |
|   | [**fprintf**](#function-fprintf) (1, " Mu\_min=%.3E) <br> |
|   | [**plot**](#function-plot) (T, props. Rho, 'r', 'LineWidth', 6) <br> |
|   | [**plot**](#function-plot) (T, PROP. Rho, 'b-', 'LineWidth', 3) <br> |
|  catch | [**printf**](#function-printf) (0, 'error') <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 1) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 2) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 3) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 4) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 5) <br> |
|  end toc | [**subplot**](#function-subplot) (rows, cols, 6) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 7) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 8) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 9) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 10) <br> |
|  diff | [**subplot**](#function-subplot) (rows, cols, 11) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 12) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 13) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 14) <br> |
|   | [**subplot**](#function-subplot) (rows, cols, 15) <br> |
|   | [**test\_constP\_H2ONaCl**](#function-test_constp_h2onacl) (100E5, 'IAPWS95') <br> |
|  end function | [**test\_constP\_H2ONaCl**](#function-test_constp_h2onacl) (P0, backendName) <br> |
|   | [**test\_constX\_H2ONaCl**](#function-test_constx_h2onacl) (1e- 6, 'IAPS84') <br> |
|  end function | [**test\_constX\_H2ONaCl**](#function-test_constx_h2onacl) (X0, backendName) <br> |
|   | [**title**](#function-title) ('phase region:matlab') <br> |
|   | [**title**](#function-title) ('Rho:matlab') <br> |
|   | [**title**](#function-title) ('H:matlab') <br> |
|   | [**title**](#function-title) ('Mu\_l:matlab') <br> |
|   | [**title**](#function-title) ('Mu\_v:matlab') <br> |
|   | [**title**](#function-title) ('phase region:xThermal:PROST') <br> |
|   | [**title**](#function-title) ('Rho:xThermal:PROST') <br> |
|   | [**title**](#function-title) ('H:xThermal:PROST') <br> |
|   | [**title**](#function-title) ('Mu\_l:xThermal:PROST') <br> |
|   | [**title**](#function-title) ('Mu\_v:xThermal:PROST') <br> |
|   | [**title**](#function-title) ('matlab - xThermal') <br> |
|   | [**title**](#function-title) ('Mu:matlab') <br> |
|   | [**title**](#function-title) () <br> |
|   | [**title**](#function-title) ('Rho:matlab - xThermal') <br> |
|   | [**title**](#function-title) ('H:matlab - xThermal') <br> |
|   | [**title**](#function-title) ('Mu:matlab - xThermal') <br> |




























## Public Attributes Documentation




### variable API 

```Objective-C
use xThermal API;
```




<hr>



### variable H 

```Objective-C
H;
```




<hr>



### variable Mu 

```Objective-C
Mu;
```




<hr>



### variable P 

```Objective-C
P;
```




<hr>



### variable PP 

```Objective-C
PP;
```




<hr>



### variable Rho 

```Objective-C
Rho;
```




<hr>



### variable T 

```Objective-C
T;
```




<hr>



### variable X 

```Objective-C
X;
```




<hr>



### variable XX 

```Objective-C
XX;
```




<hr>



### variable \_\_pad0\_\_ 

```Objective-C
P __pad0__;
```




<hr>



### variable calculate 

```Objective-C
tic calculate;
```




<hr>



### variable clf 

```Objective-C
clf;
```




<hr>



### variable colorbar 

```Objective-C
colorbar;
```




<hr>



### variable cols 

```Objective-C
cols;
```




<hr>



### variable dylib 

```Objective-C
if use shared compiled mex please add libxThermal libgsl dylib;
```




<hr>



### variable function 

```Objective-C
end function[Rho, H, Mu, props];
```




<hr>



### variable nan\_Mu 

```Objective-C
nan_Mu;
```




<hr>



### variable nan\_Mu\_PROST 

```Objective-C
nan_Mu_PROST;
```




<hr>



### variable off 

```Objective-C
hold off;
```




<hr>



### variable path 

```Objective-C
if use shared compiled mex please add libxThermal libgsl libCoolProp dylib to the dynalic library search path[PROP, props];
```




<hr>



### variable props 

```Objective-C
props;
```




<hr>



### variable rho0 

```Objective-C
rho0;
```




<hr>



### variable rows 

```Objective-C
toc plot rows;
```




<hr>



### variable xThermal 

```Objective-C
xThermal[props];
```




<hr>
## Public Functions Documentation




### function addpath 

```Objective-C
add xThermal path addpath (
    '/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/xThermal/Library/API/Matlab'
) 
```




<hr>



### function addpath 

```Objective-C
addpath (
    '/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/H2O_NaCl_EOS_Driesner_2007'
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    reshape(PROP.Reg, size(PP))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    reshape(PROP.Rho, size(PP))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    reshape(PROP.h, size(PP))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    reshape(PROP.mu_l, size(PP))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    reshape(PROP.mu_v, size(PP))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    props. phase
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    (reshape(PROP.Reg, size(PP)) - props.phase)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    (reshape(PROP.Rho, size(PP)) - props.Rho)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    (reshape(PROP.h, size(PP)) - props.H)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    (reshape(PROP.mu_l, size(PP)) - props.Mu_l)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    HH,
    PP,
    (reshape(PROP.mu_v, size(PP)) - props.Mu_v)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    reshape(PROP.Reg, size(TT))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    reshape(PROP.Rho, size(TT))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    reshape(PROP.h, size(TT))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    reshape(PROP.mu_l, size(TT))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    reshape(PROP.mu_v, size(TT))
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    props. phase
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    (reshape(PROP.Reg, size(TT)) - props.phase)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    (reshape(PROP.Rho, size(TT)) - props.Rho)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    (reshape(PROP.h, size(TT)) - props.H)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    (reshape(PROP.mu_l, size(TT)) - props.Mu_l)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    (XX),
    TT,
    (reshape(PROP.mu_v, size(TT)) - props.Mu_v)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    Rho
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    H
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    log10(Mu)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    props. Rho
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    log10(props.Mu)
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    (Rho - props.Rho)./ Rho
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    (H - props.H)./ H
) 
```




<hr>



### function contourf 

```Objective-C
contourf (
    TT,
    PP/ 1E5,
    (Mu - props.Mu)./ Mu
) 
```




<hr>



### function figure 

```Objective-C
figure (
    200
) 
```




<hr>



### function figure 

```Objective-C
figure (
    100
) 
```




<hr>



### function fprintf 

```Objective-C
fprintf (
    1,
    " Mu_min=%.3E
) 
```




<hr>



### function plot 

```Objective-C
plot (
    T,
    props. Rho,
    'r',
    'LineWidth',
    6
) 
```




<hr>



### function plot 

```Objective-C
plot (
    T,
    PROP. Rho,
    'b-',
    'LineWidth',
    3
) 
```




<hr>



### function printf 

```Objective-C
catch printf (
    0,
    'error'
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    1
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    2
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    3
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    4
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    5
) 
```




<hr>



### function subplot 

```Objective-C
end toc subplot (
    rows,
    cols,
    6
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    7
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    8
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    9
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    10
) 
```




<hr>



### function subplot 

```Objective-C
diff subplot (
    rows,
    cols,
    11
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    12
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    13
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    14
) 
```




<hr>



### function subplot 

```Objective-C
subplot (
    rows,
    cols,
    15
) 
```




<hr>



### function test\_constP\_H2ONaCl 

```Objective-C
test_constP_H2ONaCl (
    100E5,
    'IAPWS95'
) 
```




<hr>



### function test\_constP\_H2ONaCl 

```Objective-C
end function test_constP_H2ONaCl (
    P0,
    backendName
) 
```




<hr>



### function test\_constX\_H2ONaCl 

```Objective-C
test_constX_H2ONaCl (
    1e- 6,
    'IAPS84'
) 
```




<hr>



### function test\_constX\_H2ONaCl 

```Objective-C
end function test_constX_H2ONaCl (
    X0,
    backendName
) 
```




<hr>



### function title 

```Objective-C
title (
    'phase region:matlab'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Rho:matlab'
) 
```




<hr>



### function title 

```Objective-C
title (
    'H:matlab'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu_l:matlab'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu_v:matlab'
) 
```




<hr>



### function title 

```Objective-C
title (
    'phase region:xThermal:PROST'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Rho:xThermal:PROST'
) 
```




<hr>



### function title 

```Objective-C
title (
    'H:xThermal:PROST'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu_l:xThermal:PROST'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu_v:xThermal:PROST'
) 
```




<hr>



### function title 

```Objective-C
title (
    'matlab - xThermal'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu:matlab'
) 
```




<hr>



### function title 

```Objective-C
title () 
```




<hr>



### function title 

```Objective-C
title (
    'Rho:matlab - xThermal'
) 
```




<hr>



### function title 

```Objective-C
title (
    'H:matlab - xThermal'
) 
```




<hr>



### function title 

```Objective-C
title (
    'Mu:matlab - xThermal'
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/Matlab/test_prop_HPX.m`

