

# Class xThermal::IAPWS95::cIAPWS95



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) **>** [**cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md)



_Class of IAPWS-95 formula of H2O, which inherits from_ [_**IAPWS\_IF97::cIAPWS\_IF97**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) _._

* `#include <IAPWS95.h>`



Inherits the following classes: [xThermal::cxThermal](classxThermal_1_1cxThermal.md)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) | [**m\_constants**](#variable-m_constants)  <br> |


## Public Attributes inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](classxThermal_1_1cxThermal.md#variable-m_isshowprogressbar)  <br> |






























## Public Functions

| Type | Name |
| ---: | :--- |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-14) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & T\_K, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-24) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-34) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-14) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_K, [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br>_Solve saturated pressure, liquid density and vapor density on phase boundary of boiling curve by given temperature. Using_ [_**P\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-p_sat_estimate-12) _,_[_**Rho\_Liquid\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-rho_liquid_sat_estimate) _and_[_**Rho\_Vapor\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-rho_vapor_sat_estimate) _calculate initial values of_\(p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime}\) _, respectively; then use GSL solve three nonlinear quations with three unknowns (_\(p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime}\) _) according to the phase-equilibrium condition (See Table 3, pp. 11 in IAPWS-95)._ |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-24) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-34) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Mu**](#function-mu) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate dynamic viscosity for give T,p._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P\_Sat\_estimate**](#function-p_sat_estimate-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_K) <br>_Vapor-pressure equation. See Eq. 2.5 (pp. 399) in IAPWS-95-art._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_Sat\_estimate**](#function-p_sat_estimate-22) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_K, [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & dpdT) <br>_Vapor-pressure and its derivative equation. See Eq. 2.5, 2.5a (pp. 399) in IAPWS-95-art._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P\_delta\_tau**](#function-p_delta_tau) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho**](#function-rho) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_K, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, std::string method="bisection") <br>_Calculate density_ \(\rho\) _by given temperature and pressure._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho\_Liquid\_Sat\_estimate**](#function-rho_liquid_sat_estimate) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Saturated liquid density. See Eq. 2.6 (pp. 399) in IAPWS-95-art._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho\_Vapor\_Sat\_estimate**](#function-rho_vapor_sat_estimate) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Saturated vapor density. See Eq. 2.7 (pp. 399) in IAPWS-95-art._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_Sat\_estimate**](#function-t_sat_estimate) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Solve nonlinear equation of_ [_**P\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-p_sat_estimate-12) _._ |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](#function-updatestate_hpx) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Update state using T, p as independent variables._  |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-22) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Verification\_Mu**](#function-verification_mu) () <br>_Computer-program verification given by mu2008._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**\_enthalpy**](#function-_enthalpy-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_K, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, std::string method="bisection") <br>_Calculate specific enthalpy by given temperature and pressure._  |
|   | [**cIAPWS95**](#function-ciapws95) () <br> |
|  [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_HPX**](#function-findphaseregion_hpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () <br> |
| virtual std::string | [**name**](#function-name) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_o**](#function-phi_o) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau, [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) & phio, [**unsigned**](classxThermal_1_1xThermalError.md) [**int**](classxThermal_1_1xThermalError.md) update\_derivatives=[**Update\_phi\_all**](group__BITMASK__phi.md#define-update_phi_all)) <br>_The ideal-gas part_ \(\phi^o\) _of the dimensionless Helmholtz free energy and its derivatives. See Table 4 in IAPWS-95._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r**](#function-phi_r-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Implementation of residual part_ \(\phi^r\) _of the dimensionless Helmholtz free energy: Eq. 6 in IAPWS-95 and its derivatives, see Table 5 in IAPWS-95._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r**](#function-phi_r-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau, [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) & phir) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r**](#function-phi_r-33) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r\_d**](#function-phi_r_d-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Calculate partial derivative of residual part_ \(\left[ \frac{\partial \phi^r}{\partial \delta} \right]_{\tau}\) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r\_d**](#function-phi_r_d-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r\_dd**](#function-phi_r_dd-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Calculate partial derivative of residual part_ \(\left[ \frac{\partial^2 \phi^r}{\partial \delta^2} \right]_{\tau}\) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r\_dd**](#function-phi_r_dd-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r\_dt**](#function-phi_r_dt-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Calculate partial derivative of residual part_ \(\frac{\partial^2 \phi^r}{\partial \delta \partial \tau}\) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r\_dt**](#function-phi_r_dt-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r\_t**](#function-phi_r_t-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Calculate partial derivative of residual part_ \(\left[ \frac{\partial \phi^r}{\partial \tau} \right]_{\delta}\) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r\_t**](#function-phi_r_t-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**phi\_r\_tt**](#function-phi_r_tt-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & delta, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & tau) <br>_Calculate partial derivative of residual part_ \(\left[ \frac{\partial^2 \phi^r}{\partial \tau^2} \right]_{\delta}\) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**phi\_r\_tt**](#function-phi_r_tt-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; delta, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; tau, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () <br> |
|   | [**~cIAPWS95**](#function-ciapws95) () <br> |


## Public Functions inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_T\_props**](classxThermal_1_1cxThermal.md#function-boiling_t_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_p\_props**](classxThermal_1_1cxThermal.md#function-boiling_p_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
|  [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**Rho\_lookup**](classxThermal_1_1cxThermal.md#function-rho_lookup) ([**double**](classxThermal_1_1xThermalError.md) & Rho\_estimate, [**double**](classxThermal_1_1xThermalError.md) & Rho\_min, [**double**](classxThermal_1_1xThermalError.md) & Rho\_max, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Only lookup density for given T,P,X from loaded LUT. This information can be a good guess of Rho for IAPWS95 solution of Rho._  |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](classxThermal_1_1cxThermal.md#function-t_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](classxThermal_1_1cxThermal.md#function-tmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](classxThermal_1_1cxThermal.md#function-tmin) () = 0<br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Triangulation**](classxThermal_1_1cxThermal.md#function-triangulation) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x\_poly, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y\_poly, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) pointInMesh, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) dxdy, [**TriMesh**](structxThermal_1_1TriMesh.md) & trimesh, [**bool**](classxThermal_1_1xThermalError.md) normalize=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](classxThermal_1_1cxThermal.md#function-ttriple) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & H, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Update state using T, p as independent variables._  |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](classxThermal_1_1cxThermal.md#function-createlut_2d-12) ([**double**](classxThermal_1_1xThermalError.md) xy\_min, [**double**](classxThermal_1_1xThermalError.md) xy\_max, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](classxThermal_1_1cxThermal.md#function-createlut_2d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](classxThermal_1_1cxThermal.md#function-createlut_3d-12) ([**double**](classxThermal_1_1xThermalError.md) xyz\_min, [**double**](classxThermal_1_1xThermalError.md) xyz\_max, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](classxThermal_1_1cxThermal.md#function-createlut_3d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) zmin, [**double**](classxThermal_1_1xThermalError.md) zmax, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
|   | [**cxThermal**](classxThermal_1_1cxThermal.md#function-cxthermal) () <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](classxThermal_1_1cxThermal.md#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
|  [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) | [**getLutInfo**](classxThermal_1_1cxThermal.md#function-getlutinfo) (std::string file\_lut, [**bool**](classxThermal_1_1xThermalError.md) printSummary=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**const**](classxThermal_1_1xThermalError.md) std::map&lt; [**int**](classxThermal_1_1xThermalError.md), [**propInfo**](structxThermal_1_1propInfo.md) &gt; & | [**get\_UpdateWhichProps**](classxThermal_1_1cxThermal.md#function-get_updatewhichprops) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut**](classxThermal_1_1cxThermal.md#function-get_dim_lut) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut\_lookup**](classxThermal_1_1cxThermal.md#function-get_dim_lut_lookup) () <br> |
| virtual [**int**](classxThermal_1_1xThermalError.md) | [**get\_num\_threads**](classxThermal_1_1cxThermal.md#function-get_num_threads) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT**](classxThermal_1_1cxThermal.md#function-get_plut) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT\_lookup**](classxThermal_1_1cxThermal.md#function-get_plut_lookup) () <br> |
|  std::vector&lt; T &gt; | [**linspace**](classxThermal_1_1cxThermal.md#function-linspace-12) (T xmin, T xmax, std::size\_t n, [**bool**](classxThermal_1_1xThermalError.md) isLogScale=[**false**](classxThermal_1_1xThermalError.md)) <br>_Make a linearly spaced vector of points._  |
|  std::vector&lt; T &gt; | [**linspace**](classxThermal_1_1cxThermal.md#function-linspace-22) (T xmindxxmax) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**loadLUT**](classxThermal_1_1cxThermal.md#function-loadlut) ([**const**](classxThermal_1_1xThermalError.md) std::string & filename, [**bool**](classxThermal_1_1xThermalError.md) printStatus=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-24) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-34) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-44) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup\_only**](classxThermal_1_1cxThermal.md#function-lookup_only-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup\_only**](classxThermal_1_1cxThermal.md#function-lookup_only-22) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  T | [**max\_vector**](classxThermal_1_1cxThermal.md#function-max_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  T | [**mean\_vector**](classxThermal_1_1cxThermal.md#function-mean_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  std::vector&lt; [**size\_t**](classxThermal_1_1xThermalError.md) &gt; | [**meshgrid**](classxThermal_1_1cxThermal.md#function-meshgrid) (T xmindxxmax, T ymindyymax, std::vector&lt; T &gt; & XX, std::vector&lt; T &gt; & YY) <br> |
|  T | [**min\_vector**](classxThermal_1_1cxThermal.md#function-min_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](classxThermal_1_1cxThermal.md#function-molar_mass) () = 0<br> |
| virtual std::string | [**name**](classxThermal_1_1cxThermal.md#function-name) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**normalizePhaseBoundaries**](classxThermal_1_1cxThermal.md#function-normalizephaseboundaries) ([**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](classxThermal_1_1cxThermal.md#function-p_critical) () = 0<br> |
| virtual std::string | [**phase\_name**](classxThermal_1_1cxThermal.md#function-phase_name) ([**PhaseRegion**](namespacexThermal.md#enum-phaseregion) phase\_index) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](classxThermal_1_1cxThermal.md#function-pmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](classxThermal_1_1cxThermal.md#function-pmin) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](classxThermal_1_1cxThermal.md#function-rhomass_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomolar\_critical**](classxThermal_1_1cxThermal.md#function-rhomolar_critical) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_binary**](classxThermal_1_1cxThermal.md#function-save_lut_to_binary) (std::string filename) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_vtk**](classxThermal_1_1cxThermal.md#function-save_lut_to_vtk) (std::string filename, [**bool**](classxThermal_1_1xThermalError.md) isNormalizeXYZ=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**set\_num\_threads**](classxThermal_1_1cxThermal.md#function-set_num_threads) ([**int**](classxThermal_1_1xThermalError.md) num\_threads) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**showProgressBar**](classxThermal_1_1cxThermal.md#function-showprogressbar) ([**bool**](classxThermal_1_1xThermalError.md) isShow) <br> |
|  T | [**sum\_vector**](classxThermal_1_1cxThermal.md#function-sum_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeLine2VTU**](classxThermal_1_1cxThermal.md#function-writeline2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Y, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Z, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeMeshGrid2VTK**](classxThermal_1_1cxThermal.md#function-writemeshgrid2vtk) ([**const**](classxThermal_1_1xThermalError.md) std::string & vtkFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x, [**const**](classxThermal_1_1xThermalError.md) std::string & xTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y, [**const**](classxThermal_1_1xThermalError.md) std::string & yTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & z, [**const**](classxThermal_1_1xThermalError.md) std::string & zTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & props, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**propInfo**](structxThermal_1_1propInfo.md) &gt; & propsInfo, [**bool**](classxThermal_1_1xThermalError.md) isNormalize=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writePhaseBoundaries2VTU**](classxThermal_1_1cxThermal.md#function-writephaseboundaries2vtu) (std::string outputPath, [**const**](classxThermal_1_1xThermalError.md) [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries, [**double**](classxThermal_1_1xThermalError.md) scale\_T=1, [**double**](classxThermal_1_1xThermalError.md) scale\_p=1, [**double**](classxThermal_1_1xThermalError.md) scale\_X=1) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeTriMesh2Txt**](classxThermal_1_1cxThermal.md#function-writetrimesh2txt) ([**const**](classxThermal_1_1xThermalError.md) [**TriMesh**](structxThermal_1_1TriMesh.md) & mesh, std::string path\_out=".") <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeXXYYZZ2VTU**](classxThermal_1_1cxThermal.md#function-writexxyyzz2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & XX, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & YY, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & ZZ, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual  | [**~cxThermal**](classxThermal_1_1cxThermal.md#function-cxthermal) () <br> |




## Public Static Functions inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  INDEX\_FLUID | [**validateFluid**](classxThermal_1_1cxThermal.md#function-validatefluid) ([**const**](classxThermal_1_1xThermalError.md) std::string & fluidName) <br> |


















































## Public Attributes Documentation




### variable m\_constants 

```C++
CONSTENTS_Thermo xThermal::IAPWS95::cIAPWS95::m_constants;
```




<hr>
## Public Functions Documentation




### function Boiling\_T [1/4]

```C++
void xThermal::IAPWS95::cIAPWS95::Boiling_T (
    const  double & P,
    double & T_K,
    double & rho_l,
    double & rho_v
) 
```




<hr>



### function Boiling\_T [2/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_T (
    const  double & p
) 
```



Calculate boiling temperature [K] of water for a given p [Pa] 

**Parameters:**


* `p` 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-13)


<hr>



### function Boiling\_T [3/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_T (
    const  double & p,
    double & rho_l,
    double & rho_v
) 
```



Calculate both boiling temperature and density of liquid and vapor phase. 

**Parameters:**


* `p` [Pa] 
* `rho_l` [kg/m3] 
* `rho_v` [kg/m3] 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-23)


<hr>



### function Boiling\_T [4/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_T (
    const  double & p,
    ThermodynamicProperties & props
) 
```



Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-33)


<hr>



### function Boiling\_p [1/4]

_Solve saturated pressure, liquid density and vapor density on phase boundary of boiling curve by given temperature. Using_ [_**P\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-p_sat_estimate-12) _,_[_**Rho\_Liquid\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-rho_liquid_sat_estimate) _and_[_**Rho\_Vapor\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-rho_vapor_sat_estimate) _calculate initial values of_\(p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime}\) _, respectively; then use GSL solve three nonlinear quations with three unknowns (_\(p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime}\) _) according to the phase-equilibrium condition (See Table 3, pp. 11 in IAPWS-95)._
```C++
void xThermal::IAPWS95::cIAPWS95::Boiling_p (
    const  double & T_K,
    double & P,
    double & rho_l,
    double & rho_v
) 
```





**Note:**

The pressure unit in the phase-equilibrium condition(Maxwell criterion) is Pa only if [**xThermal::R**](namespacexThermal.md#variable-r) unit is J/kg/K !!! see Table 3 (pp.11) IAPWS-95. I set [**xThermal::R**](namespacexThermal.md#variable-r) in unit of J/kg/K. 




**Parameters:**


* `T_K` 
* `P` 
* `rho_l` 
* `rho_v` 




        

<hr>



### function Boiling\_p [2/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_p (
    const  double & T
) 
```



Calculate boiling pressure [Pa] of water for a given T [K] 

**Parameters:**


* `T` 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-13)


<hr>



### function Boiling\_p [3/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_p (
    const  double & T,
    double & rho_l,
    double & rho_v
) 
```



Calculate both boiling temperature and density of liquid and vapor phase. 

**Parameters:**


* `T` [K] 
* `rho_l` [kg/m3] 
* `rho_v` [kg/m3] 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-23)


<hr>



### function Boiling\_p [4/4]

```C++
virtual double xThermal::IAPWS95::cIAPWS95::Boiling_p (
    const  double & T,
    ThermodynamicProperties & props
) 
```



Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-33)


<hr>



### function Mu 

_Calculate dynamic viscosity for give T,p._ 
```C++
double xThermal::IAPWS95::cIAPWS95::Mu (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

[Pa s] 





        

<hr>



### function P\_Sat\_estimate [1/2]

_Vapor-pressure equation. See Eq. 2.5 (pp. 399) in IAPWS-95-art._ 
```C++
double xThermal::IAPWS95::cIAPWS95::P_Sat_estimate (
    const  double & T_K
) 
```





**Note:**

This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. IAPWS-95-art. 




**Warning:**

The unit of this function depends on unit of H2O::P\_c, I set critical pressure H2O::P\_c to Pa, so this function return pressure in Pa. 




**Parameters:**


* `T_K` [K] 



**Returns:**

double [Pa] 





        

<hr>



### function P\_Sat\_estimate [2/2]

_Vapor-pressure and its derivative equation. See Eq. 2.5, 2.5a (pp. 399) in IAPWS-95-art._ 
```C++
void xThermal::IAPWS95::cIAPWS95::P_Sat_estimate (
    const  double & T_K,
    double & p,
    double & dpdT
) 
```





**Parameters:**


* `T_K` 
* `p` 
* `dpdT` 



**Returns:**

double 





        

<hr>



### function P\_delta\_tau 

```C++
double xThermal::IAPWS95::cIAPWS95::P_delta_tau (
    const  double & delta,
    const  double & tau
) 
```




<hr>



### function Rho 

_Calculate density_ \(\rho\) _by given temperature and pressure._
```C++
double xThermal::IAPWS95::cIAPWS95::Rho (
    const  double & T_K,
    const  double & P,
    std::string method="bisection"
) 
```



![Image](comparOtherCodes/TP/Rho.svg)
 

**Note:**

The python package **pyIAPWS95** and **pyIAPWS97** can be found at [https://github.com/jjgomera/iapws](https://github.com/jjgomera/iapws), which implements both IAPWS95 and IAPWS-IF97 formula. **freesteam**([http://freesteam.sourceforge.net](http://freesteam.sourceforge.net)) and **PROST** ([http://fluidos.etsii.upm.es/faculty/Jaime\_Carpio/Fumatas\_negas](http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas)) implement IAPWS-IF97 and IAPWS84, respectively. 




**Parameters:**


* `T_K` [K] 
* `P` [Pa] 
* `method` 



**Returns:**

double 





        

<hr>



### function Rho\_Liquid\_Sat\_estimate 

_Saturated liquid density. See Eq. 2.6 (pp. 399) in IAPWS-95-art._ 
```C++
double xThermal::IAPWS95::cIAPWS95::Rho_Liquid_Sat_estimate (
    const  double & T
) 
```





**Note:**

This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. IAPWS-95-art. 




**Parameters:**


* `T_K` [K] 



**Returns:**

double [kg/m3] 





        

<hr>



### function Rho\_Vapor\_Sat\_estimate 

_Saturated vapor density. See Eq. 2.7 (pp. 399) in IAPWS-95-art._ 
```C++
double xThermal::IAPWS95::cIAPWS95::Rho_Vapor_Sat_estimate (
    const  double & T
) 
```





**Note:**

This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. IAPWS-95-art. 




**Parameters:**


* `T_K` [K] 



**Returns:**

double [kg/m3] 





        

<hr>



### function T\_Sat\_estimate 

_Solve nonlinear equation of_ [_**P\_Sat\_estimate**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-p_sat_estimate-12) _._
```C++
double xThermal::IAPWS95::cIAPWS95::T_Sat_estimate (
    const  double & P
) 
```





**Parameters:**


* `P` 



**Returns:**

double 





        

<hr>



### function T\_critical 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::T_critical () 
```



Get the triple point temperature in K 


        
Implements [*xThermal::cxThermal::T\_critical*](classxThermal_1_1cxThermal.md#function-t_critical)


<hr>



### function Tmax 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::Tmax () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmax*](classxThermal_1_1cxThermal.md#function-tmax)


<hr>



### function Tmin 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::Tmin () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmin*](classxThermal_1_1cxThermal.md#function-tmin)


<hr>



### function Ttriple 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::Ttriple () 
```



Get the maximum pressure in Pa 


        
Implements [*xThermal::cxThermal::Ttriple*](classxThermal_1_1cxThermal.md#function-ttriple)


<hr>



### function UpdateState\_HPX 

```C++
virtual void xThermal::IAPWS95::cIAPWS95::UpdateState_HPX (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X=0
) 
```



Implements [*xThermal::cxThermal::UpdateState\_HPX*](classxThermal_1_1cxThermal.md#function-updatestate_hpx-14)


<hr>



### function UpdateState\_TPX [1/2]

_Update state using T, p as independent variables._ 
```C++
virtual void xThermal::IAPWS95::cIAPWS95::UpdateState_TPX (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Calculate the critical molar density in \(mol/m^3\) 




**Parameters:**


* `T` [K] 
* `p` [Pa] 




        
Implements [*xThermal::cxThermal::UpdateState\_TPX*](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14)


<hr>



### function UpdateState\_TPX [2/2]

```C++
ThermodynamicProperties xThermal::IAPWS95::cIAPWS95::UpdateState_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function Verification\_Mu 

_Computer-program verification given by mu2008._ 
```C++
void xThermal::IAPWS95::cIAPWS95::Verification_Mu () 
```





**Returns:**







        

<hr>



### function \_enthalpy [1/2]

_Calculate specific enthalpy by given temperature and pressure._ 
```C++
double xThermal::IAPWS95::cIAPWS95::_enthalpy (
    const  double & T_K,
    const  double & P,
    std::string method="bisection"
) 
```



![Image](comparOtherCodes/TP/H.svg)
 

**Note:**

The python package **pyIAPWS95** and **pyIAPWS97** can be found at [https://github.com/jjgomera/iapws](https://github.com/jjgomera/iapws), which implements both IAPWS95 and IAPWS-IF97 formula. **freesteam**([http://freesteam.sourceforge.net](http://freesteam.sourceforge.net)) and **PROST** ([http://fluidos.etsii.upm.es/faculty/Jaime\_Carpio/Fumatas\_negas](http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas)) implement IAPWS-IF97 and IAPWS84, respectively. 




**Parameters:**


* `T_K` 
* `P` 
* `method` 



**Returns:**

double 





        

<hr>



### function cIAPWS95 

```C++
xThermal::IAPWS95::cIAPWS95::cIAPWS95 () 
```




<hr>



### function findPhaseRegion\_HPX 

```C++
PhaseRegion xThermal::IAPWS95::cIAPWS95::findPhaseRegion_HPX (
    const  double & H,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function findPhaseRegion\_TPX 

```C++
virtual PhaseRegion xThermal::IAPWS95::cIAPWS95::findPhaseRegion_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Return the molar mass in kg/mol 


        
Implements [*xThermal::cxThermal::findPhaseRegion\_TPX*](classxThermal_1_1cxThermal.md#function-findphaseregion_tpx)


<hr>



### function molar\_mass 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::molar_mass () 
```



Return the critical mass density in \(kg/m^3\) 


        
Implements [*xThermal::cxThermal::molar\_mass*](classxThermal_1_1cxThermal.md#function-molar_mass)


<hr>



### function name 

```C++
inline virtual std::string xThermal::IAPWS95::cIAPWS95::name () 
```



Name of the model 


        
Implements [*xThermal::cxThermal::name*](classxThermal_1_1cxThermal.md#function-name)


<hr>



### function p\_critical 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::p_critical () 
```



Return the critical temperature in K 


        
Implements [*xThermal::cxThermal::p\_critical*](classxThermal_1_1cxThermal.md#function-p_critical)


<hr>



### function phi\_o 

_The ideal-gas part_ \(\phi^o\) _of the dimensionless Helmholtz free energy and its derivatives. See Table 4 in IAPWS-95._
```C++
void xThermal::IAPWS95::cIAPWS95::phi_o (
    const  double & delta,
    const  double & tau,
    HelmholtzEnergy_dimensionless & phio,
    unsigned  int update_derivatives=Update_phi_all
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Parameters:**


* `phio` \(\phi^o\) 
* `update_derivatives` Combination of the bitmask value: [**Update\_phi\_d**](group__BITMASK__phi.md#define-update_phi_d), [**Update\_phi\_dd**](group__BITMASK__phi.md#define-update_phi_dd), [**Update\_phi\_t**](group__BITMASK__phi.md#define-update_phi_t), [**Update\_phi\_tt**](group__BITMASK__phi.md#define-update_phi_tt), [**Update\_phi\_dt**](group__BITMASK__phi.md#define-update_phi_dt), [**Update\_phi\_all**](group__BITMASK__phi.md#define-update_phi_all) 




        

<hr>



### function phi\_r [1/3]

_Implementation of residual part_ \(\phi^r\) _of the dimensionless Helmholtz free energy: Eq. 6 in IAPWS-95 and its derivatives, see Table 5 in IAPWS-95._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r (
    const  double & delta,
    const  double & tau
) 
```




\[\begin{matrix}
    \phi^r & = & \sum\limits_{i=1}^{7} n_i \delta^{d_i}\tau^{t_i} + \sum\limits_{i=8}^{51} n_i \delta^{d_i} e^{-\delta^{c_i}}  \tau^{t_i} + \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} + \sum\limits_{i=55}^{56} n_i \delta \Delta^{b_i}\psi
    , \left(
    \Delta  =  \theta^2 + B_i [\color{red}{(\delta -1)^2}]^{a_i} ,
    \theta  =  (1-\tau) + A_i [\color{red}{(\delta -1)^2}]^{1/2\beta_i} ,
    \psi  =  e^{-C_i (\delta -1)^2 - D_i (\tau - 1)^2}
    \right)
    \end{matrix}\]



![Image](IAPWS96_phi/phir_der_T300K.svg)
 

**Note:**

The red dashed line is calculated from a python package of IAPWS pyIAPWS. 




**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r [2/3]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r (
    const  double & delta,
    const  double & tau,
    HelmholtzEnergy_dimensionless & phir
) 
```




<hr>



### function phi\_r [3/3]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function phi\_r\_d [1/2]

_Calculate partial derivative of residual part_ \(\left[ \frac{\partial \phi^r}{\partial \delta} \right]_{\tau}\) _._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r_d (
    const  double & delta,
    const  double & tau
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r\_d [2/2]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r_d (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function phi\_r\_dd [1/2]

_Calculate partial derivative of residual part_ \(\left[ \frac{\partial^2 \phi^r}{\partial \delta^2} \right]_{\tau}\) _._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r_dd (
    const  double & delta,
    const  double & tau
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r\_dd [2/2]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r_dd (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function phi\_r\_dt [1/2]

_Calculate partial derivative of residual part_ \(\frac{\partial^2 \phi^r}{\partial \delta \partial \tau}\) _._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r_dt (
    const  double & delta,
    const  double & tau
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r\_dt [2/2]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r_dt (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function phi\_r\_t [1/2]

_Calculate partial derivative of residual part_ \(\left[ \frac{\partial \phi^r}{\partial \tau} \right]_{\delta}\) _._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r_t (
    const  double & delta,
    const  double & tau
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r\_t [2/2]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r_t (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function phi\_r\_tt [1/2]

_Calculate partial derivative of residual part_ \(\left[ \frac{\partial^2 \phi^r}{\partial \tau^2} \right]_{\delta}\) _._
```C++
double xThermal::IAPWS95::cIAPWS95::phi_r_tt (
    const  double & delta,
    const  double & tau
) 
```





**Parameters:**


* `delta` \(\delta = \rho / \rho_{critical}\) is the reduced density. 
* `tau` \(\tau = T_{critical}/T\) is the inverse reduced temperature. 



**Note:**

The unit of temperature is \(K\). The unit of density if \(kg/m^3\) 




**Returns:**

double 





        

<hr>



### function phi\_r\_tt [2/2]

```C++
void xThermal::IAPWS95::cIAPWS95::phi_r_tt (
    const std::vector< double > delta,
    const std::vector< double > tau,
    std::vector< double > & res
) 
```




<hr>



### function pmax 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::pmax () 
```



Get the minimum pressure in Pa 


        
Implements [*xThermal::cxThermal::pmax*](classxThermal_1_1cxThermal.md#function-pmax)


<hr>



### function pmin 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::pmin () 
```



Get the maximum temperature in K 


        
Implements [*xThermal::cxThermal::pmin*](classxThermal_1_1cxThermal.md#function-pmin)


<hr>



### function rhomass\_critical 

```C++
inline virtual double xThermal::IAPWS95::cIAPWS95::rhomass_critical () 
```



Return the critical pressure in Pa 


        
Implements [*xThermal::cxThermal::rhomass\_critical*](classxThermal_1_1cxThermal.md#function-rhomass_critical)


<hr>



### function ~cIAPWS95 

```C++
xThermal::IAPWS95::cIAPWS95::~cIAPWS95 () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.h`

