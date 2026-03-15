
# Class Hierarchy

This inheritance list is sorted roughly, but not completely, alphabetically:


* **class** [**ABSTRACT::abstract**](classABSTRACT_1_1abstract.md)     
    * **class** [**IMPLEMENT\_A::implement\_A**](classIMPLEMENT__A_1_1implement__A.md) 
    * **class** [**IMPLEMENT\_AB::implement\_AB**](classIMPLEMENT__AB_1_1implement__AB.md) 
    * **class** [**IMPLEMENT\_B::implement\_B**](classIMPLEMENT__B_1_1implement__B.md) 
* **class** [**LOOKUPTABLE\_FOREST::LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md) _Pass dimension and data type to the class._ 
* **class** [**MultiProgressBar**](classMultiProgressBar.md) 
* **class** [**cExample**](classcExample.md) 
* **class** [**cExample\_generic**](classcExample__generic.md) 
* **class** [**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) _Top abstract class of the thermodynamic model._     
    * **class** [**xThermal::H2ONaCl::cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md) _Class of_ \(H_2O-NaCl\) _EOS._
    * **class** [**xThermal::IAPWS95::cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md) _Class of IAPWS-95 formula of H2O, which inherits from_ [_**IAPWS\_IF97::cIAPWS\_IF97**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) _._
    * **class** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97**](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) 
    * **class** [**xThermal::NaCl::cNaCl**](classxThermal_1_1NaCl_1_1cNaCl.md) _Class of NaCl EOS._ 
    * **class** [**xThermal::PROST::cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md) 
* **struct** [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md) 
* **struct** [**LOOKUPTABLE\_FOREST::LeafQuad**](structLOOKUPTABLE__FOREST_1_1LeafQuad.md) 
* **struct** [**LOOKUPTABLE\_FOREST::NonLeafQuad**](structLOOKUPTABLE__FOREST_1_1NonLeafQuad.md) 
* **struct** [**LOOKUPTABLE\_FOREST::PropsData**](structLOOKUPTABLE__FOREST_1_1PropsData.md) 
* **struct** [**LOOKUPTABLE\_FOREST::Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) 
* **struct** [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md) 
* **struct** [**LOOKUPTABLE\_FOREST::RMSD\_RefineCriterion**](structLOOKUPTABLE__FOREST_1_1RMSD__RefineCriterion.md) _Property refinement criterion, minimum RMSD of a quadran, if the RMSD of a property in a quadran grater than this criterion, it will be refined._ 
* **struct** [**PROPS**](structPROPS.md) 
* **struct** [**badsubseg**](structbadsubseg.md) 
* **struct** [**badtriang**](structbadtriang.md) 
* **struct** [**behavior**](structbehavior.md) 
* **struct** [**event**](structevent.md) 
* **struct** [**flipstacker**](structflipstacker.md) 
* **struct** [**memorypool**](structmemorypool.md) 
* **struct** [**mesh**](structmesh.md) 
* **struct** [**osub**](structosub.md) 
* **struct** [**otri**](structotri.md) 
* **struct** [**params**](structparams.md) _See_ [https://nlopt.readthedocs.io/en/latest/NLopt\_Tutorial/](https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/) __\(f=(2x + 0)^3 - (-x + 1)^3 = 0\) _._
* **struct** [**quadratic\_params**](structquadratic__params.md) 
* **struct** [**rparams**](structrparams.md) 
* **struct** [**splaynode**](structsplaynode.md) 
* **struct** [**triangulateio**](structtriangulateio.md) 
* **struct** [**xThermal::COLOR**](structxThermal_1_1COLOR.md) 
* **struct** [**xThermal::CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) _为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_ 
* **struct** [**xThermal::DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) 
* **struct** [**xThermal::H2ONaCl::Coeffs\_Estimate\_CriticalT**](structxThermal_1_1H2ONaCl_1_1Coeffs__Estimate__CriticalT.md) _Poly fit coefficient for the estimate critical temperature._ 
* **struct** [**xThermal::H2ONaCl::Coeffs\_Viscosity**](structxThermal_1_1H2ONaCl_1_1Coeffs__Viscosity.md) 
* **struct** [**xThermal::H2ONaCl::Params\_Inversion\_PTX**](structxThermal_1_1H2ONaCl_1_1Params__Inversion__PTX.md) _Parameters for pressure and temperature calculation in VLH zone._ 
* **struct** [**xThermal::H2ONaCl::Params\_P2CriticalT**](structxThermal_1_1H2ONaCl_1_1Params__P2CriticalT.md) _Parameters for inversion of critical temperature from pressure._ 
* **struct** [**xThermal::H2ONaCl::Params\_PX2T\_VL**](structxThermal_1_1H2ONaCl_1_1Params__PX2T__VL.md) 
* **struct** [**xThermal::H2ONaCl::Table4\_Driesner2007a\_CriticalCurve**](structxThermal_1_1H2ONaCl_1_1Table4__Driesner2007a__CriticalCurve.md) 
* **struct** [**xThermal::H2ONaCl::Table6\_Pressure\_VLH**](structxThermal_1_1H2ONaCl_1_1Table6__Pressure__VLH.md) 
* **struct** [**xThermal::H2ONaCl::Table7\_XL\_VL**](structxThermal_1_1H2ONaCl_1_1Table7__XL__VL.md) _Parameters for liquid composition,_ \(X_{NaCl}^{VL,liq}\) _, on V+L coexistence surface. V+L liquid branch. See Table 7 of Driesner2007Part1._
* **struct** [**xThermal::H2ONaCl::Table8\_VaporComposition**](structxThermal_1_1H2ONaCl_1_1Table8__VaporComposition.md) 
* **struct** [**xThermal::Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) 
* **struct** [**xThermal::IAPWS95::Coeff\_phi\_o**](structxThermal_1_1IAPWS95_1_1Coeff__phi__o.md) 
* **struct** [**xThermal::IAPWS95::Coeff\_phi\_r**](structxThermal_1_1IAPWS95_1_1Coeff__phi__r.md) _Numerical values of the coefficients and parameters of the residual part of the dimensionless Helmholtz free energy. Eq.6 in IAPWS-95_ [_**cIAPWS95::phi\_r**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-phi_r-13) _._
* **struct** [**xThermal::IAPWS95::Constants\_Viscosity2008\_Water**](structxThermal_1_1IAPWS95_1_1Constants__Viscosity2008__Water.md) _Constants for water viscosity calculation, see mu2008._ 
* **struct** [**xThermal::IAPWS95::HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) _Dimensionless Helmholtz free energy_ \(\phi = f/(RT)\) _and its partial derivatives._
* **struct** [**xThermal::IAPWS95::HelmholtzEnergy\_dimensionless\_SinglePhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__SinglePhase.md) 
* **struct** [**xThermal::IAPWS95::HelmholtzEnergy\_dimensionless\_TwoPhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__TwoPhase.md) 
* **struct** [**xThermal::IAPWS95::Params\_HP2RhoT**](structxThermal_1_1IAPWS95_1_1Params__HP2RhoT.md) 
* **struct** [**xThermal::IAPWS95::Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md) _Data struct for phase equilibrium solving. It is used in_ [_**func\_PhaseEquilibrium**_](namespacexThermal_1_1IAPWS95.md#function-func_phaseequilibrium) _,_cIAPWS95::Boiling\_P _and_cIAPWS95::Boiling\_T _._
* **struct** [**xThermal::IAPWS95::Params\_TP2Rho**](structxThermal_1_1IAPWS95_1_1Params__TP2Rho.md) 
* **struct** [**xThermal::IAPWS95::Params\_T\_Sat\_estimate**](structxThermal_1_1IAPWS95_1_1Params__T__Sat__estimate.md) 
* **struct** [**xThermal::IAPWS95::STRUCT\_delta\_TwoPhase**](structxThermal_1_1IAPWS95_1_1STRUCT__delta__TwoPhase.md) 
* **struct** [**xThermal::IAPWS95::State**](structxThermal_1_1IAPWS95_1_1State.md) 
* **struct** [**xThermal::IAPWS\_IF97::GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) _Dimensionless Gibbs free energy_ \(\gamma = g/(RT)\) _and its partial derivatives._
* **struct** [**xThermal::IAPWS\_IF97::State\_Region1**](structxThermal_1_1IAPWS__IF97_1_1State__Region1.md) 
* **struct** [**xThermal::IAPWS\_IF97::State\_Region2**](structxThermal_1_1IAPWS__IF97_1_1State__Region2.md) 
* **struct** [**xThermal::IAPWS\_IF97::State\_Region5**](structxThermal_1_1IAPWS__IF97_1_1State__Region5.md) 
* **struct** [**xThermal::Line**](structxThermal_1_1Line.md) 
* **struct** [**xThermal::Line\_slice**](structxThermal_1_1Line__slice.md) 
* **struct** [**xThermal::NaCl::Coeff\_H**](structxThermal_1_1NaCl_1_1Coeff__H.md) _Parameters for enthalpies. See table 5 of Driesner2007Part2._ 
* **struct** [**xThermal::NaCl::Coeff\_Rho**](structxThermal_1_1NaCl_1_1Coeff__Rho.md) _Parameters for halite and liquid NaCl volumetric properties. See table 3 of Driesner2007Part2._ 
* **struct** [**xThermal::NaCl::H\_halite\_ref**](structxThermal_1_1NaCl_1_1H__halite__ref.md) 
* **struct** [**xThermal::PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) 
* **struct** [**xThermal::PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) 
* **struct** [**xThermal::PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) _Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is_ \(kg/m^3\) _, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0]. This idea comes from OpenFOAM._
* **struct** [**xThermal::Point**](structxThermal_1_1Point.md) 
* **struct** [**xThermal::Point\_slice**](structxThermal_1_1Point__slice.md) 
* **struct** [**xThermal::Polygon\_slice**](structxThermal_1_1Polygon__slice.md) 
* **struct** [**xThermal::Surface**](structxThermal_1_1Surface.md) 
* **struct** [**xThermal::ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) 
* **struct** [**xThermal::ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) 
* **struct** [**xThermal::ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) 
* **struct** [**xThermal::ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md) _Data struct of a thermodynamic property._ 
* **struct** [**xThermal::TriMesh**](structxThermal_1_1TriMesh.md) 
* **struct** [**xThermal::propInfo**](structxThermal_1_1propInfo.md) _Information of a thermodynamic property._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region1**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region1.md) 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region2a**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2a.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2a, Eq. (22). See Table 20 of IF97._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region2b**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2b.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2b, Eq. (23). See Table 21 of IF97._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region2c**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2c.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2c, Eq. (24). See Table 22 of IF97._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region3a**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region3a.md) _Coefficients and exponents of the backward equation_ \(T_{3a}(p,h)\) _for subregion 3a in its dimensionless form, Eq. (2) in IF97-Region3._
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Backward\_T\_PH\_Region3b**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region3b.md) _Coefficients and exponents of the backward equation_ \(T_{3a}(p,h)\) _for subregion 3b in its dimensionless form, Eq. (3) in IF97-Region3._
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Boundary3ab**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Boundary3ab.md) 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Region1**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region1.md) _Numerical values of the coefficients and exponents of the dimensionless Gibbs free energy for region 1, Eq. (7) of IF97._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Region2**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region2.md) _Numerical values of the coefficients and exponents of the ideal-gas part_ \(\gamma^o\) _o of the dimensionless Gibbs free energy for region 2, Eq. (16)a, and the residual part_\(\gamma^r\) _._
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Coeff\_Region5**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region5.md) _Numerical values of the coefficients and exponents of the residual part  r of the dimensionless Gibbs free energy for region 5. See Table 37, 38 in IF97._ 
* **struct** [**xThermal::IAPWS\_IF97::cIAPWS\_IF97::Table34**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Table34.md) _Table 34 of IF97 for saturation curve calculation._ 
* **class** **UIResponder**    
    * **class** [**AppDelegate**](interfaceAppDelegate.md) 
* **class** **<UIApplicationDelegate>**    
    * **class** [**AppDelegate**](interfaceAppDelegate.md) 
* **class** **NSObject**    
    * **class** [**CppInterface**](interfaceCppInterface.md) 
* **class** **UIViewController**    
    * **class** [**ViewController**](interfaceViewController.md) 
* **class** **std::exception**    
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
    * **class** [**xThermal::xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
        * **class** [**xThermal::xThermalError**](classxThermal_1_1xThermalError.md)     
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 
            * **class** [**xThermal::ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) 

