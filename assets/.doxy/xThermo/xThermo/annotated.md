
# Class List


Here are the classes, structs, unions and interfaces with brief descriptions:

* **namespace** [**ABSTRACT**](namespaceABSTRACT.md)     
    * **class** [**abstract**](classABSTRACT_1_1abstract.md)     
* **class** [**AppDelegate**](interfaceAppDelegate.md)     
* **class** [**CppInterface**](interfaceCppInterface.md) 
* **namespace** [**IMPLEMENT\_A**](namespaceIMPLEMENT__A.md)     
    * **class** [**implement\_A**](classIMPLEMENT__A_1_1implement__A.md)     
* **namespace** [**IMPLEMENT\_AB**](namespaceIMPLEMENT__AB.md)     
    * **class** [**implement\_AB**](classIMPLEMENT__AB_1_1implement__AB.md)     
* **namespace** [**IMPLEMENT\_B**](namespaceIMPLEMENT__B.md)     
    * **class** [**implement\_B**](classIMPLEMENT__B_1_1implement__B.md)     
* **namespace** [**LOOKUPTABLE\_FOREST**](namespaceLOOKUPTABLE__FOREST.md)     
    * **struct** [**FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)     
    * **struct** [**LeafQuad**](structLOOKUPTABLE__FOREST_1_1LeafQuad.md)     
    * **class** [**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md) _Pass dimension and data type to the class._     
    * **struct** [**NonLeafQuad**](structLOOKUPTABLE__FOREST_1_1NonLeafQuad.md)     
    * **struct** [**PropsData**](structLOOKUPTABLE__FOREST_1_1PropsData.md)     
    * **union** [**Quad\_data**](unionLOOKUPTABLE__FOREST_1_1Quad__data.md)     
    * **struct** [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md)     
    * **struct** [**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)     
    * **struct** [**RMSD\_RefineCriterion**](structLOOKUPTABLE__FOREST_1_1RMSD__RefineCriterion.md) _Property refinement criterion, minimum RMSD of a quadran, if the RMSD of a property in a quadran grater than this criterion, it will be refined._     
* **class** [**MultiProgressBar**](classMultiProgressBar.md)     
* **struct** [**PROPS**](structPROPS.md)     
* **namespace** [**PhaseDiagramSlice**](namespacePhaseDiagramSlice.md)     
* **namespace** [**TEST\_SWIG**](namespaceTEST__SWIG.md)     
* **class** [**ViewController**](interfaceViewController.md) 
* **struct** [**badsubseg**](structbadsubseg.md)     
* **struct** [**badtriang**](structbadtriang.md)     
* **struct** [**behavior**](structbehavior.md)     
* **class** [**cExample**](classcExample.md)     
* **class** [**cExample\_generic**](classcExample__generic.md)     
* **namespace** [**compare\_VL**](namespacecompare__VL.md)     
* **struct** [**event**](structevent.md)     
* **struct** [**flipstacker**](structflipstacker.md)     
* **namespace** [**helpfunc**](namespacehelpfunc.md)     
* **struct** [**memorypool**](structmemorypool.md)     
* **struct** [**mesh**](structmesh.md)     
* **struct** [**osub**](structosub.md)     
* **struct** [**otri**](structotri.md)     
* **struct** [**params**](structparams.md) _See_ [https://nlopt.readthedocs.io/en/latest/NLopt\_Tutorial/](https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/) __\(f=(2x + 0)^3 - (-x + 1)^3 = 0\) _._    
* **struct** [**quadratic\_params**](structquadratic__params.md)     
* **struct** [**rparams**](structrparams.md)     
* **struct** [**splaynode**](structsplaynode.md)     
* **namespace** [**std**](namespacestd.md) 
* **namespace** [**test**](namespacetest.md) 
* **namespace** [**test\_H2O**](namespacetest__H2O.md)     
* **namespace** [**test\_H2ONaCl**](namespacetest__H2ONaCl.md)     
* **struct** [**triangulateio**](structtriangulateio.md)     
* **namespace** [**xThermal**](namespacexThermal.md) _Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._    
    * **struct** [**COLOR**](structxThermal_1_1COLOR.md)     
    * **struct** [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) _为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_     
    * **struct** [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md)     
    * **namespace** [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md)     
        * **struct** [**Coeffs\_Estimate\_CriticalT**](structxThermal_1_1H2ONaCl_1_1Coeffs__Estimate__CriticalT.md) _Poly fit coefficient for the estimate critical temperature._     
        * **struct** [**Coeffs\_Viscosity**](structxThermal_1_1H2ONaCl_1_1Coeffs__Viscosity.md)     
        * **struct** [**Params\_Inversion\_PTX**](structxThermal_1_1H2ONaCl_1_1Params__Inversion__PTX.md) _Parameters for pressure and temperature calculation in VLH zone._     
        * **struct** [**Params\_P2CriticalT**](structxThermal_1_1H2ONaCl_1_1Params__P2CriticalT.md) _Parameters for inversion of critical temperature from pressure._     
        * **struct** [**Params\_PX2T\_VL**](structxThermal_1_1H2ONaCl_1_1Params__PX2T__VL.md)     
        * **struct** [**Table4\_Driesner2007a\_CriticalCurve**](structxThermal_1_1H2ONaCl_1_1Table4__Driesner2007a__CriticalCurve.md)     
        * **struct** [**Table6\_Pressure\_VLH**](structxThermal_1_1H2ONaCl_1_1Table6__Pressure__VLH.md)     
        * **struct** [**Table7\_XL\_VL**](structxThermal_1_1H2ONaCl_1_1Table7__XL__VL.md) _Parameters for liquid composition,_ \(X_{NaCl}^{VL,liq}\) _, on V+L coexistence surface. V+L liquid branch. See Table 7 of Driesner2007Part1._    
        * **struct** [**Table8\_VaporComposition**](structxThermal_1_1H2ONaCl_1_1Table8__VaporComposition.md)     
        * **class** [**cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md) _Class of_ \(H_2O-NaCl\) _EOS._    
    * **struct** [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md)     
    * **namespace** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md)     
        * **struct** [**Coeff\_phi\_o**](structxThermal_1_1IAPWS95_1_1Coeff__phi__o.md)     
        * **struct** [**Coeff\_phi\_r**](structxThermal_1_1IAPWS95_1_1Coeff__phi__r.md) _Numerical values of the coefficients and parameters of the residual part of the dimensionless Helmholtz free energy. Eq.6 in IAPWS-95_ [_**cIAPWS95::phi\_r**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-phi_r-13) _._    
        * **struct** [**Constants\_Viscosity2008\_Water**](structxThermal_1_1IAPWS95_1_1Constants__Viscosity2008__Water.md) _Constants for water viscosity calculation, see mu2008._     
        * **struct** [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) _Dimensionless Helmholtz free energy_ \(\phi = f/(RT)\) _and its partial derivatives._    
        * **struct** [**HelmholtzEnergy\_dimensionless\_SinglePhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__SinglePhase.md)     
        * **struct** [**HelmholtzEnergy\_dimensionless\_TwoPhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__TwoPhase.md)     
        * **struct** [**Params\_HP2RhoT**](structxThermal_1_1IAPWS95_1_1Params__HP2RhoT.md)     
        * **struct** [**Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md) _Data struct for phase equilibrium solving. It is used in_ [_**func\_PhaseEquilibrium**_](namespacexThermal_1_1IAPWS95.md#function-func_phaseequilibrium) _,_cIAPWS95::Boiling\_P _and_cIAPWS95::Boiling\_T _._    
        * **struct** [**Params\_TP2Rho**](structxThermal_1_1IAPWS95_1_1Params__TP2Rho.md)     
        * **struct** [**Params\_T\_Sat\_estimate**](structxThermal_1_1IAPWS95_1_1Params__T__Sat__estimate.md)     
        * **struct** [**STRUCT\_delta\_TwoPhase**](structxThermal_1_1IAPWS95_1_1STRUCT__delta__TwoPhase.md)     
        * **struct** [**State**](structxThermal_1_1IAPWS95_1_1State.md)     
        * **class** [**cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md) _Class of IAPWS-95 formula of H2O, which inherits from_ [_**IAPWS\_IF97::cIAPWS\_IF97**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) _._    
    * **namespace** [**IAPWS\_IF97**](namespacexThermal_1_1IAPWS__IF97.md)     
        * **struct** [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) _Dimensionless Gibbs free energy_ \(\gamma = g/(RT)\) _and its partial derivatives._    
        * **struct** [**State\_Region1**](structxThermal_1_1IAPWS__IF97_1_1State__Region1.md)     
        * **struct** [**State\_Region2**](structxThermal_1_1IAPWS__IF97_1_1State__Region2.md)     
        * **struct** [**State\_Region5**](structxThermal_1_1IAPWS__IF97_1_1State__Region5.md)     
        * **class** [**cIAPWS\_IF97**](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md)     
    * **struct** [**Line**](structxThermal_1_1Line.md)     
    * **struct** [**Line\_slice**](structxThermal_1_1Line__slice.md)     
    * **namespace** [**NaCl**](namespacexThermal_1_1NaCl.md)     
        * **struct** [**Coeff\_H**](structxThermal_1_1NaCl_1_1Coeff__H.md) _Parameters for enthalpies. See table 5 of Driesner2007Part2._     
        * **struct** [**Coeff\_Rho**](structxThermal_1_1NaCl_1_1Coeff__Rho.md) _Parameters for halite and liquid NaCl volumetric properties. See table 3 of Driesner2007Part2._     
        * **struct** [**H\_halite\_ref**](structxThermal_1_1NaCl_1_1H__halite__ref.md)     
        * **class** [**cNaCl**](classxThermal_1_1NaCl_1_1cNaCl.md) _Class of NaCl EOS._     
    * **namespace** [**PROST**](namespacexThermal_1_1PROST.md)     
        * **class** [**cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md)     
    * **struct** [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md)     
    * **struct** [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md)     
    * **struct** [**PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) _Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is_ \(kg/m^3\) _, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0]. This idea comes from OpenFOAM._    
    * **struct** [**Point**](structxThermal_1_1Point.md)     
    * **struct** [**Point\_slice**](structxThermal_1_1Point__slice.md)     
    * **struct** [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md)     
    * **struct** [**Surface**](structxThermal_1_1Surface.md)     
    * **struct** [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md)     
    * **struct** [**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md)     
    * **struct** [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md)     
    * **struct** [**ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md) _Data struct of a thermodynamic property._     
    * **struct** [**TriMesh**](structxThermal_1_1TriMesh.md)     
    * **class** [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md)     
    * **class** [**cxThermal**](classxThermal_1_1cxThermal.md) _Top abstract class of the thermodynamic model._     
    * **struct** [**propInfo**](structxThermal_1_1propInfo.md) _Information of a thermodynamic property._     
    * **class** [**xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)     
    * **class** [**xThermalError**](classxThermal_1_1xThermalError.md)     
* **struct** [**Coeff\_Backward\_T\_PH\_Region1**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region1.md)     
* **struct** [**Coeff\_Backward\_T\_PH\_Region2a**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2a.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2a, Eq. (22). See Table 20 of IF97._     
* **struct** [**Coeff\_Backward\_T\_PH\_Region2b**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2b.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2b, Eq. (23). See Table 21 of IF97._     
* **struct** [**Coeff\_Backward\_T\_PH\_Region2c**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region2c.md) _Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2c, Eq. (24). See Table 22 of IF97._     
* **struct** [**Coeff\_Backward\_T\_PH\_Region3a**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region3a.md) _Coefficients and exponents of the backward equation_ \(T_{3a}(p,h)\) _for subregion 3a in its dimensionless form, Eq. (2) in IF97-Region3._    
* **struct** [**Coeff\_Backward\_T\_PH\_Region3b**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Backward__T__PH__Region3b.md) _Coefficients and exponents of the backward equation_ \(T_{3a}(p,h)\) _for subregion 3b in its dimensionless form, Eq. (3) in IF97-Region3._    
* **struct** [**Coeff\_Boundary3ab**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Boundary3ab.md)     
* **struct** [**Coeff\_Region1**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region1.md) _Numerical values of the coefficients and exponents of the dimensionless Gibbs free energy for region 1, Eq. (7) of IF97._     
* **struct** [**Coeff\_Region2**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region2.md) _Numerical values of the coefficients and exponents of the ideal-gas part_ \(\gamma^o\) _o of the dimensionless Gibbs free energy for region 2, Eq. (16)a, and the residual part_\(\gamma^r\) _._    
* **struct** [**Coeff\_Region5**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Coeff__Region5.md) _Numerical values of the coefficients and exponents of the residual part  r of the dimensionless Gibbs free energy for region 5. See Table 37, 38 in IF97._     
* **struct** [**Table34**](structxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97_1_1Table34.md) _Table 34 of IF97 for saturation curve calculation._     

