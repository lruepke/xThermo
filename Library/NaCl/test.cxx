/**
 * @file test.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Test implementation of H2O
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "NaCl.h"

using namespace std;
using namespace xThermal;

void test_constants( xThermal::cxThermal* thermo);
void test_prop_TP(xThermal::cxThermal* thermo, double T, double P);

int main()
{
    double T0 = 900+273.15, P0=30E5;
  
    STATUS("Test NaCl: IAPS84");
    NaCl::cNaCl nacl("IAPS84");
    // test_constants(&nacl);
    double Hs = nacl.H_Solid(594.60+273.15, 388.14744434E5);
    exit(0);
    test_prop_TP(&nacl, T0, P0);

    STATUS("Test NaCl: IAPWS95");
    NaCl::cNaCl nacl2("IAPWS95");
    // test_constants(&nacl2);
    test_prop_TP(&nacl2, T0, P0);

    return 0;
}
 
void test_prop_TP(xThermal::cxThermal* thermo, double T, double P)
{
    ThermodynamicProperties props;
    thermo->UpdateState_TPX(props, T, P);
    cout<<"Input T="<<T<<" [K], p="<<P<<" [Pa]"<<endl;
    cout<<"Phase: "<<thermo->phase_name(props.phase)<<endl;
    cout<<"rho [kg/m^3]: "<<props.Rho<<endl;
    cout<<"h [J/kg]: "<<props.H<<endl;
    cout<<"cp [J/kg/K]: "<<props.Cp<<endl;
}

void test_constants( xThermal::cxThermal* thermo)
{
    cout<<"Tmin [K]: "<<thermo->Tmin()<<endl;
    cout<<"Tmax [K]: "<<thermo->Tmax()<<endl;
    cout<<"pmin [Pa]: "<<thermo->pmin()<<endl;
    cout<<"pmax [Pa]: "<<thermo->pmax()<<endl;
    cout<<"T_triple [K]: "<<thermo->Ttriple()<<endl;
    cout<<"Molar mass [kg/mol]: "<<thermo->molar_mass()<<endl;
}