

# File test.cxx

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**NaCl**](dir_99a45c5f26c2a3d925e70f340945fe86.md) **>** [**test.cxx**](NaCl_2test_8cxx.md)

[Go to the documentation of this file](NaCl_2test_8cxx.md)


```C++

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
```


