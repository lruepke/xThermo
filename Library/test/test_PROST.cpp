#include <iostream>
#include "thermo.h"
#include "steam4.h"
#include "iaps.h"
#include <cstdio>

#include "IAPWS95.h" // used to compare with PROST
#include "IAPWS95_CoolProp.h"

void test_basicTP();
void test_basicHP();
void test_derivative();
int main()
{
    // test_basicTP();
    // test_basicHP();
    test_derivative();
    return 0;
}
void test_derivative()
{
    double T = 373.15, P = 300E5;
    //---------------------- PROST calculation ----------------------
    double d, dp, ds, dh;
    dp = 1.0e-8;
    ds = 1.0e-8;
    dh = 1.0e-8;
    int dev = 2; // first derivative
    Prop *prop = newProp('t', 'p', dev);
    d = 0.0;
    water_tp(T,P,d,dp,prop);
    double dPdT = prop->dp->T_Cd, dPdRho = prop->dp->d_CT, rho = prop->d;
    double dRhodT = -dPdT/dPdRho;
    // -------- print result of RPOST ------
    // dumpProp(stdout, prop);
    printf("===== Derivative test of PROST in TP space =====\n");
    printf("T=%f deg.C, P=%f bar\n- Rho = %f kg/m3\n- Cp = %f J/kg/K\n- H = %f MJ/kg\n- Mu = %f Pa s\n", T-273.15, P/1E5, prop->d, prop->cp, prop->h/1E6, viscos(prop));
    printf("Derivatives: \n");
    printf("- dPdT: %f\n", dPdT);
    printf("- dPdd: %f\n", dPdRho);
    printf("- Isothermal compressibility: %E\n", 1.0/(rho * dPdRho));
    printf("- Volume expansivity: %E\n", -dRhodT/rho);
    // -----------------------------------------------------------------

    // ---------------------- IAPWS95 calculation ----------------------
    xThermal::IAPWS95::cIAPWS95 iapws95;
    xThermal::ThermodynamicProperties props;
    iapws95.UpdateState_TPX(props, T, P);
    printf("\n===== Derivative test of IAPWS95 in TP space =====\n");
    printf("T=%f deg.C, P=%f bar\n- Rho = %f kg/m3\n- Cp = %f J/kg/K\n- H = %f MJ/kg\n- Mu = %f Pa s\n", T-273.15, P/1E5, props.Rho, props.Cp, props.H/1E6, props.Mu);
    printf("Derivatives: \n");
    // printf("- dPdT: %f\n", props.dPdT_Rho);
    // printf("- dPdd: %f\n", props.dPdRho_T);
    printf("- Isothermal compressibility: %E\n", props.IsothermalCompressibility);
    printf("- Volume expansivity: %E\n", props.IsobaricExpansivity);

    // ---------------------- IAPWS95 coolprop calculation ----------------------
#ifdef USE_COOLPROP
    xThermal::COOLPROP::cIAPWS95_CoolProp iapws95_coolprop;
    xThermal::ThermodynamicProperties props_coolProp;
    iapws95_coolprop.UpdateState_TPX(props_coolProp, T, P);
    printf("\n===== Derivative test of IAPWS95-COOLPROP in TP space =====\n");
    printf("T=%f deg.C, P=%f bar\n- Rho = %f kg/m3\n- Cp = %f J/kg/K\n- H = %f MJ/kg\n- Mu = %f Pa s\n", T-273.15, P/1E5, props_coolProp.Rho, props_coolProp.Cp, props_coolProp.H/1E6, props_coolProp.Mu);
    printf("Derivatives: \n");
    printf("- Isothermal compressibility: %E\n", props.IsothermalCompressibility);
    printf("- Volume expansivity: %E\n", props.IsobaricExpansivity);
#endif
    freeProp(prop);

}

void test_basicHP()
{
    double H = 441643.903950, P = 300E5;

    double t, d, dp, ds, dh;
    dp = 1.0e-8;
    ds = 1.0e-8;
    dh = 1.0e-8;
    Prop *prop = newProp('p', 'h', 1);
    d = 0.0;
    water_ph(P,H,t, d,dp,dh,prop);

    printf("===== Basic test of PROST in HP space =====\n");
    printf("H=%f MJ/kg, P=%f bar\n- T = %f deg.C\n- Rho = %f kg/m3\n- Cp = %f J/kg\n- H = %f J\n- Mu = %f Pa s\n\n", H/1E6, P/1E5, prop->T-273.15, prop->d, prop->cp, prop->h, viscos(prop));
    freeProp(prop);
}

void test_basicTP()
{
    double T = 373.15, P = 300E5;
    double d, dp, ds, dh;
    dp = 1.0e-8;
    ds = 1.0e-8;
    dh = 1.0e-8;
    Prop *prop = newProp('t', 'p', 1);
    d = 0.0;
    water_tp(T,P,d,dp,prop);
    printf("===== Basic test of PROST in TP space =====\n");
    printf("T=%f deg.C, P=%f bar\n- Rho = %f kg/m3\n- Cp = %f J/kg/K\n- H = %f MJ/kg\n- Mu = %f Pa s\n\n", T-273.15, P/1E5, prop->d, prop->cp, prop->h/1E6, viscos(prop));
    freeProp(prop);
}