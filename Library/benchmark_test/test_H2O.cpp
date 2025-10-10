#include "H2O.h"
#include "IAPWS-95.h"
#include "IAPWS-IF97.h"
#ifdef IAPWS_OTHERS
    #include "freesteam_class.h"
    #include "PROST_class.h"
    #include "CoolProp_class.h"
#endif
using namespace std;

void test_Boiling_P(double T_K);
void test_Boiling_T(double P);
double test_TP2Rho(double T_K, double P_Pa);
/**
 * @brief Compare and validate result given in Table 7 of \cite IAPWS-95
 * 
 */
void validate_SinglePhaseRegion_Table7();
/**
 * @brief Compare and validate result given in Table 8 of \cite IAPWS-95
 * 
 */
void validate_TwoPhaseRegion_Table8();
/**
 * @brief Validate ideal-gas part \f$ \phi^o \f$ and residual part \f$ \phi^r \f$ of the dimensionless Helmholtz energy and its derivatives for T=500 K and \f$ \rho = 838.025 kg/m^3\f$. See Table 6 of \cite IAPWS-95
 * 
 */
void validate_phio_phir_Table6();

void test_H(double T_K, double P);

void test_state_HP(double H, double P);
void test_state_HP_2D();
int main(int argc, char** argv)
{
    // IAPWS95::cIAPWS95 iapws95;
    // printf("Pc = %.20E\n",iapws95.Boundary_region3ab_P2H(H2O::P_c));
    // test_Boiling_P(450);
    // test_Boiling_T(0.169082693E8);
    // test_TP2Rho(300, 0.992418352E5);
    // validate_SinglePhaseRegion_Table7();
    // validate_TwoPhaseRegion_Table8();
    // validate_phio_phir_Table6();
    // test_H(200 + 273.15, 300E5);
    // test_state_HP(3900000.000, 9100000.000);
    // test_state_HP(800000.000, 29100000.000);
    // test_state_HP(700000.000, 18100000.000);
    // test_state_HP(2500000.000, 19100000.000);
    // test_state_HP(1350000.000, 8700000.000);
    // test_state_HP(CONST_IF97_H_c, H2O::P_c);
    // test_state_HP(CONST_IF97_H_c, 120E6);
    // test_state_HP(6E6, 120E6);
    // test_state_HP(0.4E6, 120E6);
    // test_state_HP(200000.000, 358586031.862);
    test_state_HP(1509418.838, 12450741.483);
    // test_state_HP_2D();
    return 0;
}

void test_state_HP(double H, double P)
{
    using namespace std;
    IAPWS95::cIAPWS95 iapws95;
    IAPWS95::State state;
    iapws95.UpdateState_HP(H, P);
    state = iapws95.getState();
    double h, hl, hv;
    iapws95.h(h, hl, hv);
    double cp, cpl, cpv;
    iapws95.cp(cp, cpl, cpv);
    printf("---- H = %f MJ/kg, P = %f bar\n", H/1E6, P/1E5);
    printf("Phase region: %d, T = %.2f K, cp = %.2f J/kg/K, rho = %.2f kg/m3, x = %.2E\n", state.phase, state.T.value, cp, state.Rho.value, state.x.value);
    #ifdef IAPWS_OTHERS
        FREESTEAM::cFreeSteam steam;
        steam.UpdateState_HP(H, P);
        printf("Freesteam, T = %.2f K, cp = %.2f J/kg/K, rho = %.2f kg/m3\n",steam.T(), steam.cp(), steam.rho());
        PROST::cPROST prost;
        prost.UpdateState_TP(steam.T(), P);
        printf("PROST, H = %f MJ/kg, rho = %f kg/m3\n", prost.h()/1E6, prost.rho());
        COOLPROP::cCoolProp coolprop;
        coolprop.UpdateState_HP(H, P);
        printf("CoolProp, T = %.2f K, cp = %.2f J/kg/K, rho = %.2f kg/m3\n",coolprop.T(), coolprop.cp(), coolprop.rho());
    #endif
    // printf("Phase region: %d, T = %.2f K, x = %.2E, hl = %f, hv=\n", state.phase, state.T.value, state.x.value);
    // WAIT("TEST");
}
void test_state_HP_2D()
{
    using namespace std;
    IAPWS95::cIAPWS95 iapws95;
    IAPWS95::State state;
    double H, P;
    // test 2D
    double Hmin = 0.1E6, Hmax = 4.5E6, Pmin = H2O::P_MIN, Pmax = 200E6;
    double dh=0.01E6, dp = 1E5;
    int np = int((Pmax-Pmin)/dp);
    int nh = int((Hmax - Hmin)/dh);
    for (int i = 0; i < np; i++)
    {
        // printf("%d\n", i);
        P = Pmin + i*dp;
        for (int j = 0; j < nh; j++)
        {
            H = Hmin + j*dh;
            // iapws95.UpdateState_HP(H , P, state);
        }
        
    }
}

void test_H(double T_K, double P)
{
    using namespace std;
    IAPWS95::cIAPWS95 iapws95;
    double h = iapws95._enthalpy(T_K, P);
    printf("T = %.2f K, P = %.2f bar, H = %.2f J/kg\n", T_K, P, h);
}

void validate_phio_phir_Table6()
{
    using namespace std;
    IAPWS95::cIAPWS95 iapws95;
    double T = 500;
    double rho = 838.025;
    double delta = rho/H2O::Rho_c;
    double tau = H2O::T_c/T;
    IAPWS95::HelmholtzEnergy_dimensionless phio_ = {0.204797733E1, 0.384236747, -0.147637878, 0.904611106E1, -0.193249185E1, 0};
    IAPWS95::HelmholtzEnergy_dimensionless phir_ = {-0.342693206E1, -0.364366650, 0.856063701, -0.581403435E1, -0.223440737E1, -0.112176915E1};
    IAPWS95::HelmholtzEnergy_dimensionless phio, phir;
    iapws95.phi_o(delta, tau, phio, Update_phi_all);
    phir.value = iapws95.phi_r(delta, tau);
    phir.d = iapws95.phi_r_d(delta, tau);
    phir.dd = iapws95.phi_r_dd(delta, tau);
    phir.t = iapws95.phi_r_t(delta, tau);
    phir.tt = iapws95.phi_r_tt(delta, tau);
    phir.dt = iapws95.phi_r_dt(delta, tau);
    printf("======= T = %.0f K, rho = %.3f kg/m3 ========\n", T, rho);
    if(phio==phio_)
    {
        cout<<"phio: "<<COLOR_GREEN<<"pass"<<COLOR_DEFAULT<<endl;
    }else
    {
        cout<<"phi0: "<<COLOR_RED<<"failed"<<COLOR_DEFAULT<<endl;
        cout<<"phio_cal: "<<phio<<endl;
        cout<<"phio_tab: "<<phio_<<endl;
        cout<<"err:      "<<phio - phio_<<endl;
    }
    if(phir==phir_)
    {
        cout<<"phir: "<<COLOR_GREEN<<"pass"<<COLOR_DEFAULT<<endl;
    }else
    {
        cout<<"phir: "<<COLOR_RED<<"failed"<<COLOR_DEFAULT<<endl;
        cout<<"phir_cal: "<<phir<<endl;
        cout<<"phir_tab: "<<phir_<<endl;
        cout<<"err:      "<<phir - phir_<<endl;
    }
    printf("===============================================\n");

}

void validate_TwoPhaseRegion_Table8()
{
    std::vector<double> T_K     = {275, 450, 625};
    std::vector<double> P       = {0.698451167E3,0.932203564E6,0.169082693E8};
    std::vector<double> rho_l   = {0.999887406E3,0.890341250E3,0.567090385E3};
    std::vector<double> rho_v   = {0.550664919E-2,0.481200360E1,0.118290280E3};
    // calculate and compare rho
    IAPWS95::cIAPWS95 iapws95;
    double err, p_cal, rho_l_cal, rho_v_cal, T_cal;
    printf("=========== Validate Boiling_P =========\n");
    for (size_t i = 0; i < T_K.size(); i++)
    {
        iapws95.Boiling_P(T_K[i], p_cal, rho_l_cal, rho_v_cal); 
        err = (fabs(p_cal - P[i])/p_cal + fabs(rho_l_cal - rho_l[i])/rho_l_cal + fabs(rho_v_cal - rho_v[i])/rho_v_cal)/3;
        if(err<1E-6)
        {
            printf("Two phase validation %ssame%s: T = %.1f, P = %.8E, rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_GREEN, COLOR_DEFAULT, T_K[i], p_cal, rho_l_cal, rho_v_cal, fabs(p_cal - P[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }else if(err<1E-3)
        {
            printf("Two phase validation %sgood%s: T = %.1f, P = %.8E, rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_BLUE, COLOR_DEFAULT, T_K[i], p_cal, rho_l_cal, rho_v_cal, fabs(p_cal - P[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }else
        {
            printf("Two phase validation %sfailed%s: T = %.1f, P = %.8E, rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_RED, COLOR_DEFAULT, T_K[i], p_cal, rho_l_cal, rho_v_cal, fabs(p_cal - P[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }
    }
    // test Boiling_T
    printf("=========== Validate Boiling_T =========\n");
    for (size_t i = 0; i < P.size(); i++)
    {
        iapws95.Boiling_T(P[i], T_cal, rho_l_cal, rho_v_cal); 
        err = (fabs(T_cal - T_K[i])/T_cal + fabs(rho_l_cal - rho_l[i])/rho_l_cal + fabs(rho_v_cal - rho_v[i])/rho_v_cal)/3;
        if(err<1E-6)
        {
            printf("Two phase validation %ssame%s: P = %.8E, T = %.1f,  rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_GREEN, COLOR_DEFAULT, P[i], T_cal, rho_l_cal, rho_v_cal, fabs(T_cal - T_K[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }else if(err<1E-3)
        {
            printf("Two phase validation %sgood%s: P = %.8E, T = %.1f, rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_BLUE, COLOR_DEFAULT, P[i], T_cal, rho_l_cal, rho_v_cal, fabs(T_cal - T_K[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }else
        {
            printf("Two phase validation %sfailed%s: P = %.8E, T = %.1f, rho_l = %.8E, rho_v: %.8E, err = %.2E, err = %.2E, err = %.2E\n", COLOR_RED, COLOR_DEFAULT, P[i], T_cal, rho_l_cal, rho_v_cal, fabs(T_cal - T_K[i]) , fabs(rho_l_cal - rho_l[i]) , fabs(rho_v_cal - rho_v[i]));
        }
    }
}

void validate_SinglePhaseRegion_Table7()
{
    std::vector<double> T_K     = {300, 300, 300, 500, 500, 500, 500, 647, 900, 900, 900};
    std::vector<double> p_Pa    = {0.992418352E5, 0.200022515E8, 0.700004704E9, 0.999679423E5, 
                                   0.999938125E6, 0.100003858E8, 0.700000405E9, 0.220384756E8, 
                                   0.100062559E6, 0.200000690E8, 0.700000006E9};
    std::vector<double> rho     = {0.9965560E3, 0.1005308E4, 0.1188202E4, 0.4350000, 0.4532000E1, 
                                   0.8380250E3, 0.1084564E4, 0.3580000E3, 0.2410000, 0.5261500E2, 
                                   0.8707690E3};
    // calculate and compare rho
    double err, rho_cal;
    for (size_t i = 0; i < T_K.size(); i++)
    {
        rho_cal = test_TP2Rho(T_K[i], p_Pa[i]);
        err = fabs(rho_cal - rho[i]);
        if(err<1E-6)
        {
            printf("Single phase validation %ssame%s: T = %.1f K, P = %.8E Pa, rho = %.8E kg/m3, err = %.8E kg/m3\n", COLOR_GREEN, COLOR_DEFAULT, T_K[i], p_Pa[i], rho_cal, err);
        }else
        {
            printf("Single phase validation %sgood%s: T = %.1f K, P = %.8E Pa, rho = %.8E kg/m3, err = %.8E kg/m3\n", COLOR_RED, COLOR_DEFAULT, T_K[i], p_Pa[i], rho_cal, err);
        }
    }
    
}

double test_TP2Rho(double T_K, double P_Pa)
{
    IAPWS95::cIAPWS95 iapws95;
    double rho;
    // rho = iapws95.Rho(T_K, P_Pa);
    rho = iapws95.Rho(T_K, P_Pa, "bisection");
    // printf("T=%f C, P = %f bar, Rho = %f kg/m3\n", T_K - 273.15, P_Pa/1E5, rho);
    // rho = iapws95.Rho_fdf(T_K, P_Pa);
    
    return rho;
}

void test_Boiling_P(double T_K)
{
    IAPWS95::cIAPWS95 iapws95;
    double P, rho_l, rho_v;
    iapws95.Boiling_P(T_K, P, rho_l, rho_v); 
    printf("T=%f C, P = %.10f bar, Rho_l = %.10f kg/m3, Rho_v: %.10f kg/m3\n", T_K - 273.15, P/1E5, rho_l, rho_v);
}

void test_Boiling_T(double P)
{
    IAPWS95::cIAPWS95 iapws95;
    double T_K, rho_l, rho_v;
    iapws95.Boiling_T(P, T_K, rho_l, rho_v); 
    printf("P = %.10f bar, T=%f C, Rho_l = %.10f kg/m3, Rho_v: %.10f kg/m3\n", P/1E5, T_K - 273.15, rho_l, rho_v);
}