

# File H2ONaCl.h

[**File List**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**H2ONaCl.h**](H2ONaCl_8h.md)

[Go to the documentation of this file](H2ONaCl_8h.md)


```C++

#ifndef H2ONACL_xThermal_H
#define H2ONACL_xThermal_H

#include "IAPS84.h"
#include "IAPWS95_CoolProp.h"
#include "IAPWS95.h"
#include "NaCl.h"
// #include "IAPWS-95.h"

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

namespace xThermal
{
    namespace H2ONaCl
    {
        #define ITERATION_MAX   1000    
        #define TOL_Pressure   1E-4  
        #define TOL_PTRro 1E-10 
        const double T_MIN = 273.16;
        const double P_MIN = 1E5;
        const double P_MAX = 5000E5;
        const double T_MAX = 1273.15;
        const double X_MIN = 0.0;
        const double X_MAX = 1.0;
        const double T_MAX_VLH = 1073.5662157838103;
        const double T_MIN_VLH = 380.912102015643;
        const double T_Peak_VLH = 867.782443; //594.632443 deg.C
        const double P_Peak_VLH = 39014744.433797;
        const double T_Peak_VLH_C = 867.782443 - 273.15; //594.632443 deg.C
        const double P_Peak_VLH_bar = 390.14744433797;
        struct Table4_Driesner2007a_CriticalCurve
        {
            double c[14], cA[11], d[11];
        };

        struct Coeffs_Estimate_CriticalT
        {
            const int num = 9;
            const double coeffs[9] = {3.18924462e+05, -4.83993904e+05, 3.17857136e+05, 
                                    -1.17956949e+05,  2.70680030e+04, -3.93489090e+03,
                                    3.54050856e+02, -1.80362368e+01, 3.98481737e-01};
        };
        
        struct Params_InvCriticalT;
        struct Params_Inversion_PTX;
        
        struct Table6_Pressure_VLH
        {   
            const int num = 11;
            // 0.00464, 5E-07, 16.9078, -269.148, 7632.04, -49563.6, 233119.0, -513556.0, 549708.0, -284628.0
            const double f[11] = {464.0, 0.05, 1690780.0, -26914800.0, 763204000.0, -4956360000.0, 233119.0E5, -513556.0E5, 549708.0E5, -284628.0E5, 
                                NaCl::P_Triple - (464.0 + 0.05 + 1690780.0 - 26914800.0 + 763204000.0 - 4956360000.0 + 233119.0E5 - 513556.0E5 + 549708.0E5 - 284628.0E5)};
        };
        struct Table7_XL_VL
        {
            // original
            //double h[11]={0.00168486, 0.000219379, 438.58, 18.4508, -5.6765E-10, 6.73704E-06,1.44951E-07, 384.904, 7.07477, 6.06896E-05, 0.00762859};
            const double h[11] = {1.68486E-8, 2.19379E-9, 438.58, 18.4508, -5.6765E-15,6.73704E-16, 1.44951E-17, 384.904, 7.07477, 6.06896E-15, 0.00762859};
        };
        // this j[4] array can not be in this struct, otherwise parallel computing will give wrong result, because member variable is default as shared data
        struct Table8_VaporComposition {
            const double k[16] = {-0.235694, -0.188838, 0.004, 0.0552466, 0.66918, 396.848, 45, -3.2719E-07,
                                  141.699, -0.292631, -0.00139991, 1.95965E-06, -7.3653E-10, 0.904411, 0.000769766,
                                  -1.18658E-06};
            //j[4] will be calculated it in Log10_Kprime function
        };

        struct Coeffs_Viscosity
        {
            const double a1 = -35.9858;
            const double a2 = 0.80017;
            const double b1 = 1e-6;
            const double b2 = -0.05239;
            const double b3 = 1.32936;
        };

        // struct ThermodynamicProperty
        // {
        //     State state;
        //     double S_l, S_v, S_h;
        //     double X_l, X_v;
        //     double Rho_l, Rho_v, Rho_h;
        //     double H_l, H_v, H_h;
        //     double Rho;
        // };

        class xTHERMO_VAR cH2ONaCl: public cxThermal
        {
        private:
            cxThermal* m_Water;
            NaCl::cNaCl* m_NaCl;
            CONSTENTS_Thermo m_constants_Water;
            Table4_Driesner2007a_CriticalCurve m_tab4;
            Coeffs_Estimate_CriticalT m_coeff_estimate_Tcrit;
            Table6_Pressure_VLH m_tab6;
            Table7_XL_VL m_table7;
            Table8_VaporComposition m_tab8;
            Coeffs_Viscosity m_coeffs_mu;
            std::string m_backendname;
        private:
            void createTable4_Driesner2007a(Table4_Driesner2007a_CriticalCurve& table4);
            void initialize_data();
            void X_Critical_mol(double T, double& X_crit);
            void T_Critical_Bisection(double P, double& T_crit);
            void T_Critical_estimate(double P, double& T_crit);
            double X_HaliteLiquidus_mol(const double& T, const double& P);
            double T_VLH_P0_(double P0, double Tmin, double Tmax);
            void X_HaliteSat_mol(const double& T, const double& P, double& X_L, double& X_V);
            double Log10_Kprime(const double& T, const double& P, double& P_NaCl);
            double XL_VL_mol(const double& T, const double& P);
            double XV_VL_mol(const double& T, const double& P);
            void _getPhaseRegion_node_HaliteLiquidus(const double& T, const double& P, const double& X, xThermal::PhaseRegion& phase_region, double& X_V, double& X_L);
            void _getPhaseRegion_node_CheckTcrit_H2O(const double& T, const double& P, const double& X, xThermal::PhaseRegion& phase_region, double& X_V, double& X_L);
            void _getPhaseRegion_node_XVXL_VL(const double& T, const double& P, const double& X, const double& X_crit, xThermal::PhaseRegion& phase_region, double& X_V, double& X_L);
            void n1n2_Tstar_V(const double& P, const double& X, double& n1, double& n2);
            double D_Tstar_V(const double& T, const double& P, const double& X);
            double extrapolation_V_highT(const double& T_C, const double& P_bar, const double& X_mol);
            double extrapolation_V_lowPlowT(const double& T_C, const double& P_bar, const double& X_mol);
            double _Rho_water(const double& T, const double& P, double& dRhodP, double& dRhodT, double& IsothermalCompresibility, double& IsobaricExpansivity, PhaseType phase=Liquid);
            double _Rho_water(const double& T, const double& P, PhaseType phase=Liquid);
            void q1q2_Tstar_H(const double& P_bar, const double& X_mol, double& q1, double& q2);
            void extrapolation_H_Cp_highT(const double& T_C, const double& P_bar, const double& X_mol, double& H, double& Cp);
            double extrapolation_H_highT(const double& T_C, const double& P_bar, const double& X_mol); //测试阶段使用
            double _H_water(const double& T, const double& P, PhaseType phase);
            void _H_Cp_water(const double& T, const double& P, double& H, double& Cp, PhaseType phase);
            double _Mu_water(const double& T, const double& P, PhaseType phase);
            double Sl_VL(const double& X_L, const double& X_V, const double& Rho_L, const double& Rho_V, const double& X);
            double Saturation_Phase1(const double& X_Phase1, const double& X_Phase2, const double& Rho_Phase1, const double& Rho_Phase2, const double& X);
            bool SolveSaturation_ThreePhase_HX(const double& H, const double& X, const double X_lvh[3], const double Rho_lvh[3], const double H_lvh[3], double S_lvh[3]);
            void T_VL_bisection(double P, double X, double& T, PhaseType phase);
        public:
            cH2ONaCl(std::string backend_H2O);
            cH2ONaCl(const cH2ONaCl& sw);
            ~cH2ONaCl();
        public:
            cxThermal* get_pWater() const {return m_Water;};
            double Wt2Mol(const double& X_wt);
            void Wt2Mol(const double& X_wt, double& X_mol);
            double Mol2Wt(double X_mol);
            // phase boundary related functions
            // 1. critical curve: usually the minimum valid T is T_critical of H2O
            void X_Critical(double T, double& X_crit);
            void P_Critical(double T, double& P_crit);
            void T_Critical(double P, double& T_crit);
            void P_X_Critical(double T, double& P_crit, double& X_crit);
            void T_X_Critical(double P, double& T_crit, double& X_crit);
            // 2. Pressure of VLH surface for a given T, the valid T in range [T_MIN_VLH, T_MAX_VLH]
            void P_VLH(const double& T, double& P);
            double P_VLH(const double& T);
            double dPdT_VLH(const double& T);
            double Tmax_VLH();
            double Tmin_VLH();
            void T_VLH_P0(double P0, double& T_min, double& T_max);
            void X_VLH(const double& T, const double& P, double& X_L, double& X_V);
            // 2. Halite liquidus, it is actually the same as X_L calculated from X_VLH
            double X_HaliteLiquidus(const double& T, const double& P);
            // 3. Halite saturated vapor composition in VH zone, this [TT,PP,XX] surface is the "left" bound of VH zone along X-axis direction.
            // The valid pressure depends on the given T, it is [pmin, P_VLH(T)]
            double X_VH(const double& T, const double& P);
            // 4. Saturated liquid & vapor composition on VL surface, [TT,PP,XV_VL] is the "left" bound and [TT,PP,XL_VL] is the "right" bound of the VL region along X-axis direction
            double XL_VL(const double& T, const double& P);
            double T_VL(const double& P, const double& X, PhaseType phase);
            double T_VL_L(const double& p, const double& X, const double& T_low, const double& T_high);
            double XV_VL(const double& T, const double& P);
            double T_VL_V(const double& p, const double& X, const double& T_low, const double& T_high);
            void X_VL(const double& T, const double& P, double& X_L, double& X_V);
            void prop_VL(const double& T, const double& P, const double& X, ThermodynamicProperties& prop);
            // 5. calculate phase region and related X_V, X_L for given T, p, X
            void getPhaseRegion_TPX(const double& T, const double& P, const double& X, xThermal::PhaseRegion& phase_region, double& X_V, double& X_L);
            // 6. Calculate density, specific enthalpy and viscosity for given T, P, X. The argument phase [Liquid, Vapor] is used to determine how whether to extrapolate.
            void Rho_phase(const double& T, const double& P, const double& X, double& rho, double& dRhodP, double& dRhodT, PhaseType phase);
            void H_phase(const double& T, const double& P, const double& X, double& H, PhaseType phase);
            void H_phase(const double& T, const double& P, const double& X, double& H, double& Cp, PhaseType phase);
            void Mu_phase(const double& T, const double& P, const double& X, double& Mu, PhaseType phase);
            // calculate T for given H,P,X
            double T_HPX(const double& H, const double& p, const double& X, const double& T_low=T_MIN, const double& T_high=T_MAX);
            // calculate phase region for given H,P,X
            void getPhaseRegion_HPX(const double& H, const double& p, const double& X, xThermal::PhaseRegion& phase_region, double& T, double S_lvh[3]);
            void compressibility_VL(const double& T0, const double& p0, const double& X0, const double& bulkRho, double& compressibility, double dp=1);
            void compressibility_VH(const double& T0, const double& p0, const double& X0, const double& bulkRho, double& compressibility, double dp=1);
            void compressibility_LH(const double& T0, const double& p0, const double& X0, const double& bulkRho, double& compressibility, double dp=1);
            void HminHmax_VLH_triangle(const double& H_v0, const double& H_l0, const double& H_h0, const double& X_v0, const double& X_l0, const double& X0, double& Hmin0, double& Hmax0);
        public:
            std::string name(){return "H2O-NaCl"; }
            std::string name_backend();
            // Thermodynamic constants
            double Tmin(){return std::max(T_MIN, m_Water->Tmin()); };                
            double Tmax(){return T_MAX;};                            
            double pmin(){return P_MIN; };                            
            double pmax(){return P_MAX;};                          
            double Ttriple(){throw NotImplementedError("Ttriple is not available for H2O-NaCl");};               
            double T_critical(){throw NotImplementedError("T_critical is not available for H2O-NaCl");};            
            double p_critical(){throw NotImplementedError("p_critical is not available for H2O-NaCl");};            
            double rhomass_critical(){throw NotImplementedError("rhomass_critical is not available for H2O-NaCl");};      
            double molar_mass(){throw NotImplementedError("molar_mass is not available for H2O-NaCl");};       
            // Update thermodynamic state and properties for given T,P,X
            PhaseRegion findPhaseRegion_TPX(const double& T, const double& p, const double& X);
            // - Single point mode
            void UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X);
            void UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X);
            bool UpdateState_HPX_vlh(ThermodynamicProperties& props, const double& H, const double& p, const double& X, const std::vector<double>& T1_T2);
            bool UpdateState_HPX_vl_lowXlowP(ThermodynamicProperties& props, const double& H, const double& p, const double& X);
        // Utility function
        public:
            // 7. Calculate phase boundary surface, usually used to make plot in python. In two formats, (1) deformed-linear mesh[TT,PP,XX], which can be plotted using ax.plot_surface
            DeformLinearMesh PhaseBoundary_HaliteLiquidus_DeformLinear(double p_max=2500E5, double dT = 10, double dp = 50E5);
            DeformLinearMesh PhaseBoundary_VL_DeformLinear(PhaseType VaporOrLiquid,int nT=100, int np=200);
            DeformLinearMesh PhaseBoundary_VLH_DeformLinear(int nT=100, int nX=60);
            DeformLinearMesh PhaseBoundary_VH_DeformLinear(int nT=100, int np=100);
            // 8. call function writePhaseBoundaries2VTU can save the phase boundary DeformLinearMesh to vtu file directly.
            PhaseBoundaries calc_PhaseBoundaries(std::string scale_X = "linear", double ratio_log_to_linear = 1, double Xcenter=0.01);
            // provide triangular mesh as well, the triangular mesh is generated using the Triangle code.
            TriMesh PhaseBoundary_HaliteLiquidus(std::string fmt_out="",double p_max=2500E5, double dT = 10, double dp = 50E5);
            TriMesh Triangulation(const std::vector<double>& x_poly, const std::vector<double>& y_poly, double xIn, double yIn, double dx, double dy);
            PhaseRegion_Slice Slice_constP(const double P0, size_t nPoints=500);
            PhaseRegion_Slice Slice_constT(const double T0, size_t nPoints=500);
        public: //vector version for python API
            std::vector<double> Mol2Wt(std::vector<double> X_mol);
            std::vector<double> Wt2Mol(std::vector<double> X_wt);
            void P_X_Critical(std::vector<double> T, std::vector<double>& P_crit, std::vector<double>& X_crit);
            void P_Critical(std::vector<double> T, std::vector<double>& P_crit);
            void T_Critical(std::vector<double> P, std::vector<double>& T_crit);
            void T_X_Critical(std::vector<double> P, std::vector<double>& T_crit, std::vector<double>& X_crit);
            void X_HaliteLiquidus(std::vector<double> T, std::vector<double> P, std::vector<double>& res);
            void P_VLH(std::vector<double> T, std::vector<double>& P);
            void X_VH(std::vector<double> T, std::vector<double> P, std::vector<double>& res);
            void X_VLH(std::vector<double> T, std::vector<double> P, std::vector<double>& X_L, std::vector<double>& X_V);
            void XL_VL(std::vector<double> T, std::vector<double> P, std::vector<double>& res);
            void XV_VL(std::vector<double> T, std::vector<double> P, std::vector<double>& res);
            void X_VL(std::vector<double> T, std::vector<double> P, std::vector<double>& X_L, std::vector<double>& X_V);
            std::vector<int> getPhaseRegion_TPX(std::vector<double> T, std::vector<double> P, std::vector<double> X, std::vector<double>& X_V, std::vector<double>& X_L);
            std::vector<double> Rho_phase(std::vector<double> T, std::vector<double> P, std::vector<double> X, PhaseType phase);
            std::vector<double> H_phase(std::vector<double> T, std::vector<double> P, std::vector<double> X, PhaseType phase);
            double H_phase(const double& T, const double& P, const double& X, PhaseType phase){double H; H_phase(T,P,X,H,phase); return H;};
            double Rho_phase(const double& T, const double& P, const double& X, PhaseType phase){double rho, dRhodP, dRhodT; Rho_phase(T,P,X,rho, dRhodP, dRhodT,phase); return rho;};
            void n1n2_Tstar_V(std::vector<double> P, std::vector<double> X, std::vector<double>& n1, std::vector<double>& n2);
            std::vector<double> D_Tstar_V(std::vector<double> T, std::vector<double> P, std::vector<double> X);
            void q1q2_Tstar_H(std::vector<double> P, std::vector<double> X, std::vector<double>& q1, std::vector<double>& q2);
        // public:  //for python wrapper
        //     ThermodynamicProperties UpdateState_TPX(const double& T, const double& p, const double& X);
        //     ThermodynamicProperties UpdateState_HPX(const double& H, const double& p, const double& X);
        };

        struct Params_P2CriticalT
        {
            cH2ONaCl* sw; 
            double P; 
        };
        struct Params_Inversion_PTX
        {
            cH2ONaCl* sw; 
            double P;
            double X;
            double H;
        };
        struct Params_PX2T_VL
        {
            cH2ONaCl* sw; 
            double P; 
            double X; 
            PhaseType phase; 
        };
        
    };

};

#endif
```


