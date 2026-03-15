

# File NaCl.h

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**NaCl**](dir_99a45c5f26c2a3d925e70f340945fe86.md) **>** [**NaCl.h**](NaCl_8h.md)

[Go to the documentation of this file](NaCl_8h.md)


```C++

#ifndef NaCl_xThermal_H
#define NaCl_xThermal_H

#include "IAPS84.h"
#include "IAPWS95_CoolProp.h"
#include "IAPWS95.h"

namespace xThermal
{
    namespace NaCl
    {
        double const M = 0.058443; 
        double const T_Triple_C = 800.7; 
        double const T_Triple = T_Triple_C + 273.15; 
        double const P_Triple = 50; 
        double const log10_P_Triple = log10(P_Triple);
        double const inv_T_Triple = 1.0/T_Triple;
        double const T_Triple_sqr = T_Triple*T_Triple;
        double const T_Triple_sqr_C = T_Triple_C*T_Triple_C;
        double const slope_Clapeyron = 1.0/2.4726e-7;
        struct H_halite_ref
        {
            const double H = 9.415867359e4; 
            const double T_C = 100;   
            const double P = 100E5;   
            double R0;    
            double R1; 
            double R2; 
            double R3; 
        };
        struct Coeff_Rho
        {
            const double l[6] = {2170.4, -0.24599, -9.5797E-05, 0.005727E-5, 0.002715E-5, 733.4};
            const double m[6] = {58443, 23.772, 0.018639, -1.9687E-06, -1.5259E-10, 5.5058E-13};
        };

        struct Coeff_H
        {
            double r[5] = {1148.81, 0.275774, 8.8103E-05, 
                            0,
                            5.29063E-18 - 9.63084E-21 + 6.50745E-23};
            const double r3[3] = {-0.0017099E-5, - 3.82734E-11, - 8.65455E-14};
        };
        
        class xTHERMO_VAR cNaCl : public cxThermal
        {
        // data
        private:
            cxThermal* m_Water;
            Coeff_Rho m_coeff_rho;
            Coeff_H m_coeff_H;
            H_halite_ref m_H_h_ref;
            std::string m_backendname;
            CONSTENTS_Thermo m_constants;
        // methods
        private:
            void initialize_data();
            void get_q1q2_Tstar_H(const double& P, double& q1, double& q2);
            double Tstar_H(const double& q1, const double& q2, const double& T);
            double DeltaH_fus(const double& T, const double& P);
            double _H1(const double& T_C, const double& P);
            void _calc_R_H1(H_halite_ref& coeff_haliteH_ref);
            double _H2(const double& T_C, const double& P);
            void _printCoeffs_haliteH();
            double _H_Solid(const double& T, const double& P);
        public:
            cNaCl(std::string backend_H2O);
            cNaCl(const cNaCl& salt);
            ~cNaCl();
            // xThermal::BasicState m_state;
        public:
            std::string name(){return "NaCl"; };
            std::string name_backend();
            double Tmin(){return m_constants.Tmin; };                  
            double Tmax(){return m_constants.Tmax;};                          
            double pmin(){return m_constants.pmin; };                             
            double pmax(){return m_constants.pmax;};                  
            double Ttriple(){return m_constants.Ttriple;};               
            double T_critical(){throw NotImplementedError("T_critical is not available for NaCl");};            
            double p_critical(){throw NotImplementedError("T_critical is not available for NaCl");};            
            double rhomass_critical(){throw NotImplementedError("T_critical is not available for NaCl");};      
            double molar_mass(){return m_constants.molar_mass;};       
            void UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);

        public:
            double Boiling_p(const double& T);
            double Sublimation_p(const double& T);
            double Melting_T(const double& P);
            double Melting_T_C(const double& P);
            double Melting_p(const double& T);
            double Boiling_T(const double& P);
            double Sublimation_T(const double& P);
            double P_Vapor(const double& T);
            double T_Vapor(const double& P);
            double Rho_Solid(const double& T, const double& P);
            void Rho_Solid(const double& T, const double& P, double& rho, double& dRhodP, double& dRhodT, double& kappa, double& beta);
            double Rho_Liquid(const double& T, const double& P);
            void Rho_Liquid(const double& T, const double& P, double& rho, double& dRhodP, double& dRhodT, double& kappa, double& beta);
            double H_Liquid(const double& T, const double& P);
            double H_Solid(const double& T, const double& P);
            double Cp_Solid(const double& T, const double& P);
            double Cp_Liquid(const double& T, const double& P);
        public: // vector version
            void Boiling_p(const std::vector<double> T, std::vector<double>& res);
            void Melting_T(const std::vector<double> P, std::vector<double>& res);
            void Sublimation_p(const std::vector<double> T, std::vector<double>& res);
            void Melting_p(const std::vector<double> T, std::vector<double>& res);
            void Boiling_T(const std::vector<double> P, std::vector<double>& res);
            void Sublimation_T(const std::vector<double> P, std::vector<double>& res);
            void T_Vapor(const std::vector<double> P, std::vector<double>& res);
            void P_Vapor(const std::vector<double> T, std::vector<double>& res);
            void Cp_Solid(const std::vector<double> T, const std::vector<double> P, std::vector<double>& res);
            void DeltaH_fus(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res);
            void H_Solid(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res);
            void H_Liquid(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res);
        };
    };

};

#endif
```


