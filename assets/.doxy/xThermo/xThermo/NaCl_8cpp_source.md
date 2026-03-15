

# File NaCl.cpp

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**NaCl**](dir_99a45c5f26c2a3d925e70f340945fe86.md) **>** [**NaCl.cpp**](NaCl_8cpp.md)

[Go to the documentation of this file](NaCl_8cpp.md)


```C++

#include "NaCl.h"
// ----- LUT related head filess --------
#include "LookUpTableForestI.H"
#include "interpolationI.H"
#include "AMR_LUT_RefineFuncI.H"
//---------------------------------------

namespace xThermal
{
    namespace NaCl
    {
        void cNaCl::initialize_data()
        {
            m_constants.Tmin = m_Water->Tmin();
            m_constants.Tmax = 1273.15;
            m_constants.pmin = 1E5;
            m_constants.pmax = m_Water->pmax();
            m_constants.Ttriple = T_Triple;
            // m_constants.T_critical = 0;
            // m_constants.p_critical = 0;
            // m_constants.rhomass_critical = 0;
            m_constants.molar_mass = M;
            // calculate coefficients for halite enthalpy
            _calc_R_H1(m_H_h_ref);
//            _printCoeffs_haliteH();//only call this function at the paper writing stage.
        }
        cNaCl::cNaCl(std::string backend_H2O)
        :m_Water(NULL)
        {
            m_backendname = backend_H2O;
            if(m_backendname==Name_Backend_IAPS84)
            {
                m_Water = new PROST::cIAPS84;
            }
#ifdef USE_COOLPROP
            else if(m_backendname==Name_Backend_IAPWS95_CoolProp)
            {
                m_Water = new COOLPROP::cIAPWS95_CoolProp;
            }
#endif
            else if(m_backendname==Name_Backend_IAPWS95)
            {
                m_Water = new IAPWS95::cIAPWS95;
            }
            else
            {
                std::string available_backends = "'IAPS84', 'IAPWS95'";
#ifdef USE_COOLPROP
                available_backends = "'IAPS84', 'IAPWS95' and 'IAPWS95_CoolProp'";
#endif
                throw NotImplementedError("Error in cNaCl(std::string backend_H2O). The supported H2O backend is only one of " + available_backends + ". Input name is " + m_backendname);
            }
            // initialize phase
            initialize_data();
        }
        // copy constructor
        cNaCl::cNaCl(const cNaCl& salt)
        {
            if(salt.m_backendname==Name_Backend_IAPS84)
            {
                m_Water = new PROST::cIAPS84;
            }
#ifdef USE_COOLPROP
            else if(salt.m_backendname==Name_Backend_IAPWS95_CoolProp)
            {
                m_Water = new COOLPROP::cIAPWS95_CoolProp;
            }
#endif
            else if(m_backendname==Name_Backend_IAPWS95)
            {
                m_Water = new IAPWS95::cIAPWS95;
            }else
            {
                throw NotImplementedError("Error in cNaCl(const cNaCl& salt). The supported H2O backend is only one of 'IAPS84', 'IAPWS95' and 'IAPWS95_CoolProp'. ");
            }
            // initialize phase
            initialize_data();
            // copy
            // m_state = salt.m_state;
        }

        cNaCl::~cNaCl()
        {
            if(m_Water)delete m_Water;
        }

        void cNaCl::UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X)
        {
            props.fluidName = name();
            props.T = T;
            props.p = p;
            double T_melting = Melting_T(p);
            if (T < T_melting)
            {
                props.phase = SinglePhase_S;
                Rho_Solid(T,p, props.Rho, props.dRhodP, props.dRhodT, props.IsothermalCompressibility, props.IsobaricExpansivity);
                props.H = H_Solid(T,p);
                props.Cp = Cp_Solid(T,p);
            }else
            {
                props.phase = SinglePhase_L;
                Rho_Liquid(T,p, props.Rho, props.dRhodP, props.dRhodT, props.IsothermalCompressibility, props.IsobaricExpansivity);
                props.H = H_Liquid(T,p);
                props.Cp = Cp_Liquid(T,p);
            }
        }
        
        void cNaCl::get_q1q2_Tstar_H(const double& P, double& q1, double& q2)
        {
            double P_sqr = P * P;

            q1 = 47.9048 - 9.36994E-8 * P + 6.51059E-16 * P_sqr;
            q2 = 0.241022 + 3.45087E-10 * P - 4.28356E-19 * P_sqr;
        }

        double cNaCl::Tstar_H(const double& q1, const double& q2, const double& T)
        {
            return q1 + q2*(T-273.15) + 273.15;
        }

        double cNaCl::Boiling_p(const double& T)
        {
            return pow(10.0, (log10_P_Triple + 9418.12 * (1.0 / T_Triple - 1.0 /T)));
        }

        double cNaCl::Boiling_T(const double& P)
        {
            return 1.0/(inv_T_Triple - (log10(P) - log10_P_Triple)/9418.12);
        }

        double cNaCl::Melting_T(const double& P)
        { 
            return T_Triple + 2.4726e-7*(P - P_Triple);
        }

        double cNaCl::Melting_T_C(const double& P)
        {
            return T_Triple_C + 2.4726e-7*(P - P_Triple);
        }

        double cNaCl::Melting_p(const double& T)
        {
            return P_Triple + 1.0/2.4726e-7 * (T - T_Triple);
        }

        double cNaCl::Sublimation_p(const double& T)
        {
            return pow(10, (log10_P_Triple + 11806.1 * (1.0 / T_Triple - 1.0 / T)));
        };

        double cNaCl::Sublimation_T(const double& P)
        {
            return 1.0/(inv_T_Triple - (log10(P) - log10_P_Triple)/11806.1);
        }

        double cNaCl::P_Vapor(const double& T)
        {
            return T<T_Triple ? Sublimation_p(T) : Boiling_p(T);
        }
        
        double cNaCl::T_Vapor(const double& P)
        {
            return P<P_Triple ? Sublimation_T(P) : Boiling_T(P); 
        }

        double cNaCl::Rho_Solid(const double& T, const double& P)
        {
            double T_C = T - 273.15;

            return m_coeff_rho.l[0] + m_coeff_rho.l[1]*T_C +m_coeff_rho.l[2]*T_C*T_C + (m_coeff_rho.l[3] + m_coeff_rho.l[4]*exp(T_C/m_coeff_rho.l[5]))*P;
        }

        void cNaCl::Rho_Solid(const double& T, const double& P, double& rho, double& dRhodP, double& dRhodT, double& kappa, double& beta)
        {
            double T_C = T - 273.15;
            double rho0 = m_coeff_rho.l[0] + m_coeff_rho.l[1]*T_C +m_coeff_rho.l[2]*T_C*T_C;
            double exp_T_l5 = exp(T_C/m_coeff_rho.l[5]);
            double l = m_coeff_rho.l[3] + m_coeff_rho.l[4]*exp_T_l5;
            dRhodP = l;
            dRhodT = m_coeff_rho.l[1] + 2.0*m_coeff_rho.l[2]*T_C + m_coeff_rho.l[4]/m_coeff_rho.l[5]*exp_T_l5*P;
            rho = rho0 + l*P;
            kappa = 1.0/rho * dRhodP;
            beta = -1.0/rho * dRhodT;
        }
        
        double cNaCl::Rho_Liquid(const double& T, const double& P)
        {
            double T_C = T - 273.15;

            return m_coeff_rho.m[0]/((m_coeff_rho.m[1] + m_coeff_rho.m[2]*T_C + m_coeff_rho.m[3]*T_C*T_C) * (1 - 0.1*log(1 + 10*P*(m_coeff_rho.m[4] + m_coeff_rho.m[5]*T_C))));
        }

        void cNaCl::Rho_Liquid(const double& T, const double& P, double& rho, double& dRhodP, double& dRhodT, double& kappa, double& beta)
        {
            double T_C = T - 273.15;
            double rho0 = m_coeff_rho.m[0]/((m_coeff_rho.m[1] + m_coeff_rho.m[2]*T_C + m_coeff_rho.m[3]*T_C*T_C));
            double kappa_ = (m_coeff_rho.m[4] + m_coeff_rho.m[5]*T_C);
            double one_minus_ln = (1 - 0.1*log(1 + 10*P*kappa_));
            rho = rho0/one_minus_ln;
            double rho_sqre = rho*rho;
            dRhodP = (kappa_ * rho_sqre)/(rho0 * (1 + 10*P * kappa_));
            dRhodT = - rho_sqre * one_minus_ln/m_coeff_rho.m[0] * (m_coeff_rho.m[2] + 2.0*m_coeff_rho.m[3]*T_C);
            kappa = 1.0/rho * dRhodP;
            beta = -1.0/rho * dRhodT;
        }

        double cNaCl::H_Liquid(const double& T, const double& P)
        {
            double q1, q2;
            get_q1q2_Tstar_H(P, q1, q2);
            // cout<<"T*h: "<<Tstar_H(q1, q2, T)<<", P: "<<P/1E5<<endl;
            ThermodynamicProperties props;
            m_Water->UpdateState_TPX(props, Tstar_H(q1, q2, T), P);
            // return m_Water->hmass();
            return props.H;
        }

        double cNaCl::Cp_Liquid(const double& T, const double& P)
        {
            double q1, q2;
            get_q1q2_Tstar_H(P, q1, q2);
            double Tstar = Tstar_H(q1, q2, T);
            ThermodynamicProperties props;
            m_Water->UpdateState_TPX(props, Tstar, P);
            // return m_Water->cpmass() * q2;
            return props.Cp * q2;
        }

        double cNaCl::Cp_Solid(const double& T, const double& P)
        {
            double T_C = T - 273.15;
            double T_minus_Ttriple = T - T_Triple;
            m_coeff_H.r[3] = m_coeff_H.r3[0] + m_coeff_H.r3[1]*T_C + m_coeff_H.r3[2]*T_C*T_C;

            return m_coeff_H.r[0] + 2.0*m_coeff_H.r[1]*T_minus_Ttriple + 3.0*m_coeff_H.r[2]*T_minus_Ttriple*T_minus_Ttriple + m_coeff_H.r[3]*P + m_coeff_H.r[4]*P*P;
        }

        double cNaCl::DeltaH_fus(const double& T, const double& P)
        {
            return slope_Clapeyron * T * (1.0/Rho_Liquid(T, P) - 1.0/Rho_Solid(T, P));
        }

        double cNaCl::_H_Solid(const double& T, const double& P)
        {
            double T_C = T - 273.15;
            // double H1_TP = _H1(T_C, P);
            // double H1_TtripleP = _H1(m_H_h_ref.T_C, P);
            double H1 = m_H_h_ref.R1*T_C + m_H_h_ref.R2*T_C*T_C + m_H_h_ref.R3*pow(T_C,3) - m_H_h_ref.R0;
            double H2_TP = _H2(T_C, P);
            double H2_TPtriple = _H2(T_C, m_H_h_ref.P);

            // return H1_TP - H1_TtripleP + (H2_TP - H2_TPtriple);
            return H1 + (H2_TP - H2_TPtriple);
        }
        double cNaCl::_H1(const double& T_C, const double& P)
        {
            // R1*T + R2*T^2 + R3*T^3
            // R1 = R1[0] + R1[1]*P + R1[2]*P^2
             double R1[3]={m_coeff_H.r[0] - 2.0*m_coeff_H.r[1]*T_Triple_C + 3.0*m_coeff_H.r[2]*T_Triple_sqr_C,  m_coeff_H.r3[0], m_coeff_H.r[4]};
            return (m_coeff_H.r[0] - 2.0*m_coeff_H.r[1]*T_Triple_C + 3.0*m_coeff_H.r[2]*T_Triple_sqr_C + m_coeff_H.r3[0]*P + m_coeff_H.r[4]*P*P)*T_C 
                + (m_coeff_H.r[1] - 3.0*m_coeff_H.r[2]*T_Triple_C + 0.5*m_coeff_H.r3[1]*P)*T_C*T_C 
                + (m_coeff_H.r[2] + 1.0/3.0*m_coeff_H.r3[2]*P)*pow(T_C, 3.0);
        }

        void cNaCl::_calc_R_H1(H_halite_ref& coeff_haliteH_ref)
        {
            double R1[3], R2[2], R3[2];
            R1[0] = m_coeff_H.r[0] - 2.0*m_coeff_H.r[1]*T_Triple_C + 3.0*m_coeff_H.r[2]*T_Triple_sqr_C;
            R1[1] = m_coeff_H.r3[0];
            R1[2] = m_coeff_H.r[4];
            R2[0] = m_coeff_H.r[1] - 3.0*m_coeff_H.r[2]*T_Triple_C;
            R2[1] = 0.5*m_coeff_H.r3[1];
            R3[0] = m_coeff_H.r[2];
            R3[1] = 1.0/3.0*m_coeff_H.r3[2];
            // calculate R0 = R1*T0 + R2*T0^2 + R3*T0^3
            coeff_haliteH_ref.R1 = R1[0] + R1[1]*coeff_haliteH_ref.P + R1[2]*coeff_haliteH_ref.P*coeff_haliteH_ref.P;
            coeff_haliteH_ref.R2 = R2[0] + R2[1]*coeff_haliteH_ref.P;
            coeff_haliteH_ref.R3 = R3[0] + R3[1]*coeff_haliteH_ref.P;
            coeff_haliteH_ref.R0 = coeff_haliteH_ref.R1*coeff_haliteH_ref.T_C + coeff_haliteH_ref.R2*coeff_haliteH_ref.T_C*coeff_haliteH_ref.T_C + coeff_haliteH_ref.R3*pow(coeff_haliteH_ref.T_C, 3.0);
        }
        void cNaCl::_printCoeffs_haliteH()
        {
            // print R
            printf("-- Coefficients of R for halite enthalpy\n");
            printf("R0: %.8E\n", m_H_h_ref.R0);
            printf("R1: %.8E\n", m_H_h_ref.R1);
            printf("R2: %.8E\n", m_H_h_ref.R2);
            printf("R3: %.8E\n", m_H_h_ref.R3);
        }
        double cNaCl::_H2(const double& T_C, const double& P)
        {
            double exp_TC_l5 = exp(T_C/m_coeff_rho.l[5]);
            double l = m_coeff_rho.l[3] + m_coeff_rho.l[4]*exp_TC_l5;
            double rho_h = m_coeff_rho.l[0] + m_coeff_rho.l[1]*T_C +m_coeff_rho.l[2]*T_C*T_C + l*P;
            double H2_1  = log(rho_h)/l;

            double dldT = m_coeff_rho.l[4]*exp_TC_l5/m_coeff_rho.l[5];
            double H2_2  = ((m_coeff_rho.l[1] + 2.0*m_coeff_rho.l[2]*T_C + P*dldT)/rho_h - H2_1*dldT)/l;

            return H2_1 - (T_C+273.15)*H2_2;
        }
        double cNaCl::H_Solid(const double& T, const double& P)
        {
            return _H_Solid(T, P) + m_H_h_ref.H;
        }

        void cNaCl::Boiling_p(const std::vector<double> T, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = Boiling_p(T[i]);
        }

        void cNaCl::Melting_T(const std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = Melting_T(P[i]);
        }

        void cNaCl::Sublimation_p(const std::vector<double> T, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = Sublimation_p(T[i]);
        }

        void cNaCl::Melting_p(const std::vector<double> T, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = Melting_p(T[i]);
        }

        void cNaCl::Boiling_T(const std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = Boiling_T(P[i]);
        }

        void cNaCl::Sublimation_T(const std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = Sublimation_T(P[i]);
        }

        void cNaCl::T_Vapor(const std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = T_Vapor(P[i]);
        }

        void cNaCl::P_Vapor(const std::vector<double> T, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = P_Vapor(T[i]);
        }

        void cNaCl::Cp_Solid(const std::vector<double> T, const std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = Cp_Solid(T[i], P[i]);
        }

        void cNaCl::DeltaH_fus(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res)
        {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = DeltaH_fus(T[i], P[i]);
        }
        void cNaCl::H_Solid(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res)
        {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = H_Solid(T[i], P[i]);
        }
        void cNaCl::H_Liquid(const std::vector<double>& T, const std::vector<double>& P, std::vector<double>& res)
        {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = H_Liquid(T[i], P[i]);
        }

        std::string cNaCl::name_backend() {
            return m_backendname;
        }
    };

};
```


