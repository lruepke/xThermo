/**
 * @file IAPWS95.h 
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of IAPWS-95 EOS.
 * @version 0.1
 * @date 2022-04-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef IAPWS95_xThermal_H
#define IAPWS95_xThermal_H

#include "thermo.h"
#include "IAPWS-IF97.h"

namespace xThermal
{
    #define Name_Backend_IAPWS95 "IAPWS95"
    namespace IAPWS95
    {
        /**
         * @defgroup BITMASK_phi  Bit mask for updating of derivatives of phi
         *
         * @{
         */
        /** \f$ 2^1 \f$: \f$ \left(\frac{\partial \phi}{\partial \delta} \right)_{\tau} \f$ */
        #define Update_phi_d    2
        /** \f$2^2\f$: \f$ \left(\frac{\partial^2 \phi}{\partial \delta^2} \right)_{\tau} \f$ */
        #define Update_phi_dd   4
        /** \f$2^3\f$: \f$ \left(\frac{\partial \phi}{\partial \tau} \right)_{\delta} \f$ */
        #define Update_phi_t    8
        /** \f$2^4\f$: \f$ \left(\frac{\partial^2 \phi}{\partial \tau^2} \right)_{\delta} \f$ */
        #define Update_phi_tt   16
        /** \f$2^5\f$: \f$ \left(\frac{\partial^2 \phi}{\partial \delta \partial \tau} \right) \f$ */
        #define Update_phi_dt   32
        /** Update all */
        #define Update_phi_all  Update_phi_d|Update_phi_dd|Update_phi_t|Update_phi_tt|Update_phi_dt
        /** @} */

        /**
         * @defgroup PARAMS_NONLINEAR Parameters definition for nonlinear equation solving
         *
         * @{
         */
        /** Maximum iteration step */
        #define ITERATION_MAX   1000
        /** Tolerance for pressure with unit Pa */
        #define TOL_Pressure   1E-4
        /** Tolerance for solving saturated \f$ T, P, \rho_l, \rho_v \f$ at phase boundary, or for solving \f$\rho\f$ by given [T,P], or solving T by given [P, H] */
        #define TOL_PTRro 1E-10
        /** @} */

        /**
         * @brief Dimensionless Helmholtz free energy \f$ \phi = f/(RT) \f$ and its partial derivatives.
         *
         */
        struct HelmholtzEnergy_dimensionless
        {
            double value;    /**< \f$ \phi\f$ */
            double d;        /**< \f$ \left[ \frac{\partial \phi}{\partial \delta} \right]_{\tau} \f$ */
            double dd;       /**< \f$ \left[ \frac{\partial^2 \phi}{\partial \delta^2} \right]_{\tau} \f$ */
            double t;        /**< \f$ \left[ \frac{\partial \phi}{\partial \tau} \right]_{\delta} \f$ */
            double tt;       /**< \f$ \left[ \frac{\partial^2 \phi}{\partial \tau^2} \right]_{\delta} \f$ */
            double dt;       /**< \f$ \frac{\partial^2 \phi}{\partial \delta \partial \tau} \f$ */
            friend std::ostream& operator<<(std::ostream& os, const HelmholtzEnergy_dimensionless& phi)
            {
                return os << "v: "<<phi.value<<", d: "<<phi.d<<", dd: "<<phi.dd<<", t: "<<phi.t<<", tt: "<<phi.tt<<", dt: "<<phi.dt;
            }
            bool operator== (const HelmholtzEnergy_dimensionless &phi) const
            {
                return (fabs(value-phi.value) + fabs(d-phi.d) + fabs(dd-phi.dd) + fabs(t-phi.t) + fabs(tt-phi.tt) + fabs(dt-phi.dt))/6.0 < 1E-7;
            }
            HelmholtzEnergy_dimensionless operator+ (const HelmholtzEnergy_dimensionless &phi) const
            {
                HelmholtzEnergy_dimensionless tmp = {value + phi.value, d + phi.d, dd+phi.dd, t+phi.t, tt+phi.tt, dt+phi.dt};
                return tmp;
            }
            HelmholtzEnergy_dimensionless operator- (const HelmholtzEnergy_dimensionless &phi) const
            {
                HelmholtzEnergy_dimensionless tmp = {value - phi.value, d - phi.d, dd-phi.dd, t-phi.t, tt-phi.tt, dt-phi.dt};
                return tmp;
            }
        };
        struct Coeff_phi_o
        {
            const int n2 = 5;
            const double n0_term1[3] = {-8.3204464837497, 6.6832105275932, 3.00632};
            const double n0_term2[5] = {0.012436, 0.97315, 1.2795, 0.96956, 0.24873};
            const double gamma0_term2[5] = {1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105};
            // reused variables: optimization purpose
//            double one_minus_expGammaTau[5];
        };
        /**
         * @brief Numerical values of the coefficients and parameters of the residual part of the dimensionless Helmholtz free energy. Eq.6 in \cite IAPWS-95 cIAPWS95::phi_r
         *
         */
        struct Coeff_phi_r
        {
            const int n1 = 7, n2 = 44, n3=3, n4=2;
            const double ni_term1[7]   ={0.12533547935523e-1, 0.78957634722828e1,-0.87803203303561e1,
                                         0.31802509345418,   -0.26145533859358,  -0.78199751687981e-2,
                                         0.88089493102134e-2};
            const double di_term1[7]   ={1, 1, 1, 2, 2, 3, 4};
            const double ti_term1[7]   ={-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1};
            const double ni_term2[44]  = {-0.66856572307965,    0.20433810950965,       -0.66212605039687e-4,
                                          -0.19232721156002,   -0.25709043003438,        0.16074868486251,
                                          -0.40092828925807e-1, 0.39343422603254e-6,    -0.75941377088144e-5,
                                          0.56250979351888e-3,-0.15608652257135e-4,     0.11537996422951e-8,
                                          0.36582165144204e-6,-0.13251180074668e-11,   -0.62639586912454e-9,
                                          -0.10793600908932,    0.17611491008752e-1,     0.22132295167546,
                                          -0.40247669763528,    0.58083399985759,        0.49969146990806e-2,
                                          -0.31358700712549e-1,-0.74315929710341,        0.47807329915480,
                                          0.20527940895948e-1,-0.13636435110343,        0.14180634400617e-1,
                                          0.83326504880713e-2,-0.29052336009585e-1,     0.38615085574206e-1,
                                          -0.20393486513704e-1,-0.16554050063734e-2,     0.19955571979541e-2,
                                          0.15870308324157e-3,-0.16388568342530e-4,     0.43613615723811e-1,
                                          0.34994005463765e-1,-0.76788197844621e-1,     0.22446277332006e-1,
                                          -0.62689710414685e-4,-0.55711118565645e-9,    -0.19905718354408,
                                          0.31777497330738,   -0.11841182425981};
            const double ci_term2[44]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
                                          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,
                                          6, 6};
            const double di_term2[44]  = {1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
                                          4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14,
                                          3, 6, 6, 6};
            const double ti_term2[44]  = {4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10,
                                          10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,
                                          23, 10, 50, 44, 46, 50};
            const double ni_term3[3]       = {-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4};
            const double d_term3          = 3;
            const double ti_term3[3]       = {0, 1, 4};
            const double alpha_term3      = 20;
            const double betai_term3[3]    = {150, 150, 250};
            const double gammai_term3[3]   = {1.21, 1.21, 1.25};
            const double epsilon_term3    = 1;
            const double ni_term4[2]       = {-0.14874640856724, 0.31806110878444};
            const double a_term4          = 3.5;
            const double bi_term4[2]       = {0.85, 0.95};
            const double B_term4          = 0.2;
            const double Ci_term4[2]       = {28, 32};
            const double Di_term4[2]       = {700, 800};
            const double A_term4          = 0.32;
            const double beta_term4       = 0.3;
        };

        /**
    * @brief Constants for water viscosity calculation, see \cite mu2008.
    */
        struct Constants_Viscosity2008_Water
        {
            const double T_star = 647.096;
            const double rho_star = 322.0;
            const double p_star = 22.064E5;
            const double mu_star = 1.0E-6;
            static const int N_Hi = 4;
            const double H_i[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
            static const int row_Hij = 6;
            static const int col_Hij = 7;
            const double H_ij[6][7] = {
                    {5.20094e-1,    2.22531e-1,  -2.81378e-1,   1.61913e-1,  -3.25372e-2,  0,            0},
                    {8.50895e-2,    9.99115e-1,  -9.06851e-1,   2.57399e-1,   0,           0,            0},
                    {-1.08374,      1.88797,     -7.72479e-1,   0,            0,           0,            0},
                    {-2.89555e-1,   1.26613,     -4.89837e-1,   0,            6.98452e-2,  0,           -4.35673e-3},
                    {0,             0,           -2.57040e-1,   0,            0,           8.72102e-3,   0},
                    {0,             1.20573e-1,   0,            0,            0,           0,           -5.93264e-4}
            };
            // critical-region constants
            const double x_mu = 0.068;
            const double qC = 1.0/1.9;
            const double qD = 1.0/1.1;
            const double nu = 0.630;
            const double gamma = 1.239;
            const double xi0 = 0.13;
            const double Gamma0 = 0.06;
            const double T_R_bar = 1.5;
            const double T_R = T_R_bar * T_star;
            const double p_starByRho_star = p_star/rho_star;
            const double nuBygamma = nu/gamma;
        };

        struct Params_T_Sat_estimate;
        struct Params_SolvePhaseEquilibrium;
        struct Params_TP2Rho;
        struct Params_HP2RhoT;
        struct HelmholtzEnergy_dimensionless_SinglePhase
        {
            HelmholtzEnergy_dimensionless phio, phir;
        };
        struct HelmholtzEnergy_dimensionless_TwoPhase
        {
            HelmholtzEnergy_dimensionless phio_l, phir_l, phio_v, phir_v;
        };
        struct STRUCT_delta_TwoPhase
        {
            double delta_l, delta_v;
        };
        struct State
        {
            PhaseRegion phase;
            double x     = 0;      /**< vapor mass fraction */
            double tau;
            union
            {
                HelmholtzEnergy_dimensionless_SinglePhase   singlePhase;
                HelmholtzEnergy_dimensionless_TwoPhase      twoPhase;
            }phi;
            union
            {
                double    singlePhase;
                STRUCT_delta_TwoPhase       twoPhase;
            }delta;
        };

        typedef std::map<unsigned long int, ThermodynamicProperty> MAP_INDEX2PROP;

        /**
         * @brief Class of IAPWS-95 formula of H2O, which inherits from IAPWS_IF97::cIAPWS_IF97.
         *
         */
        class cIAPWS95 : public cxThermal
        {
        private:
            Coeff_phi_o m_coeff_phio;
            Coeff_phi_r m_coeff_phir;
            IAPWS_IF97::cIAPWS_IF97 m_IF97;
            Constants_Viscosity2008_Water m_constants_mu2008;
        public:
            cIAPWS95(/* args */);
            ~cIAPWS95();
        private:
            void initialize_data();
            // void print_state_PhaseEquilibrium(size_t iter, gsl_multiroot_fsolver * s);
            // void print_state_TP2Rho(size_t iter, gsl_multiroot_fsolver * s);
            // void print_state_HP2RhoT(size_t iter, gsl_multiroot_fsolver * s);
            double Rho_Newton(const double T_K, const double P);
            double Rho_bisection(const double T_K, const double P, double Rho_guess = 322.0, double Rho_min = 1E-4, double Rho_max = 1400);
            void _enthalpy(const double& T_K, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir, double& h);
            void SinglePhase_HP(const double& H, const double& P, double& rho, double& T_K, std::string method = "bisection");
            // update state
            void UpdateState_HP(ThermodynamicProperties& props, State& state, const double& H, const double& P, std::string method = "bisection");
            void UpdateState_TP(ThermodynamicProperties& props, State& state, const double& T, const double& P);
            // Thermodynamic properties functions
            void _h(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir);
            void h(ThermodynamicProperties& props,const State& state);
            void _dhdT_P(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir);
            void dhdT_P(ThermodynamicProperties& props,const State& state);
            void _dPdRho_T(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir) const;
            void dPdRho_T(ThermodynamicProperties& props,const State& state);
            void _dPdT_Rho(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir) const;
            void dPdT_Rho(ThermodynamicProperties& props,const State& state);
            void Viscosity_H2O_IAPWS2008(const double& T, const double& Rho, double& Mu);
        public:
            CONSTENTS_Thermo m_constants;
            void phi_o(const double& delta, const double& tau, HelmholtzEnergy_dimensionless& phio, unsigned int update_derivatives=Update_phi_all);
            double phi_r(const double& delta, const double& tau);
            double phi_r_d(const double& delta, const double& tau);
            double phi_r_dd(const double& delta, const double& tau);
            double phi_r_t(const double& delta, const double& tau);
            double phi_r_tt(const double& delta, const double& tau);
            double phi_r_dt(const double& delta, const double& tau);
            void phi_r(const double& delta, const double& tau, HelmholtzEnergy_dimensionless& phir);
            double P_Sat_estimate(const double& T_K);
            void P_Sat_estimate(const double& T_K, double& p, double& dpdT);
            double T_Sat_estimate(const double& P);
            double Rho_Liquid_Sat_estimate(const double& T);
            double Rho_Vapor_Sat_estimate(const double& T);
            void Boiling_p(const double& T_K, double& P, double& rho_l, double& rho_v);
            void Boiling_T(const double& P, double& T_K, double& rho_l, double& rho_v);
            double P_delta_tau(const double& delta, const double& tau);
            double Rho(const double& T_K, const double& P, std::string method = "bisection");
            double _enthalpy(const double& T_K, const double& P, std::string method = "bisection");
            double Mu(const double& T, const double& P);
            void Verification_Mu();
        public:
            std::string name(){return Name_Backend_IAPWS95;};
            double Tmin(){return m_constants.Tmin; };                  /**< Get the minimum temperature in K */
            double Tmax(){return m_constants.Tmax;};                  /**< Get the maximum temperature in K */
            double pmin(){return m_constants.pmin; };                    /**< Get the minimum pressure in Pa */
            double pmax(){return m_constants.pmax;};                  /**< Get the maximum pressure in Pa */
            double Ttriple(){return m_constants.Ttriple;};               /**< Get the triple point temperature in K */
            double T_critical(){return m_constants.T_critical;};            /**< Return the critical temperature in K */
            double p_critical(){return m_constants.p_critical;};            /**< Return the critical pressure in Pa */
            double rhomass_critical(){return m_constants.rhomass_critical;};      /**< Return the critical mass density in \f$ kg/m^3 \f$ */
            double molar_mass(){return m_constants.molar_mass;};       /**< Return the molar mass in kg/mol */

            // update state and properties for give T,p
            PhaseRegion findPhaseRegion_TPX(const double& T, const double& p, const double& X=0);
            void UpdateState_TPX(ThermodynamicProperties& props,const double& T, const double& P, const double& X=0);
            // void UpdateState_TPX(ThermodynamicPropertiesArray& stateArray, const size_t& N, const double* T, const double* p, const double* X=NULL); //for Matlab API
            ThermodynamicProperties UpdateState_TPX(const double &T, const double &p, const double &X=0); //for Python API
            PhaseRegion findPhaseRegion_HPX(const double& H, const double& p, const double& X=0);
            void UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X=0);
            double Boiling_p(const double& T);
            double Boiling_p(const double& T, double& rho_l, double& rho_v);
            double Boiling_p(const double& T, ThermodynamicProperties& props);
            double Boiling_T(const double& p);
            double Boiling_T(const double& p, double& rho_l, double& rho_v);
            double Boiling_T(const double& p, ThermodynamicProperties& props);
        public:
            // vector version for python API
//            void P_Sat_estimate(std::vector<double> T_K, std::vector<double>& res);
//            void Rho_Liquid_Sat_estimate(std::vector<double> T_K, std::vector<double>& res);
//            void Rho_Vapor_Sat_estimate(std::vector<double> T_K, std::vector<double>& res);
//            void Boiling_P(const std::vector<double> T_K, std::vector<double>& P, std::vector<double>& rho_l, std::vector<double>& rho_v);
//            void Boiling_T(const std::vector<double> P, std::vector<double>& T_K, std::vector<double>& rho_l, std::vector<double>& rho_v);
//            void P_delta_tau(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r_d(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r_dd(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r_t(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r_tt(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
            void phi_r_dt(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res);
//            void Rho(const std::vector<double> T_K, const std::vector<double> P, std::vector<double>& res, std::string method = "bisection");
//            void _enthalpy(const std::vector<double> T_K, const std::vector<double> P, std::vector<double>& res, std::string method = "bisection");
//            void getProperties_HP(std::vector<double> H, std::vector<double> P, unsigned long int Request_properties, std::vector<double>& properties, std::vector<double>& properties_l, std::vector<double>& properties_v, std::vector<unsigned long int>& propIndexReturn);
//            void getProperties_TP(std::vector<double> T, std::vector<double> P, unsigned long int Request_properties, std::vector<double>& properties, std::vector<double>& properties_l, std::vector<double>& properties_v, std::vector<unsigned long int>& propIndexReturn);
//            ThermodynamicProperty getPropInfo(unsigned long int which);
        };


        /**
         * @brief Data struct for phase equilibrium solving. It is used in \link func_PhaseEquilibrium \endlink , \link cIAPWS95::Boiling_P \endlink and \link cIAPWS95::Boiling_T \endlink
         *
         */
        struct Params_SolvePhaseEquilibrium
        {
            cIAPWS95* iapws;
            double RT; /**< xThermal::R times T: \f$RT \f$ [J/kg] */
            /**
             * @brief One of \f$ \tau \f$ and \f$p\f$, which is used for solving saturated pressure and saturated temperature, respectively.
             * \note Because for solving phase equilibrium, only P or T is given, so here I use an union to store it.
             */
            union
            {
                double tau; /**< \f$T_c/T \f$ */
                double P;  /**< Pressure [Pa] */
            }TorP;
            SOLVE_SATURATED_PorT Solve_PorT; /**< Flag of solving \f$ P_{sat} \f$ or \f$ T_{sat} \f$, the value will be one of #SOLVE_SATURATED_P, #SOLVE_SATURATED_T */
        };
        struct Params_TP2Rho
        {
            cIAPWS95* iapws;
            double T_K, tau, p, RhocRT;
        };
        struct Params_HP2RhoT
        {
            cIAPWS95* iapws;
            double h, p;
        };
        struct Params_T_Sat_estimate
        {
            cIAPWS95* iapws;
            double p;
        };

    };

};

#endif