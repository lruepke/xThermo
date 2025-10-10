/**
 * @file H2ONaCl.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of H2O-NaCl model.
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "H2ONaCl.h"
// ----- LUT related head filess --------
#include "LookUpTableForestI.H"
#include "interpolationI.H"
#include "AMR_LUT_RefineFuncI.H"
//---------------------------------------

namespace xThermal
{
    namespace H2ONaCl
    {
        using namespace std;
        /**
        * @brief Initialize all necessary data (member variables), include all constants and coefficients.
        */
        void cH2ONaCl::initialize_data()
        {
            createTable4_Driesner2007a(m_tab4);
            // initialize thermodynamic constants of Water
            m_constants_Water.Tmin = m_Water->Tmin();
            m_constants_Water.Tmax = m_Water->Tmax();
            m_constants_Water.pmin = m_Water->pmin();
            m_constants_Water.pmax = m_Water->pmax();
            m_constants_Water.Ttriple = m_Water->Ttriple();
            m_constants_Water.T_critical = m_Water->T_critical();
            m_constants_Water.p_critical = m_Water->p_critical();
            m_constants_Water.rhomass_critical = m_Water->rhomass_critical();
            m_constants_Water.molar_mass = m_Water->molar_mass();
            m_constants_Water.rhomolar_critical = m_Water->rhomolar_critical();
        }
        cH2ONaCl::cH2ONaCl(std::string backend_H2O)
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
                throw NotImplementedError("Error in cH2ONaCl(std::string backend_H2O). The supported H2O backend is only one of 'IAPS84', 'IAPWS95' and 'IAPWS95_CoolProp'. Input name is "+backend_H2O);
            }
            m_NaCl = new NaCl::cNaCl(m_backendname);
            // initialization
            initialize_data();
        }
        // copy constructor
        cH2ONaCl::cH2ONaCl(const cH2ONaCl& sw)
        {
            m_backendname = sw.m_backendname;
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
            else
            {
                std::string available_backends = "'IAPS84', 'IAPWS95'";
#ifdef USE_COOLPROP
                available_backends = "'IAPS84', 'IAPWS95' and 'IAPWS95_CoolProp'";
#endif
                throw NotImplementedError("Error in cH2ONaCl(const cH2ONaCl& sw). The supported H2O backend is only one of " + available_backends + ". Input name is "+sw.m_backendname);
            }
            m_NaCl = new NaCl::cNaCl(m_backendname);
            // initialization
            initialize_data();
        }

        std::string cH2ONaCl::name_backend()
        {
            return m_backendname;
        }

        cH2ONaCl::~cH2ONaCl()
        {
            if(m_NaCl)delete m_NaCl;
            if(m_Water)delete m_Water;
        }

        /**
         * @brief Parameters for the critical curve Critical, see table 4 of \cite Driesner2007Part1
         * 
         * \note Change 500 deg.C to 773.15 K to make temperature with unit [K]. Times \f$ 10^5 \f$ to \f$ c_{i}\f$ to make critical pressure with unit [Pa]
         * 
         * \warning The critical pressure (\f$ 2.2054915\times 10^2 \f$ bar) and temperature(373.976\f$ ^{\circ}\f$C) of H2O used in Driesner & Heinrich(2007) formula comes from IAPS-84 EOS, which is different from that of IAPWS95_CoolProp (m_constants_Water.T_critical, H2O::P_c). The difference between them are \f$ \Delta T=0.0299 ^{\circ} \f$C, \f$ \Delta P = 9085 \f$ Pa.
         * 
         * @param table4 
         */
        void cH2ONaCl::createTable4_Driesner2007a(Table4_Driesner2007a_CriticalCurve& table4)
        {
            // Table 4 of Driesner and Heinrich(2007)
            double c[14] = {-2.36, 0.128534, -0.023707, 0.00320089, -0.000138917, 
                            1.02789E-07, -4.8376E-11, 2.36, -0.0131417, 0.00298491,
                            -0.000130114, 581.0101000000000, 3.405487809994650, -0.000488336};// c[11] and c[12] are calculated below
            double cA[11] = {1, 1.5, 2, 2.5, 3, 4, 5, 1, 2, 2.5, 3};
            // =============== Change c*1E5 to make pressure with unit Pa ==============
            for(int i=0;i<14;i++)c[i]=c[i]*1E5; 
            // =========================================================================
            // Interesting!!! If calculate c11,c12 very actually, and the result will have up to 0.0003 difference for VL of V+L surface at (1.000000e+03 deg.C,2.161670e+03 bar) see mmc4 of Driesner(2007b)
            // the pre-calculate value above is from Falko's matlab code: calc_boundary_prop.m, this result close to both Falko's and mmc4 table.
//            // c[11] (c12 in Driesner and Heinrich(2007)) is the value of P_crit at 500 deg.C, calculated from eq. 5b.
//            // c[12] (c13) is the first temperature derivative of eq. 5b at 500 deg.C
//            for (size_t i = 7; i < 11; i++)
//            {
//                c[11] += c[i] * pow(773.15 - m_constants_Water.T_critical, cA[i]); //the second term of eq. 5b in Drisner and Heinrich (2007)
//                c[12] += c[i] * cA[i] * pow(773.15 - m_constants_Water.T_critical, cA[i] - 1); //the first temperature derivative of eq. 5b
//            }
//            c[11] = m_constants_Water.p_critical + c[11];
            double d[11] = {8E-05, 1E-05, -1.37125E-07, 9.46822E-10, -3.50549E-12, 6.57369E-15, 
                            -4.89423E-18, 7.77761E-2, 2.7042E-4, -4.244821E-07, 2.580872E-10};
            // copy 
            for(int i=0;i<14; i++)table4.c[i]    = c[i];
            for(int i=0;i<11; i++)table4.cA[i]   = cA[i];
            for(int i=0;i<11; i++)table4.d[i]    = d[i];
        }

        /**
         * @brief Convert mass fraction of NaCl to molar fraction. 
         * 
         * \f{align}
         * X_{mol} = \frac{X_{wt}/M_{NaCl}}{X_{wt}/M_{NaCl} + (1-X_{wt})/M_{H_2O}}
         * \f}
         * 
         * @param X_wt [0,1]
         * @return double [0,1]
         */
        double cH2ONaCl::Wt2Mol(const double& X_wt)
        {
            double X_mol;
            Wt2Mol(X_wt, X_mol);
            return X_mol;
        };
        void cH2ONaCl::Wt2Mol(const double& X_wt, double& X_mol)
        {
            X_mol = X_wt/m_NaCl->molar_mass()/(X_wt/m_NaCl->molar_mass() + (1-X_wt)/m_constants_Water.molar_mass);
        }
        std::vector<double> cH2ONaCl::Wt2Mol(std::vector<double> X_wt)
        {
            std::vector<double> x_mol;
            for (size_t i = 0; i < X_wt.size(); i++)
            {
                x_mol.push_back(Wt2Mol(X_wt[i]));
            }
            return x_mol;
        }

        /**
         * @brief Convert molar fraction of NaCl to mass fraction.
         * 
         * \f{align}
         * X_{wt} = \frac{X_{mol}M_{NaCl}}{X_{mol}M_{NaCl} + (1-X_{mol})M_{H_2O}}
         * \f}
         * 
         * @param X_mol [0,1]
         * @return double [0,1]
         */
        double cH2ONaCl::Mol2Wt(double X_mol)
        {
            return m_NaCl->molar_mass() * X_mol / (m_NaCl->molar_mass() * X_mol + (1 - X_mol) * m_constants_Water.molar_mass);
        };

        std::vector<double> cH2ONaCl::Mol2Wt(std::vector<double> X_mol)
        {
            std::vector<double> x_wt;
            for (size_t i = 0; i < X_mol.size(); i++)
            {
                x_wt.push_back(Mol2Wt(X_mol[i]));
            }
            return x_wt;
        };

        /**
         * @brief Critical composition as a function of temperature, see Eq. 7a,b of \cite Driesner2007Part1.
         * 
         * \warning Equation 7b and 7a gives different values at (T=600 \f$ ^{\circ} C\f$). When T=600deg.C, Eq.7a gives \f$X_{crit}\f$ = 7.778263e-02 (mol fraction) and Eq.7b gives \f$X_{crit}\f$ = 7.777610e-02 (mol fraction) which is the value in supplement file mmc2 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param X_crit [mol/mol] 
         */
        void cH2ONaCl::X_Critical_mol(double T, double& X_crit)
        {
            X_crit = 0;
            if(T < m_constants_Water.T_critical && T>=Tmin())
            {
                X_crit = 0;
            }else if (T>=m_constants_Water.T_critical && T<873.15){
                for (size_t i = 0; i < 7; i++){
                    X_crit += m_tab4.d[i]*pow(T - m_constants_Water.T_critical, i+1); //eq. 7a
                }
            }else if(T >= 873.15 && T <= Tmax()){
                for (size_t i = 7; i < 11; i++){
                    X_crit += m_tab4.d[i]*pow(T - 873.15, i-7); //eq. 7b
                }
            }else
            {
                if(!isnan(T))std::cout<<WARN_COUT<<"T: "<<T<<" K out of temperature range: ["<<Tmin()<<", "<<Tmax()<<"], in void X_Critical_mol(double T, double& X_crit)"<<std::endl;
            }
        }

        /**
         * @brief Critical composition as a function of temperature, see Eq. 7a,b of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param X_crit [kg/kg] 
         */
        void cH2ONaCl::X_Critical(double T, double& X_crit)
        {
            X_Critical_mol(T, X_crit);
            // convert mole fraction to mass fraction
            X_crit = Mol2Wt(X_crit);
        }

        /**
         * @brief Critical pressure as a function of temperature, see Eq. 5a-c of \cite Driesner2007Part1.
         * \note Change 500 deg.C to 773.15 K to make temperature with unit [K].
         * 
         * \warning The critical pressure (\f$ 2.2054915\times 10^2 \f$ bar) and temperature(373.976\f$ ^{\circ}\f$C = 647.126K) of H2O used in Driesner & Heinrich(2007) formula comes from IAPS-84 EOS(see PROST code \cite prost), which is different from that of IAPWS95_CoolProp (m_constants_Water.T_critical, H2O::P_c). The difference between them are \f$ \Delta T=0.0299 ^{\circ} \f$C, \f$ \Delta P = 9085 \f$ Pa.
         * 
         * @param T [K]
         * @param P_crit [Pa] 
         */
        void cH2ONaCl::P_Critical(double T, double& P_crit)
        {
            P_crit=0;
            if(T < m_constants_Water.T_critical && T>=Tmin()){
                for (size_t i = 0; i < 7; i++){
                    P_crit += m_tab4.c[i]*pow(m_constants_Water.T_critical - T, m_tab4.cA[i]); //eq. 5a. Note: this result is a little bit different with IAPWS95_CoolProp boiling curve.
                }
                P_crit+=m_constants_Water.p_critical;
            }else if(T >= m_constants_Water.T_critical && T <= 773.15){
                for (size_t i = 7; i < 11; i++){
                    P_crit += m_tab4.c[i]*pow(T - m_constants_Water.T_critical, m_tab4.cA[i]); //eq. 5b
                }
                P_crit+=m_constants_Water.p_critical;
            }else if(T > 773.15 && T <= Tmax()){
                for (size_t i = 11; i < 14; i++){
                    P_crit += m_tab4.c[i]*pow(T - 773.15, i-11); //eq. 5c
                }
            }else
            {
                if(!isnan(T))std::cout<<WARN_COUT<<"T: "<<T<<" K out of temperature range: ["<<Tmin()<<", "<<Tmax()<<"], in void cH2ONaCl::P_Critical(double T, double& P_crit)"<<std::endl;
            }
        }

        /**
        * @brief Solve critical T and X using bisection method for given Temperature.
        * @param P [Pa]
        * @param T_crit [K]
        */
        void cH2ONaCl::T_Critical(double P, double& T_crit)
        {
            T_Critical_Bisection(P, T_crit);
        }

        /**
         * @brief Roughly estimate critical temperature from given pressure.
         * 
         * Use a poly fit function to explicitly express critical temperature as a function of pressure.
         * The forward calculation of critical pressure using \link P_Critical \endlink
         * 
         * \image html H2ONaCl/H2ONaCl_CriticalCurve_log10P_T.svg Estimate critical temperature by given pressure width=80%
         * 
         * @param P [Pa]
         * @param T_crit [K] 
         */
        void cH2ONaCl::T_Critical_estimate(double P, double& T_crit)
        {
            double y = log10(P);
            T_crit = m_coeff_estimate_Tcrit.coeffs[0];
            for (int i = 1; i < m_coeff_estimate_Tcrit.num; i++)
            {
                T_crit += m_coeff_estimate_Tcrit.coeffs[i]*pow(y, i);
            }
            // make sure the estimated value in the valid range.
            if(T_crit>=Tmax())T_crit=Tmax()-5;
            else if(T_crit<=Tmin())T_crit=Tmin()+5;
        }

        /**
         * @brief Construct function of \f$ P_{critical} = f(T) \f$, see Eq. 5a-c of \cite Driesner2007Part1.
         * 
         * @param x Temperature [K]
         * @param params 
         * @param f Pressure [Pa]
         * @return int 
         */
        double func_P2CriticalT(double T, void *params)
        {
            Params_P2CriticalT* param = (Params_P2CriticalT*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P; 
            double P_crit;
            sw->P_Critical(T, P_crit);
            return P_crit - p;
        }

        /**
         * @brief Calculate critical temperature by given pressure.
         * 
         * This is the inversion of \link P_Critical \endlink
         * 
         * @param P [Pa]
         * @param T_crit [K]
         */
        void cH2ONaCl::T_Critical_Bisection(double P, double& T_crit)
        {
            double x_lo = Tmin(), x_hi = Tmax();
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;
            
            gsl_function F;
            Params_P2CriticalT params={this, P};

            F.function = &func_P2CriticalT;
            F.params = &params;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                T_crit = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/1E4);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
                // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_crit);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            if(!((status == GSL_SUCCESS)))
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                printf("P = %.3E Pa\n", P);
                ERROR("Fatal error in void cH2ONaCl::T_Critical(double P, double& T_crit)");
            }
            gsl_root_fsolver_free (s);
        }

        double func_TL_VL(double T, void *params)
        {
            auto* param = (Params_PX2T_VL*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P;
            double X = param->X;
            double X_vl;
            X_vl = (param->phase == Vapor ? sw->XV_VL(T, p) : sw->XL_VL(T, p) );
            return X_vl - X;
        }
        /**
        * Calculate T (inversion) of V+L liquid branch for given X and P
        * @param P [Pa]
        * @param X [kg/kg]
        * @param T [K]
        */
        void cH2ONaCl::T_VL_bisection(double P, double X, double& T_res, PhaseType phase)
        {
            double x_lo=Tmin(), x_hi=Tmax();
            // double T_crit0, Xl_vlh_min, Xv_vlh_min;
            // if (P>P_Peak_VLH)
            // {
            //     T_Critical(P, x_lo);
            //     x_hi = Tmax();
            // }else if(P>=pmin() && P<=P_Peak_VLH)
            // {
            //     double T1, T2;
            //     T_VLH_P0(P, T1, T2);
            //     double Xl1_vlh, Xl2_vlh, Xv1_vlh, Xv2_vlh;
            //     X_VLH(T1, P, Xl1_vlh, Xv1_vlh);
            //     X_VLH(T2, P, Xl2_vlh, Xv2_vlh);
            //     if (X<=Xl1_vlh)
            //     {
            //         if(P<m_Water->p_critical())
            //         {
            //             double T_boil = m_Water->Boiling_T(P);
            //             double T1, T2;
            //             T_VLH_P0(P, T1, T2);
            //             x_lo = T_boil;
            //         } else
            //         {
            //             T_Critical(P, x_lo);
            //         }
            //         x_hi = T1;
            //     } else if (X>=Xl2_vlh)
            //     {
            //         x_lo = T2;
            //         x_hi = Tmax();
            //     } else
            //     {
            //         char errorinfo[200];
            //         sprintf(errorinfo,"The input parameters are wrong in TL_VL_bisection(double P, double X, double& T_res): P=%f bar, X = %.f wt%% NaCl", P/1E5, X*100.0);
            //         ERROR(string(errorinfo));
            //     }
            // } else
            // {
            //     char errorinfo[200];
            //     sprintf(errorinfo,"The input P in TL_VL_bisection(double P, double X, double& T_res) out of valid range [%f, %f]bar: P=%f bar, X = %.f wt%% NaCl", pmin()/1E5, pmax()/1E5, P/1E5, X*100.0);
            //     ERROR(string(errorinfo));
            // }
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;

            gsl_function F;
            Params_PX2T_VL params={this, P, X, phase};

            F.function = &func_TL_VL;
            F.params = &params;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo-273.15, x_hi-273.15,T_res-273.15);
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                T_res = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_PTRro);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
                // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo-273.15, x_hi-273.15,T_res-273.15);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            if(status != GSL_SUCCESS)
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                printf("P = %.3E Pa, X = %.3E wt%% NaCl\n", P, X*100.0);
                ERROR("Fatal error in void cH2ONaCl::T_VL_LiquidBranch_bisection(double P, double X, double& T_res)");
            }
            gsl_root_fsolver_free (s);
        }

        /**
        * Calculate temperature of liquid branch for given P and X. It is used for HPX calculation in the low salinity region and P<Pcrit_H2O
        * @param P [Pa]
        * @param X [kg/kg]
        * @return
        */
        double cH2ONaCl::T_VL(const double& P, const double& X, PhaseType phase)
        {
            double T_vl;
            T_VL_bisection(P,X, T_vl, phase);
            return T_vl;
        }


        /**
         * @brief Calculate critical pressure(equation 5a-c of \cite Driesner2007Part1) and salinity(equation 7a,b of \cite Driesner2007Part1) by given temperature.
         * 
         * - Benchmark test
         * 
         * If use the critical point of H2O in Driesner(2007a)\cite Driesner2007Part1, this function will get completely the same result with <a href="https://www.sciencedirect.com/science/article/pii/S0016703707002955?via%3Dihub#aep-e-component-id41">Electronic Annex EA-2</a> published with \cite Driesner2007Part1.
         * 
         * \image html H2ONaCl/H2ONaCl_CriticalCurve.svg (a) Critical curve and (b) difference with different critical point value of H2O width=40%
         * 
         * @param T [K]
         * @param P_crit [Pa]
         * @param X_crit [kg/kg]
         */
        void cH2ONaCl::P_X_Critical(double T, double& P_crit, double& X_crit)
        {
            // calculate critical pressure
            P_Critical(T, P_crit);
            // calculate critical salinity
            X_Critical(T, X_crit);
        } 
        
        /**
         * @brief Calculate critical temperature and composition by given pressure.
         * 
         * This is the inversion of \link P_X_Critical \endlink
         * 
         * @param P [Pa]
         * @param T_crit [K]
         * @param X_crit [kg/kg]
         */
        void cH2ONaCl::T_X_Critical(double P, double& T_crit, double& X_crit)
        {
            T_Critical(P, T_crit);
            X_Critical(T_crit, X_crit);
        }

        /**
         * @brief Halite-saturated liquid composition, see equation 8 of \cite Driesner2007Part1.
         * 
         * \note Change coefficient e to make pressure with unit [Pa]
         * 
         * \image html H2ONaCl/HaliteLiquidus3D.svg Liquid composition for the halite liquidus. (a) Full range temperature–pressure dependence, sets of symbols are for the same pressures as sets of lines; (b) pressure dependence at [25,50,100] deg.C ; (c) 3D view. width=80%
         * 
         * \see Fig.7a,b of \cite Driesner2007Part1. \link P_VLH \endlink
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [mol/mol] 
         */
        double cH2ONaCl::X_HaliteLiquidus_mol(const double& T, const double& P)
        {
            const double P_squre = P*P;
            // Table 5 of Driesner and Heinrich(2007)
            double e[6] = {0.0989944 + 3.30796E-11 * P - 4.71759E-20 * P_squre, 
                        0.00947257 - 8.6646E-11 * P + 1.69417E-19 * P_squre, 
                        0.610863 - 1.51716E-10 * P + 1.1929E-18 * P_squre, 
                        -1.64994 + 0.000203441E-5 * P - 6.46015E-18 * P_squre, 
                        3.36474 - 0.000154023E-5 * P + 8.17048E-18 * P_squre, 
                        1
                            };
            for (size_t i = 0; i < 5; i++)
            {
                e[5] -= e[i];
            }
            const double TbyT_hm = (T-273.15)/m_NaCl->Melting_T_C(P);
            double X_Liquids = 0;
            for (size_t i = 0; i < 6; i++)
            {
                X_Liquids += e[i] * pow(TbyT_hm, i);
            }
            if(X_Liquids>1)X_Liquids=1; //When T>NaCl::T_triple, the above equation calculate result will greater than 1, but the X_Liquidus should be 1. This is important for XV_VL calculation, i.e., Eq.15 log10_XL_P_NaCl value calculation.
            return X_Liquids;
        }

        /**
         * @brief Halite-saturated liquid composition, see equation 8 of \cite Driesner2007Part1.
         * 
         * \note Change coefficient e to make pressure with unit [Pa]
         * 
         * \image html H2ONaCl/HaliteLiquidus3D.svg Liquid composition for the halite liquidus. (a) Full range temperature–pressure dependence, sets of symbols are for the same pressures as sets of lines; (b) pressure dependence at [25,50,100] deg.C ; (c) 3D view. width=80%
         * 
         * \see Fig.7a,b of \cite Driesner2007Part1. \link P_VLH \endlink
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/kg] 
         */
        double cH2ONaCl::X_HaliteLiquidus(const double& T, const double& P)
        {
            double X_Liquids = Mol2Wt(X_HaliteLiquidus_mol(T,P));
            if(X_Liquids>1)X_Liquids=1; //make sure X_Liquids in range of [0,1]
            return X_Liquids;
        }

        /**
         * @brief Vapor-Liquid-Halite coexistence pressure as a function of temperature.
         * 
         * See Eq. 10 of \cite Driesner2007Part1.
         * 
         * \image html H2ONaCl/H2ONaCl_P_VLH.svg Pressure at vapor + liquid + halite coexistence. (a) Full range; (b) Liquid composition;  (c) low temperature, logarithmic pressure scale; (d) 3D view width=100%
         * 
         * \see Fig. 9,10 of \cite Driesner2007Part1. \link X_HaliteLiquidus \endlink and \link X_VLH \endlink
         * 
         * @param T [K]
         * @param P [Pa]
         */
        void cH2ONaCl::P_VLH(const double& T, double& P)
        {
            P = 0;
            double T_by_Ttriple_NaCl = (T - 273.15)/NaCl::T_Triple_C; //Temperatue of the original formula with unit [deg.C]
            for (int i = 0; i < m_tab6.num; i++)
            {
                P += m_tab6.f[i]*pow(T_by_Ttriple_NaCl, i);
            }
        }
        
        double cH2ONaCl::P_VLH(const double& T)
        {
            double P = 0;
            P_VLH(T,P); return P;
        }
        /**
        * Calculate derivative \f$ \frac{dP}{dT} \f$ along VLH saturation curve.
        * @param T
        * @return
        */
        double cH2ONaCl::dPdT_VLH(const double& T)
        {
            double dPdT = 0;
            double T_by_Ttriple_NaCl = (T - 273.15)/NaCl::T_Triple_C; //Temperatue of the original formula with unit [deg.C]
            for (int i = 1; i < m_tab6.num; i++)
            {
                dPdT += i*m_tab6.f[i]*pow(T_by_Ttriple_NaCl, i-1);
            }
            return dPdT/NaCl::T_Triple_C;
        }
        double func_T_VLH(double T, void *params)
        {
            Params_Inversion_PTX* param = (Params_Inversion_PTX*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P; 
            double P_vlh;
            sw->P_VLH(T, P_vlh);
            return P_vlh - p;
        }
        
        /**
         * @brief Calculate one temperature boundary (root) of VLH zone by given pressure.
         * 
         * @param P0 [Pa] valid range [#P_MIN, #P_Peak_VLH]
         * @param Tmin [K]
         * @param Tmax [K]
         * @param T_ [K]
         */
        double cH2ONaCl::T_VLH_P0_(double P0, double Tmin, double Tmax)
    {
        double T_ = 0;
        double x_lo = Tmin, x_hi = Tmax;
        int status;
        int iter = 0;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        
        gsl_function F;
        Params_Inversion_PTX params;
        params.sw = this;
        params.P = P0;

        F.function = &func_T_VLH;
        F.params = &params;

        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        // printf ("using %s method\n", gsl_root_fsolver_name (s));
        // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
        do
        {
            iter++;
            status = gsl_root_fsolver_iterate (s);
            T_ = gsl_root_fsolver_root (s);
            x_lo = gsl_root_fsolver_x_lower (s);
            x_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/1E4);
            // if (status == GSL_SUCCESS) printf ("Converged:\n");
            // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_crit);
        }
        while (status == GSL_CONTINUE && iter < ITERATION_MAX);

        if(!((status == GSL_SUCCESS)))
        {
            printf ("status = %s\n\n", gsl_strerror (status));
            ERROR("Fatal error in double cH2ONaCl::T_VLH_P0");
        }
        gsl_root_fsolver_free (s);
        return T_;
    }

        /**
         * @brief Calculate the maximum temperature of the V+L+H zone. 
         * 
         * It should be the root of \f$ P_{LVH} = P_{min} \f$, see also Eq. 10 of \cite Driesner2007Part1 and \link P_VLH \endlink
         * 
         * \note This function is only used to calculate the constant value.
         * 
         * @return double [K]
         */
        double cH2ONaCl::Tmax_VLH()
        {
            return T_VLH_P0_(pmin(), T_Peak_VLH, Tmax());
        }
        double cH2ONaCl::Tmin_VLH()
        {
            return T_VLH_P0_(pmin(), Tmin(), T_Peak_VLH);
        }

        /**
         * @brief Calculate temperature boundary of VLH zone by given pressure P0.
         *
         * \warning Only valid in pressure range [H2ONaCl::P_MIN, H2ONaCl::P_Peak_VLH]
         *
         * \note The pressure value should be checked in the function, for example to avoid program crash happening when P>P_Peak_VLH,
         * while the if statement check could decrease program performance if this function is called in large for loops.
         * Therefore, please be careful the valid pressure range!
         * 
         * @param P0 [Pa]
         * @param Tmin [K]
         * @param Tmax [K]
         */
        void cH2ONaCl::T_VLH_P0(double P0, double& T_min, double& T_max)
        {
            T_min = T_VLH_P0_(P0, Tmin(), T_Peak_VLH);
            T_max = T_VLH_P0_(P0, T_Peak_VLH, Tmax());

        }

        /**
         * @brief Calculate \f$ log_{10}K^{\prime} \f$ in Eq. 14 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param P_NaCl [Pa] Vapor pressure of NaCl
         * @return double 
         */
        double cH2ONaCl::Log10_Kprime(const double& T, const double& P, double& P_NaCl)
        {
            double T_C = T - 273.15;
            double j[4] ={0,0,0,0};
            // 1. calculate parameters j_i in the Table 8 of Driesner(2007a)
            j[0] = m_tab8.k[0] + m_tab8.k[1] * exp(-m_tab8.k[2] * T_C);
            j[1] = m_tab8.k[4] + (m_tab8.k[3] - m_tab8.k[4]) / (1 + exp((T_C - m_tab8.k[5]) / m_tab8.k[6])) + m_tab8.k[7] * pow(T_C + m_tab8.k[8], 2.0);
            for (size_t i = 0; i < 4; i++)
            {
                j[2] += m_tab8.k[i+9]*pow(T_C,i);
            }
            for (size_t i = 0; i < 3; i++)
            {
                j[3] += m_tab8.k[i+13]*pow(T_C,i);
            }
            // 2. calculate vapor pressure of pure NaCl
            P_NaCl = m_NaCl->P_Vapor(T);  //[Pa]
            // 3. calculate critical pressure 
            double P_crit = 0;  // [Pa]
            P_Critical(T, P_crit);
            // 4. Calculate normalized pressure: Eq. 16
            double P_normalized = (P - P_NaCl)/(P_crit - P_NaCl); //always <1 in the valid T,P range of V+H zone.
            // if(P_normalized>=1.0)std::cout<<"P_normalized: "<<P_normalized<<std::endl;
            if(P_normalized>1)P_normalized=1;
            // assert(P_normalized<=1.0);

            double one_minus_P_normalized = 1 - P_normalized; //used in eq. 17
            // 5. eq. 17
            double log10_K_bar = 1 + j[0]*pow(one_minus_P_normalized, j[1])
                                        + j[2]*one_minus_P_normalized
                                        + j[3]*pow(one_minus_P_normalized, 2)
                                        - (1 + j[0] + j[2] + j[3])*pow(one_minus_P_normalized, 3);
            // 6. eq. 15              
            double log10_XL_P_NaCl = log10(X_HaliteLiquidus_mol(T, P_NaCl)); //used in eq. 15
            double log10_K_prime = log10_K_bar * (log10(P_NaCl/P_crit) - log10_XL_P_NaCl) + log10_XL_P_NaCl; //eq. 15
            return log10_K_prime;
        }

        /**
         * @brief Halite-saturated vapor and/or liquid composition, see equation 8,9,14-17 of \cite Driesner2007Part1.
         * 
         * \see \link X_HaliteLiquidus \endlink
         * 
         * The correlations for distribution coefficient \f$ K \f$ used to calculate the halite-saturated vapor composition is defined as (eq. 9 of ref. \cite Driesner2007Part1),
         * \f{equation}
         * \frac{X_{NaCl, sat}^{L, metastable}}{X_{NaCl, sat}^V} = K
         * \f}
         * A modified distribution coefficient \f$ K^{\prime} \f$ (this is actually used to compute halite-saturated vapor composition) is defined as (eq. 14 of ref. \cite Driesner2007Part1), 
         * 
         * \f{equation}
         * log_{10}K^{\prime} = log_{10}\left( \frac{x_l}{x_v/(\frac{P_{NaCl}}{P})} \right) = log_{10} \left( \frac{x_l}{x_v} \right) + log_{10} \left( \frac{P_{NaCl}}{P} \right)
         * \f}
         * where \f$ x_l \f$ is calculated from #X_HaliteLiquidus, \f$ P_{NaCl} \f$ is the boling pressure or sublimation pressure depends on temperature, it can be expressed as,
         * 
         * \f{equation}
         *      P_{NaCl} = \left\{ \begin{matrix}
         *      log_{10}(P_{NaCl, liquid}),& T > T_{triple, NaCl} \\ \\[1ex]
         *      log_{10}(P_{NaCl, halite}),& \text{else}\\ \\[1ex]
         *      \end{matrix}\right.
         * \f}
         * \f$ P_{NaCl} \f$ can be calculated from \link NaCl::cNaCl::P_Vapor \endlink. And \f$ log_{10}K^{\prime} \f$ can be computed from following equation (eq. 15 of ref. \cite Driesner2007Part1),
         * 
         * \f{equation}
         * log_{10}\bar K = \frac{log_{10}K^{\prime} - log_{10}\left( X_{NaCl, sat}^L \right)_{P_{NaCl}}}{log_{10}\left( \frac{P_{NaCl}}{P_{crit}} \right) - log_{10}\left( X_{NaCl, sat}^L \right)_{P_{NaCl}}} = 1 + j_0(1-\bar P)^{j_1} + j_2(1-\bar P) + j_3(1-\bar P)^{2} - (1 + j_0 + j_2 + j_3)(1-\bar P)^{3}
         * \f}
         * where \f$ P_{crit} \f$ is calculated from #P_X_Critical, \f$ X_{NaCl, sat}^L \f$ is calculated from #X_HaliteLiquidus (T, \f$ P_{NaCl} \f$), \f$ \bar P = (P - P_{NaCl})/(P_{crit} - P_{NaCl})\f$ is the normalized pressure, \f$ j_i (i=0, 1, 2,, 3) \f$ is calculated from Table 8 of reference \cite Driesner2007Part1. 
         * 
         * \note The normalized pressure \f$ \bar P <1 \f$ in the V+H zone, because the critical curve always above the V+H zone. See \link P_VLH \endlink.
         * 
         * \warning The valid pressure and temperature range of this function are [#T_MIN, #T_MAX_VLH] and [#P_MIN, #P_Peak_VLH], respectively.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X_L [mol/mol] X_L is just a place holder for VH zone, but it is useful for VLH zone.
         * @param X_V [mol/mol] Vapor composition at halite saturated zone: VH, VHL.
         * 
         * @return double [mol/mol]
         */
        void cH2ONaCl::X_HaliteSat_mol(const double& T, const double& P, double& X_L, double& X_V)
        {
            double P_NaCl;
            double log10_K_prime = Log10_Kprime(T,P, P_NaCl);
            // 7. eq. 14
            double log10_XLbyXV = log10_K_prime - log10(P_NaCl/P); //eq. 14
            X_L = X_HaliteLiquidus_mol(T, P); //metastable extension  of the halite liquidus to low pressur (P<P_VLH), still use Eq.8.
            X_V = X_L/pow(10, log10_XLbyXV); //eq. 14
        }

        /**
         * @brief Halite-saturated vapor composition at V+H coexistence, see equation 9 and Fig. 8 of \cite Driesner2007Part1.
         * 
         * \image html H2ONaCl/H2ONaCl_X_HaliteSatVapor.svg Halite-saturated vapor composition. (a) Isotherms; (b) Isobars. width=70%
         * 
         * \warning The Fig.8a of Driesner(2007) looks shift 1E-1 to the left.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/kg] 
         */
        double cH2ONaCl::X_VH(const double& T, const double& P)
        {
            double X_L = 0, X_V = 0;
            X_HaliteSat_mol(T,P, X_L, X_V);
            return Mol2Wt(X_V);
        }

        /**
         * @brief Calculate salinity of vapor and liquid at VLH coexistence.
         * 
         * Use the same formula with halite-saturated vapor composition at VH coexistence, see \link X_VH \endlink and \link X_HaliteSat_mol \endlink. Where \link X_HaliteSat_mol \endlink can calculate both \f$ X_l \f$ and \f$ X_v \f$ at halite-saturated domain.
         * 
         * \image html H2ONaCl/H2ONaCl_X_VLH.svg Composition at VLH coexistence. (a) Liquid composition; (b) Vapor composition. size=80%
         * 
         * \see Fig. 10 and 11 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X_L [kg/kg]
         * @param X_V [kg/kg]
         */
        void cH2ONaCl::X_VLH(const double& T, const double& P, double& X_L, double& X_V)
        {
            X_HaliteSat_mol(T,P, X_L, X_V);
            X_L = Mol2Wt(X_L);
            X_V = Mol2Wt(X_V);
        }
        
        /**
         * @brief Calculate liquid composition on liquid branch of V+L coexistence.
         * 
         * See Eqn. 11 and table 7 of \cite Driesner2007Part1.
         * 
         * \image html H2ONaCl/H2ONaCl_X_VL.svg Isothermal sections of the V+L surface. Composition (mole fraction) is given on a logarithmic scale. width=100%
         *
         *
         * \see Fig. 12 of \cite Driesner2007Part1.
         *
         * \warning For T smaller than critical point of H2O \f$ T<T_{H2O} \f$, the valid maximum pressure of this function,
         * i.e. the upper boundary of VL surface, is the boiling pressure of H2O (\f$ p_{H2O} \f$) for a given T,
         * not the dummy critical pressure(\f$ p_{crit} \f$) given by Eq.5c of \cite Driesner2007Part1.
         * Because the \f$ p_{H2O} \f$ always lower than the \f$ p_{crit} \f$,
         * and the salinity at \f$ p_{H2O} \f$ equal to zero is introduced,
         * so if one calculate Xl or Xv for pressure greater than \f$ p_{H2O} \f$ (it is not on the VL surface), the result must be negative.
         * Actually, it doesn't make sense to calculate this.
         * Therefore, for phase boundary calculation (e.g. \link PhaseBoundary_VL_DeformLinear \endlink), have to set a correct upper pressure (must be \f$ p_{H2O}(T) \f$) of the surface.
         *
         * * How to calculate \f$g_0\f$ and the dummy \f$ X_{crit} \f$
         *
         * Use two constrains:
         *
         * (1) \f$ \color{red}{X_{NaCl}^{VL,liq} = X_{NaCl,sat}^{L}} \f$ (\link X_HaliteLiquidus_mol \endlink) at the lower pressure bound \f$ P_{VLH} \f$ (\link P_VLH \endlink);
         *
         * (2) \f$ \color{red}{X_{NaCl}^{VL,liq} = 0} \f$ at the upper pressure bound \f$ p_{H2O} \f$ (e.g., \link cIAPS84::Boiling_p \endlink);
         *
         * \f{align}
         * \left.
         * \begin{matrix}
         *  {\color{red}{X_{NaCl,sat}^{L}}} & = & X_{crit} + g_0\sqrt{\color{blue}{P_{crit} - P_{vlh}}} + g_1({\color{blue}{P_{crit} - P_{vlh}}}) + g_2({\color{blue}{P_{crit} - P_{vlh}}})^2 \\
         *  {\color{red}{0}} & = & X_{crit} + g_0\sqrt{\color{green}{P_{crit} - P_{H_2O}}} + g_1({\color{green}{P_{crit} - P_{H_2O}}}) + g_2({\color{green}{P_{crit} - P_{H_2O}}})^2
         * \end{matrix}
         * \right \} \Longrightarrow
         * \left\{
         * \begin{matrix}
         * g_0 & = & \frac{{\color{red}{X_{NaCl,sat}^{L}}} + g_1({\color{green}{P_{crit} - P_{H_2O}}}) + g_2({\color{green}{P_{crit} - P_{H_2O}}})^2 -g_1({\color{blue}{P_{crit} - P_{vlh}}}) - g_2({\color{blue}{P_{crit} - P_{vlh}}})^2 }{\sqrt{\color{blue}{P_{crit} - P_{vlh}}} - \sqrt{\color{green}{P_{crit} - P_{H_2O}}}} \\
         * X_{crit} & = & - g_0\sqrt{\color{green}{P_{crit} - P_{H_2O}}} - g_1({\color{green}{P_{crit} - P_{H_2O}}}) - g_2({\color{green}{P_{crit} - P_{H_2O}}})^2
         * \end{matrix}
         * \right.
         * \f}
         *
         * @param T [K]
         * @param P [Pa]
         * @return double [mol/mol]
         */
        double cH2ONaCl::XL_VL_mol(const double& T, const double& P)
        {
            double T_C = T - 273.15;
            double g0 = 0;
            double g1 = m_table7.h[1] + (m_table7.h[0] - m_table7.h[1])/(1.0 + exp((T_C - m_table7.h[2])/m_table7.h[3])) + m_table7.h[4] * T_C * T_C;
            double g2 = m_table7.h[6] + (m_table7.h[5] - m_table7.h[6])/(1.0 + exp((T_C - m_table7.h[7])/m_table7.h[8])) + m_table7.h[9]*exp(-m_table7.h[10] * T_C);

            double P_crit = 0;
            P_Critical(T, P_crit);

            double p_vlh = 0, P_boiling_H2O = 0, X_vlh = 0;
            if(T<NaCl::T_Triple)
            {
                p_vlh = P_VLH(T);
                X_vlh = X_HaliteLiquidus_mol(T, p_vlh);
            }else if(T == NaCl::T_Triple)
            {
                p_vlh = NaCl::P_Triple;
                X_vlh = 1;
            }else
            {
                p_vlh = 0;
                X_vlh = 1;
            }
            double P_crit_minus_P = P_crit - P;
            if (P_crit_minus_P<0) //this value will be used in Eq. 11 as sqrt(P_crit_minus_P), so have to check its sign. Actually, this function is only valid when P_crit_minus_P>=0
            {
                P_crit_minus_P = fabs(P_crit_minus_P);
                if(P_crit_minus_P>0.001) //Pa
                {
                    // throw OutOfRangeError("The input pressure of function XL_VL_mol(const double& T, const double& P) is greater than critical pressure, it's not valid. T="+std::to_string(T-273.15)+"deg.C, p="+std::to_string(P/1E5)+"bar, P_crit-P="+std::to_string(P_crit-P)+"Pa");
                    #if DEBUG
                    WARNING("The input pressure of function XL_VL_mol(const double& T, const double& P) is greater than critical pressure, it's not valid. T="+std::to_string(T-273.15)+"deg.C, p="+std::to_string(P/1E5)+"bar, P_crit="+std::to_string(P_crit/1E5)+"bar, P_crit-P="+std::to_string(P_crit-P)+"Pa");
                    #endif
                }
            }
            double P_crit_minus_Pvlh = P_crit - p_vlh;
            double X_VL_LiquidBranch = 0;
            if(T < m_constants_Water.T_critical)
            {

                P_boiling_H2O = m_Water->Boiling_p(T);
                if(P >= P_boiling_H2O) //If the input pressure very close to boling pressure of H2O, then Xl=0. For example, in mmc4 table of Driesner(2007), T=3.200000e+02, P=1.127932e+02, then X=0. If do not use this check, this function will calculate a negative Xl, and then the density and other properties will be wrong.
                {
                    X_VL_LiquidBranch = 0;
                    if (P_boiling_H2O - P < 100.0)  // OK, 100 Pa difference would be caused by some numerical errors(e.g. 四舍五入), for example T=3.200000e+02, P=1.127932e+02 in mmc4 table.
                    {
                        #if DEBUG
                        static bool warning_is_shown = false;
                        if(!warning_is_shown) {
                            WARNING("The input pressure of function XL_VL_mol(const double& T, const double& P) is greater than P_boiling_H2O in sub-critical region of H2O, it's not valid. T=" +
                                    std::to_string(T - 273.15) + "deg.C, p=" + std::to_string(P / 1E5) +
                                    "bar, P_boiling_H2O=" + std::to_string(P_boiling_H2O / 1E5) +
                                    "bar, P_boiling_H2O-P=" + std::to_string(P_boiling_H2O - P) + "Pa");
                            STATUS("This warning information will only shown once");
                            warning_is_shown = true;
                        }
                        #endif
                    } else //Otherwise, I have to stop the program, one could make a big mistake to use this function, have to let user check the input pressure range.
                    {
                        throw ValueError ("The input pressure of function XL_VL_mol(const double& T, const double& P) is greater(>100 Pa) than P_boiling_H2O in sub-critical region of H2O, it's not valid. Please check the input parameters and the manual using this function. T="+std::to_string(T-273.15)+"deg.C, p="+std::to_string(P/1E5)+"bar, P_boiling_H2O="+std::to_string(P_boiling_H2O/1E5)+"bar, P_boiling_H2O-P="+std::to_string(P_boiling_H2O-P)+"Pa");
                    }

                } else
                {
                    double P_crit_minus_P_H2O = P_crit - P_boiling_H2O;
                    double g1g2_P_H2O = g1*P_crit_minus_P_H2O + g2*pow(P_crit_minus_P_H2O, 2.0);
                    double g1g2_P_VLH = g1*P_crit_minus_Pvlh + g2*pow(P_crit_minus_Pvlh, 2.0);
                    double P_crit_minus_P_H2O_sqr = sqrt(P_crit_minus_P_H2O);
                    g0 = (X_vlh + g1g2_P_H2O - g1g2_P_VLH)/(sqrt(P_crit_minus_Pvlh) - P_crit_minus_P_H2O_sqr);
                    double X_crit = -g0*P_crit_minus_P_H2O_sqr  - g1g2_P_H2O;
                    X_VL_LiquidBranch = X_crit + g0* sqrt(P_crit_minus_P) + g1*P_crit_minus_P + g2*pow(P_crit_minus_P, 2.0);
                }
            }else
            {
                double X_crit = 0;
                X_Critical_mol(T, X_crit);
                g0 = (X_vlh - X_crit - g1 * (P_crit_minus_Pvlh) - g2 * pow(P_crit_minus_Pvlh, 2.0)) / sqrt(P_crit_minus_Pvlh);
                X_VL_LiquidBranch = X_crit + g0 * sqrt(P_crit_minus_P) + g1 * P_crit_minus_P + g2 * pow(P_crit_minus_P, 2.0);
            }
        
            return X_VL_LiquidBranch;
        }

        /**
         * @brief Calculate liquid composition on liquid branch of V+L coexistence.
         * 
         * \see \link XL_VL_mol \endlink
         * 
         * \image html H2ONaCl/H2ONaCl_X_VL.svg Isothermal sections of the V+L surface. Composition (mole fraction) is given on a logarithmic scale. width=100%
         * 
         * \see Fig. 12 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/kg]
         */
        double cH2ONaCl::XL_VL(const double& T, const double& P)
        {
            return Mol2Wt(XL_VL_mol(T,P));
        }

        double func_T_VL_L(double T, void *params)
        {
            Params_Inversion_PTX* param = (Params_Inversion_PTX*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P;
            double X = param->X;
            double xl_vl;
            xl_vl = sw->XL_VL(T,p);
            return xl_vl - X;
        }
        /**
         * Calculate temperature on surface VL-liquid branch for given p, X. Used in HPX calculation.
         * @param p [Pa]
         * @param X [kg/kg]
         * @return
         */
        double cH2ONaCl::T_VL_L(const double& p, const double& X, const double& Tmin, const double& Tmax)
        {
            double T_ = 0;
            double x_lo = Tmin, x_hi = Tmax;
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;

            gsl_function F;
            Params_Inversion_PTX params;
            params.sw = this;
            params.P = p;
            params.X = X;

            F.function = &func_T_VL_L;
            F.params = &params;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                T_ = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/1E4);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
                printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            if(!((status == GSL_SUCCESS)))
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                ERROR("Fatal error in double cH2ONaCl::T_VLH_P0");
            }
            gsl_root_fsolver_free (s);
            return T_;
        }

        /**
         * @brief Calculate vapor composition on vapor branch of V+L coexistence.
         * 
         * \image html H2ONaCl/H2ONaCl_X_VL.svg Isothermal sections of the V+L surface. Composition (mole fraction) is given on a logarithmic scale. width=100%
         * 
         * \see Fig. 12 of \cite Driesner2007Part1. 
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [mol/mol] 
         */
        double cH2ONaCl::XV_VL_mol(const double& T, const double& P)
        {
            double X_VL_LiquidBranch = XL_VL_mol(T,P);
            double P_NaCl;
            double log10_K_prime = Log10_Kprime(T,P, P_NaCl);
            double K = pow(10, log10_K_prime - log10(P_NaCl/P));
            double X_VL_VaporBranch = X_VL_LiquidBranch/K; //eq. 14
            return X_VL_VaporBranch;
        }
        
        /**
         * @brief Calculate vapor composition on vapor branch of V+L coexistence. 
         * 
         * \image html H2ONaCl/H2ONaCl_X_VL.svg Isothermal sections of the V+L surface. Composition (mole fraction) is given on a logarithmic scale. width=100%
         * 
         * \see Fig. 12 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/kg]
         */
        double cH2ONaCl::XV_VL(const double& T, const double& P)
        {
            return Mol2Wt(XV_VL_mol(T,P));
        }

        double func_T_VL_V(double T, void *params)
        {
            Params_Inversion_PTX* param = (Params_Inversion_PTX*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P;
            double X = param->X;
            double xl_vl;
            xl_vl = sw->XV_VL(T,p);
            return xl_vl - X;
        }
        /**
         * Inverse T at VL->Vapor branch surface for give p,X
         * @param p [Pa]
         * @param X [kg/kg]
         * @param Tmin [K]
         * @param Tmax [K]
         * @return
         */
        double cH2ONaCl::T_VL_V(const double& p, const double& X, const double& Tmin, const double& Tmax)
        {
            double T_ = 0;
            double x_lo = Tmin, x_hi = Tmax;
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;

            gsl_function F;
            Params_Inversion_PTX params;
            params.sw = this;
            params.P = p;
            params.X = X;

            F.function = &func_T_VL_V;
            F.params = &params;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                T_ = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/1E4);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
                printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            if(!((status == GSL_SUCCESS)))
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                ERROR("Fatal error in double cH2ONaCl::T_VLH_P0");
            }
            gsl_root_fsolver_free (s);
            return T_;
        }

        /**
         * @brief Calculate \f$ X_{V}, X_{L} \f$ on V+L surface in one step.
         *
         * Because \link XV_VL \endlink always need to call \link XL_VL \endlink first, so this function can save computing time if both \f$ X_{V} \f$ and \f$ X_{L} \f$ are needed.
         *
         * 
         * \image html H2ONaCl/H2ONaCl_X_VL.svg Isothermal sections of the V+L surface. Composition (mole fraction) is given on a logarithmic scale. width=100%
         * 
         * \see Fig. 12 of \cite Driesner2007Part1.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X_L [kg/kg]
         * @param X_V [kg/kg]
         */
        void cH2ONaCl::X_VL(const double& T, const double& P, double& X_L, double& X_V)
        {
            double X_VL_LiquidBranch = XL_VL_mol(T,P);
            double P_NaCl;
            double log10_K_prime = Log10_Kprime(T,P, P_NaCl);
            double X_VL_VaporBranch = X_VL_LiquidBranch/pow(10, log10_K_prime) * P_NaCl/P; //eq. 14
            if(P<=P_VLH(T))
            {
                X_VL_VaporBranch = X_HaliteLiquidus_mol(T,P)/X_VL_LiquidBranch * X_VL_VaporBranch;
            }
            X_L = Mol2Wt(X_VL_LiquidBranch);
            X_V = Mol2Wt(X_VL_VaporBranch);
        }

        /**
         * @brief Calculate phase region and composition of liquid and vapor start from the node of "Hallite liquidus" in the flow chart shown in \link getPhaseRegion \endlink (edge in red).
         * 
         * This is a function that is used repeatedly in \link getPhaseRegion \endlink.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X Bulk salinity. [kg/kg]
         * @param phase_region Phase region index
         * @param X_V Salinity of vapor phase if it exists. [kg/kg] 
         * @param X_L Salinity of liquid phase if it exists. [kg/kg] 
         */
        void cH2ONaCl::_getPhaseRegion_node_HaliteLiquidus(const double& T, const double& P, const double& X, PhaseRegion& phase_region, double& X_V, double& X_L)
        {
            double X_LH = X_HaliteLiquidus(T,P);
            if(X < X_LH)
            {
                phase_region = SinglePhase_L;
                X_V = 0; X_L = X;
            }else
            {
                double T_hm = m_NaCl->Melting_T(P);
                if(T>=T_hm)
                {
                    phase_region = SinglePhase_L;
                    X_V = 0; X_L = X;
                } else
                {
                    phase_region = TwoPhase_LH;
                    X_V = 0; X_L = X_LH;
                }
            }
        }
        
        /**
         * @brief Calculate phase region and composition of liquid and vapor start from the node of "check T<T_crit(H2O)" in the flow chart shown in \link getPhaseRegion \endlink (edge in red).
         * 
         * This is a function that is used repeatedly in \link getPhaseRegion \endlink.
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X Bulk salinity. [kg/kg]
         * @param phase_region Phase region index
         * @param X_V Salinity of vapor phase if it exists. [kg/kg] 
         * @param X_L Salinity of liquid phase if it exists. [kg/kg] 
         */
        void cH2ONaCl::_getPhaseRegion_node_CheckTcrit_H2O(const double& T, const double& P, const double& X, PhaseRegion& phase_region, double& X_V, double& X_L)
        {
            if(T <= m_constants_Water.T_critical)
            {
                phase_region = SinglePhase_V;
                X_V = X; X_L = 0;
            }else
            {
                phase_region = SinglePhase_L;
                X_V = 0; X_L = X;
            }
        }
        
        /**
         * @brief Calculate phase region and composition of liquid and vapor start from the node of "X_l, X_v of V+L" in the flow chart shown in \link getPhaseRegion \endlink (edge in red).
         * 
         * This is a function that is used repeatedly in \link getPhaseRegion \endlink.
         *
         * \warning How to deal with the case of (X<X_crit, P>P_crit, T>T_crit_h2o) ? In \cite Driesner2007Part1 (e.g., Fig.2), this region is described as "F", but in the matlab code, this region is classified as "Vapor".
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X Bulk salinity. [kg/kg]
         * @param X_crit Calculated critical composition. [kg/kg]
         * @param phase_region Phase region index
         * @param X_V Salinity of vapor phase if it exists. [kg/kg] 
         * @param X_L Salinity of liquid phase if it exists. [kg/kg] 
         */
        void cH2ONaCl::_getPhaseRegion_node_XVXL_VL(const double& T, const double& P, const double& X, const double& X_crit, PhaseRegion& phase_region, double& X_V, double& X_L)
        {
            double x_v, x_l;
            X_VL(T, P, x_l, x_v);
            // (1.2.1) compare X and [x_v, x_l]
            if((X >= x_v) && (X <= x_l))
            {
                phase_region = TwoPhase_VL;
                X_L = x_l; X_V = x_v;
            }else if(X > x_l)
            {
                // (1.2.2) compare with halite liquidus
                _getPhaseRegion_node_HaliteLiquidus(T, P, X, phase_region, X_V, X_L);
            }else  //X<x_v
            {
                // (1.2.3) compare T with critical T of H2O
                // _getPhaseRegion_node_CheckTcrit_H2O(T,P,X,phase_region, X_V, X_L);
                phase_region = SinglePhase_V;
                X_L = 0; X_V = X;
            }
        }
    
        /**
         * @brief Calculate phase region index and composition of vapor and liquid phase if the phase exists.
         * 
         * \dotfile getPhaseRegion.dot
         * 
         * @param T [K]
         * @param P [Pa]
         * @param X Bulk salinity. [kg/kg]
         * @param phase_region Phase region index
         * @param X_V Salinity of vapor phase if it exists. [kg/kg] 
         * @param X_L Salinity of liquid phase if it exists. [kg/kg] 
         * 
         * \image html H2ONaCl/H2ONaCl_isothermal_PhaseRegion.svg Isothermal phase regions. width=100%
         * 
         */
        void cH2ONaCl::getPhaseRegion_TPX(const double& T, const double& P, const double& X, PhaseRegion& phase_region, double& X_V, double& X_L)
        {
            // first check if it is pure water.
            if(X==0)
            {
                ThermodynamicProperties props;
                m_Water->UpdateState_TPX(props, T, P);
                phase_region = props.phase;
                X_V = 0; X_L = 0;
                return;
            }

            // 1. calculate critical pressure and critical component.
            double P_crit, X_crit;
            P_X_Critical(T, P_crit, X_crit);
            // (1) compare with critical pressure
            if (P > P_crit)
            {
                double T_halite_melting = m_NaCl->Melting_T(P);
                // (1.1) compare with melting T of NaCl
                if(T > T_halite_melting)
                {
                    phase_region = SinglePhase_L;
                    X_V = 0; X_L = X;
                }else
                {
                    // (1.1.1) compare with halite liquidus
                    _getPhaseRegion_node_HaliteLiquidus(T, P, X, phase_region, X_V, X_L);
                }
            }else
            {
                // (1.2) compare T with Tmax of VLH zone
                if (T>T_MAX_VLH)
                {
                    _getPhaseRegion_node_XVXL_VL(T,P,X,X_crit, phase_region, X_V, X_L);
                }else
                {
                    double P_vlh = P_VLH(T);
                    if(P > P_vlh)
                    {
                        _getPhaseRegion_node_XVXL_VL(T,P,X,X_crit, phase_region, X_V, X_L);
                    }else // if(P < P_vlh) ignore VLH region in PTX space
                    {
                        double X_vh = X_VH(T, P);
                        if(X < X_vh)
                        {
                            // _getPhaseRegion_node_CheckTcrit_H2O(T,P,X,phase_region, X_V, X_L);
                            phase_region = SinglePhase_V;
                            X_L = 0; X_V = X;
                        }else
                        {
                            phase_region = TwoPhase_VH;
                            X_V = X_vh; X_L = 0;
                        }
                    }
                    // ignore Three phase in PTX space
                    // else
                    // {
                    //     phase_region = ThreePhase_V_L_H;
                    //     X_VLH(T, P, X_L, X_V);
                    // }
                }
                
            }
        }

        /**
         * @brief Calculate water density
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/m3]
         */
        double cH2ONaCl::_Rho_water(const double& T, const double& P, double& dRhodP, double& dRhodT, double& IsothermalCompresibility, double& IsobaricExpansivity, PhaseType phase)
        {
            ThermodynamicProperties props;
            m_Water->UpdateState_TPX(props, T, P);
            // DEBUG: check if return two phase when input as T,P; Make the logics easier, this function always return single phase and its properties because the input is T and P rather than H,P, so two phase is not allowed, if it happens, need to modify related backends of H2O
            if (props.phase==TwoPhase_VL_Water)
            {
                throw NotImplementedError("In function H2ONaCl::_Rho_water, input T,P but in two phase region, should modify UpdateState_TPX function of each backend of H2O and avoid this case happens. H2O EOS is  "+m_Water->name()+", phase = "+std::to_string(props.phase) + ", name "+ phase_name(props.phase) + ", input T="+std::to_string(T)+" K, p="+std::to_string(P));
            }
            dRhodT = props.dRhodT;
            dRhodP = props.dRhodP;
            IsothermalCompresibility = props.IsothermalCompressibility;
            IsobaricExpansivity = props.IsobaricExpansivity;
            return props.Rho;

            // // std::cout<<"boilp: "<<m_Water->Boiling_p(T)<<" "<<P<<" "<<m_Water->Boiling_p(T)-P<<std::endl;
            // ThermodynamicProperties props;
            // double rho_l, rho_v;
            // if(T<=m_constants_Water.T_critical && fabs(m_Water->Boiling_p(T, rho_l, rho_v)-P)<(1E-6*P)) // CoolProp: if p-T close to boiling p, UpdateState_TP will throw errors. If dp<1E-4 of given P, it is regarded as two phase in coolprop. See line 1992 of Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp
            // {
            //     switch (phase) {
            //         case Liquid:
            //             // return m_Water->props_sat_liquid().Rho;
            //             return rho_l;
            //             break;
            //         case Vapor:
            //             // return m_Water->props_sat_vapor().Rho;
            //             return rho_v;
            //             break;
            //     }
            // } else
            // {
            //     m_Water->UpdateState_TPX(props, T, P);
            //     if (props.phase==TwoPhase_VL_Water) //This is only for PROST, because PROST check the TWO phase region, if the P,T condition in two phase region, the prop->d and prop->h will both equal 0
            //     {
            //         m_Water->Boiling_p(T, props);
            //         switch (phase) {
            //             case Liquid:
            //                 // return m_Water->props_sat_liquid().Rho;
            //                 dRhodT = props.dRhodT_l;
            //                 dRhodP = props.dRhodP_l;
            //                 IsothermalCompresibility = props.IsothermalCompressibility_l;
            //                 IsobaricExpansivity = props.IsobaricExpansivity_l;
            //                 return props.Rho_l;
            //                 break;
            //             case Vapor:
            //                 // return m_Water->props_sat_vapor().Rho;
            //                 dRhodT = props.dRhodT_v;
            //                 dRhodP = props.dRhodP_v;
            //                 IsothermalCompresibility = props.IsothermalCompressibility_v;
            //                 IsobaricExpansivity = props.IsobaricExpansivity_v;
            //                 return props.Rho_v;
            //                 break;
            //         }
            //     }
            //     dRhodT = props.dRhodT;
            //     dRhodP = props.dRhodP;
            //     IsothermalCompresibility = props.IsothermalCompressibility;
            //     IsobaricExpansivity = props.IsobaricExpansivity;
            //     return props.Rho;
            // }


        }

        double cH2ONaCl::_Rho_water(const double& T, const double& P, PhaseType phase)
        {
            double tmp_drhodp, tmp_drhodt, tmp_kappa, tmp_beta;
            return _Rho_water(T, P, tmp_drhodp, tmp_drhodt, tmp_kappa, tmp_beta, phase);
        }

        /**
         * @brief Calculate specific enthalpy of water (H2O)
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [J/kg]
         */
        double cH2ONaCl::_H_water(const double& T, const double& P, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(T<=m_constants_Water.T_critical && fabs(m_Water->Boiling_p(T, props)-P)<(1E-6*P)) // !!! CoolProp: if p-T close to boiling p, UpdateState_TP will throw errors.
            {
                switch (phase) {
                    case Liquid:
                        // return m_Water->props_sat_liquid().Rho;
                        return props.H_l;
                        break;
                    case Vapor:
                        // return m_Water->props_sat_vapor().Rho;
                        return props.H_v;
                        break;
                }
            } else
            {
                m_Water->UpdateState_TPX(props, T, P);
                if (props.phase==TwoPhase_VL_Water)
                {
                    m_Water->Boiling_p(T, props);
                    switch (phase) {
                        case Liquid:
                            // return m_Water->props_sat_liquid().H;
                            return props.H_l;
                            break;
                        case Vapor:
                            // return m_Water->props_sat_vapor().H;
                            return props.H_v;
                            break;
                    }
                }
                // return m_Water->hmass();
                return props.H;
            }
        }
        void cH2ONaCl::_H_Cp_water(const double& T, const double& P, double& H, double& Cp, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(T<=m_constants_Water.T_critical && fabs(m_Water->Boiling_p(T, props)-P)<(1E-6*P)) // !!! CoolProp: if p-T close to boiling p, UpdateState_TP will throw errors.
            {
                switch (phase) {
                    case Liquid:
                        Cp = props.Cp_l;
                        H = props.H_l;
                        break;
                    case Vapor:
                        Cp = props.Cp_v;
                        H = props.H_v;
                        break;
                }
            } else
            {
                m_Water->UpdateState_TPX(props, T, P);
                if (props.phase==TwoPhase_VL_Water)
                {
                    Cp = 0;
                    m_Water->Boiling_p(T, props);
                    switch (phase) {
                        case Liquid:
                            // return m_Water->props_sat_liquid().H;
                            H = props.H_l;
                            break;
                        case Vapor:
                            // return m_Water->props_sat_vapor().H;
                            H = props.H_v;
                            break;
                    }
                }
                Cp = props.Cp;
                H = props.H;
            }
        }
        /**
         * Calculate water viscosity for give T,P.
         * @param T [K]
         * @param P [Pa]
         * @param phase [Liquid|Vapor]
         * @return
         */
        double cH2ONaCl::_Mu_water(const double& T, const double& P, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(T<=m_constants_Water.T_critical && fabs(m_Water->Boiling_p(T, props)-P)<(1E-6*P)) // !!! CoolProp: if p-T close to boiling p, UpdateState_TP will throw errors.
            {
                switch (phase) {
                    case Liquid:
                        // return m_Water->props_sat_liquid().Rho;
                        return props.Mu_l;
                        break;
                    case Vapor:
                        // return m_Water->props_sat_vapor().Rho;
                        return props.Mu_v;
                        break;
                }
            } else
            {
                m_Water->UpdateState_TPX(props, T, P);
                if (props.phase==TwoPhase_VL_Water)
                {
                    m_Water->Boiling_p(T, props);
                    switch (phase) {
                        case Liquid:
                            // return m_Water->props_sat_liquid().H;
                            return props.Mu_l;
                            break;
                        case Vapor:
                            // return m_Water->props_sat_vapor().H;
                            return props.Mu_v;
                            break;
                    }
                }
                // return m_Water->hmass();
                return props.Mu;
            }
        }
        
        /**
         * @brief Calculate coefficient \f$n_1, n_2 \f$ for \f$T^{*} \f$ of density calculation (Eq. 4 of \cite Driesner2007Part2).
         * 
         * \f$ n_1, n_2 \f$ can be calculated from Eq. 9 of \cite Driesner2007Part2. 
         * While, the coefficients \f$ \color{red}{n_{10}, n_{12}, n_{20}, n_{23}} \f$ are not given in Tab. 4 of \cite Driesner2007Part2.
         * So, we need to calculate these four coefficients using some constraint conditions. 
         * 
         * * (1) \f$ n_{10} \f$ can be calculated by combing Eq.11 of \cite Driesner2007Part2 and condition \f$ X = 1 \f$: \f$ {\color{red}{n_{10}}} = n_{1,X_{NaCl}=1} \f$
         * * (2) Use condition \f$ T^* = T \f$ when \f$ X=0 \f$, so \f$ n_1 = n_{10} + n_{11} + n_{12} = 0 \f$. Therefore, \f$ {\color{red}{n_{12}}} = -n_{10} - n_{11} \f$ where \f$ n_{11} \f$ is given in the Tab. 4. 
         * * (3) when \f$ X = 0 \f$, \f$ n_2 = n_{20} + n_{21}\sqrt{n_{22}} = 1 \f$. Therefore, \f$ {\color{red}{n_{20}}} = 1 - n_{21}\sqrt{n_{22}}  \f$, where \f$ n_{21}, n_{22} \f$ are given in the Tab.4.
         * * (4) when \f$ X = 1 \f$, \f$ n_2 = n_{20} + n_{21}\sqrt{n_{22} + 1} + n_{23} = n_{2,X_{NaCl}=1} \f$. Therefore, \f$ {\color{red}{n_{23}}} = n_{2,X_{NaCl}=1} - n_{20} - n_{21}\sqrt{n_{22} + 1} \f$
         * 
         * \image html H2ONaCl/H2ONaCl_n1n2_Tstar_V.svg Composition and pressure-dependence of the parameters n1 (a), n2 (b) and D(T) (c) width=90%
         * 
         * @param P [bar]
         * @param X [mol/mol]
         * @param n1 [deg.C]
         * @param n2 [-]
         */
        void cH2ONaCl::n1n2_Tstar_V(const double& P, const double& X, double& n1, double& n2)
        {
            double one_minus_X = 1-X;
            double P_sqrt = sqrt(P), PP = P*P, PPP=P*PP;

            double n11 = -54.2958 - 45.7623 * exp(-0.000944785 * P);
            double n1_XNaCl = 330.47 + 0.942876 * P_sqrt + 0.0817193 * P - 2.47556E-08 * PP + 3.45052E-10 * PPP; //eq. 11
            double n10 = n1_XNaCl;                          //eq. 10 when X_NaCl=1
            double n12 = -n11 - n10;                        //eq. 10 when X_NaCl=0
            n1 = n10 + n11 * one_minus_X + n12 * one_minus_X*one_minus_X;

            double n21 = -2.6142 - 0.000239092 * P;
            double n22 = 0.0356828 + 4.37235E-06 * P + 2.0566E-09 * PP;
            double n2_XNaCl = -0.0370751 + 0.00237723 * P_sqrt + 5.42049E-05 * P + 5.84709E-09 * PP - 5.99373E-13 * PPP; //eq. 12
            double n20 = 1 - n21 * sqrt(n22);                       //eq. 10 when X_NaCl=0
            double n23 = n2_XNaCl - n20 - n21 * sqrt(1 + n22);      //eq. 10 when X_NaCl=1
            n2 = n20 + n21 * sqrt(X + n22) + n23 * X;
        }

        /**
         * @brief Calculate the deviation function \f$ D(T) \f$ in Eq.13 of \cite Driesner2007Part2.
         * 
         * \image html H2ONaCl/H2ONaCl_n1n2_Tstar_V.svg Composition and pressure-dependence of the parameters n1 (a), n2 (b) and D(T) (c) width=90%
         * 
         * @param T [deg.C]
         * @param P [bar]
         * @param X [mol/mol]
         * @return double [deg.C]
         */
        double cH2ONaCl::D_Tstar_V(const double& T, const double& P, const double& X)
        {
            double n300 = 7606640/pow(P + 472.051, 2.0);
            double n301 = -50 - 86.1446 * exp(-0.000621128 * P);
            double n302 = 294.318 * exp(-0.00566735 * P);
            double n310 = -0.0732761 * exp(-0.0023772 * P) - 5.2948E-05 * P;
            double n311 = -47.2747 + 24.3653 * exp(-0.00125533 * P);
            double n312 = -0.278529 - 0.00081381 * P;
            double n30 = n300 * (exp(n301 * X) - 1) + n302 * X; //eq. 15
            double n31 = n310 * exp(n311 * X) + n312 * X; //eq. 16
            return n30 * exp(n31 * T); //eq. 14
        }

        /**
         * @brief Calculate extrapolation density at high temperature region.
         * 
         * See eq.18 of \cite Driesner2007Part2. Valid T-P-X range: T >600deg.C, P<#P_Peak_VLH, X>0.2 mole fraction.
         * 
         * How to calculate \f$ o_3, o_4, o_5\f$ in eq. 18 ? First, need to calculate three values: \f$ V_{solution}(P=P_{max}[VLH]), \left( \frac{dV_{solution}}{d P} \right)_{P=P_{max}[VLH]}, V_{solution}(P=1000bar) \f$, where \f$ P_{max}[VLH] \f$ in unit bar.
         * 
         * * 1. \f$ {\color{blue}{\Delta V}} = V_{solution}(P=1000bar) - V_{solution}(P=P_{max}[VLH]) = \color{red}{o_4}(ln^{1000+1000} - ln^{P_{max}[VLH] + 1000}) + \color{red}{o_5}(1000 - P_{max}[VLH]) \f$
         * * 2. \f$ {\color{blue}{dVdP}} = \left(\frac{dV_{solution}}{d P} \right)_{P=P_{max}[VLH]} =  \frac{\color{red}{o_4}}{P_{max}[VLH] + 1000} + \color{red}{o_5} \f$
         * 
         * Combine the above two equations, we can get \f$ o_4 \f$ and \f$ o_5 \f$ :
         * 
         * \f{equation}
         *      o_4 = \frac{ {\color{blue}{\Delta V}} - {\color{blue}{dVdP}(1000 - P_{max}[VLH])}}{ln^{1000+1000} - ln^{P_{max}[VLH] + 1000} - \frac{1000 - P_{max}[VLH]}{P_{max}[VLH] + 1000}}
         * \f}
         * 
         * \f{equation}
         *      o_5 = {\color{blue}{dVdP}} - \frac{\color{red}{o_4}}{P_{max}[VLH] + 1000}
         * \f}
         * 
         * \f{equation}
         *      o_3 = V_{solution}(P=P_{max}[VLH]) - o_4ln^{P_{max}[VLH] + 1000} - o_5P_{max}[VLH]
         * \f}
         * 
         * @param T_C [deg.C]
         * @param P_bar [bar]
         * @param X_mol [mol/mol]
         * @return double 
         */
        double cH2ONaCl::extrapolation_V_highT(const double& T_C, const double& P_bar, const double& X_mol)
        {
            double n1, n2, Tstar;
            double dP = 1; //bar
            // extrapolation in high temperature but low pressure region
            // 1.1 calculate rho at peak of VLH zone
            n1n2_Tstar_V(P_Peak_VLH_bar, X_mol, n1, n2);
            Tstar = n1 + n2*T_C; // deg.C
            double rho_water_P390 = _Rho_water(Tstar+273.15, P_Peak_VLH);
            double V_water_P390 = m_constants_Water.molar_mass/rho_water_P390; // m^3/mol
            // printf("Tstar390: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);
            // 1.2 dRhodP at PmaxVLH
            n1n2_Tstar_V(P_Peak_VLH_bar + dP, X_mol, n1, n2);
            Tstar = n1 + n2*T_C; // deg.C
            double rho_water_P390_dP = _Rho_water(Tstar+273.15, (P_Peak_VLH_bar+dP)*1E5);
            double V_water_P390_dP = m_constants_Water.molar_mass/rho_water_P390_dP; // m^3/mol
            double dVdP = (V_water_P390_dP - V_water_P390)/dP;
            // printf("Tstar400: %.8E, dVdP: %.8E\n",Tstar,dVdP);
            // 1.3 rho at P=1000 bar
            n1n2_Tstar_V(1000, X_mol, n1, n2);
            Tstar = n1 + n2*T_C; // deg.C
            double rho_water_1000bar = _Rho_water(Tstar+273.15, 1000E5);
            double V_water_1000bar = m_constants_Water.molar_mass/rho_water_1000bar; // m^3/mol
            // printf("Tstar1000: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);

            // 1.4 delta V
            double delta_V = V_water_1000bar - V_water_P390;
            double delta_P = 1000 - P_Peak_VLH_bar;
            double P_Peak_VLH_bar_plus1000 = P_Peak_VLH_bar + 1000;
            // 2.1 calculate o4
            double o4 = (delta_V - dVdP*delta_P)/(log(2000) - log(P_Peak_VLH_bar_plus1000) - delta_P/P_Peak_VLH_bar_plus1000);
            double o5 = dVdP - o4/P_Peak_VLH_bar_plus1000;
            double o3 = V_water_P390 - o4*log(P_Peak_VLH_bar_plus1000) - o5*P_Peak_VLH_bar;
            // printf("rho390: %.4E, rho400: %.4E, rho1000: %.4E\n",rho_water_P390,rho_water_P390_dP,rho_water_1000bar);
            // std::cout<<"T: "<<T_C<<", P: "<<P_bar<<", X: "<<X_mol<<", o3: "<<o3<<", o4: "<<o4<<", o5: "<<o5<<", v_water390: "<<V_water_P390<<std::endl;
            // printf("v390: %.4E, v400: %.4E, v1000: %.4E\n",V_water_P390, V_water_P390_dP,V_water_1000bar);
            return o3 + o4*log(P_bar + 1000) + o5*P_bar; //kg/m^3
        }
        
        /**
         * @brief Calculate extrapolation density at low-T low-P region.
         * 
         * See eq.17 of \cite Driesner2007Part2. Valid pressure range: P < 15 bar.
         * 
         * @param Tstar_C [deg.C]
         * @param P_bar [bar]
         * @param X_mol [mol/mol]
         * @return double 
         */
        double cH2ONaCl::extrapolation_V_lowPlowT(const double& Tstar_C, const double& P_bar, const double& X_mol)
        {
            // printf("Tstar[C]: %.6E, P[bar]: %.6E, X[mol]: %.6E\n",Tstar_C, P_bar, X_mol);
            // 1. calculate critical density of water
            double rho_water_crit_l, rho_water_crit_v, T_water_crit;
            T_water_crit = m_Water->Boiling_T(P_bar*1E5, rho_water_crit_l, rho_water_crit_v);
            double V_water_crit_l = m_constants_Water.molar_mass/rho_water_crit_l; // m^3/mol
            // 2. calculate V(Tcrit+dT)
            double dT = -1;
            double rho_water_l_dT = _Rho_water(T_water_crit+dT, P_bar*1E5);
            double V_water_l_dT = m_constants_Water.molar_mass/rho_water_l_dT; // m^3/mol
            double dVdT = (V_water_l_dT - V_water_crit_l)/dT; // m^3/mol/degC
            // 3. calculate o0, o1, o2
            double logPbar = log10(P_bar);
            double Tstar_sqre = Tstar_C*Tstar_C;
            double o2 =2.0125e-13 + 3.29977e-15*exp(-4.31279*logPbar) - 1.17748e-13*logPbar + 7.58009e-14*logPbar*logPbar; // m^3/mol/(degC)^3. Note that the unit of o2 in Table 4 of Driesner(2007) is cm^3/mol, therefore need to multiply 1E-6 to convert it to m^3/mol
            double o1 = dVdT - 3*o2*Tstar_sqre; // m^3/mol/degC
            double o0 = V_water_crit_l - o1*(T_water_crit-273.15) - o2*pow(T_water_crit-273.15, 3.0); // m^3/mol
            // printf("Tcrit[K]: %.6E, Tcrit[C]: %.6E, rho_l_crit: %.6E\n", T_water_crit, T_water_crit-273.15, rho_water_crit_l);
            // printf("Rho_crit: %.6E, Rho_dT: %.6E\n",rho_water_crit_l, rho_water_l_dT);
            // printf("dVdT: %.6E, o0: %.6E, o1: %.6E\n",dVdT, o0, o1);
            // printf("Vextra: %.6E\n",o0 + o1*Tstar_C);
            return o0 + o1*Tstar_C + o2*Tstar_sqre*Tstar_C;
        }

        /**
         * @brief Calculate density of liquid or vapor phase.
         * 
         * The density of H2ONaCl can be calculated from water density, see Eq. 7 of \cite Driesner2007Part2.
         * \f$ V_{solution}(P,T,X) = V_{water}(T^*,P) \f$, so \f$ \frac{M_{solution}}{\rho_{solution}} = \frac{M_{water}}{\rho_{water}} \f$.
         * 
         * Therefore, \f$ \rho_{solution} = \frac{M_{solution}}{M_{H_2O}} \rho_{H_2O} \f$, where \f$ M_{solution} = M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol} \f$
         * 
         * \note Density difference calculated from IAPWS95_CoolProp and IAPS84 will cause H2ONaCl density difference up to 32 kg/m^3, test the following points:
         * 
         * 
            | T[c]         | P[bar]       | X[mole fraction] | rho[\cite Driesner2007Part2] | rho[xThermal] |
            |--------------|--------------|------------------|------------------------------|--------------|
            | 8.500000e+02 | 2.000000e+03 | 9.996030e-01     | 1.608149e+03                 | 1.602652e+03 |
            | 8.600000e+02 | 2.500000e+03 | 9.937537e-01     | 1.615147e+03                 | 1.605205e+03 |
            | 9.200000e+02 | 5.000000e+03 | 9.915095e-01     | 1.651829e+03                 | 1.619722e+03 |
            | 9.243300e+02 | 5.000000e+03 | 1.000000e+00     | 1.656120e+03                 | 1.623869e+03 |

            For example, T=920 deg.C, P = 5000 bar. The calculated \f$ T^* \f$ for density is 1279.92 deg.C. Density of water at (T=1279.92 deg.C, P=5000 bar) is 502.235 kg/m^3 (IAPWS95_CoolProp), 512.153 kg/m^3 (IAPS84). The coefficient \f$ \frac{M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol}}{M_{H_2O}} = 3.225 \f$ at this P-T-X condition, therefore the final density different will be \f$ 3.225\times(512.153-502.235)=32 \f$ (This P-T condition exceeds the valid range of IAPWS-IF97).

            \todo Estimate density difference range based on the density difference between IAPWS95_CoolProp and IAPS84.

            \todo How to deal with properties at high-T low-P region, e.g. sw.Rho_phase(1.700000e+02+273.15,7.914706e+00*1E5,0,Rhol, H2ONaCl::Liquid); see mmc4 of Driesner(2007b).

             \todo Merge Rho_phase of liquid and vapor together to calculate them simultaneously, because property of saturated vapor and liquid always needed simultaneous, it could save computing time in this way.
            \bug T=273.16, P=20E5, X=0.1, Tstar = n1 + n2*T_C + D will be negative.

        * @param T [K]
        * @param P [Pa]
        * @param X [kg/kg] salinity of liquid or salinity of vapor, not bulk salinity
        * @param Tstar [K] 
        * @param rho [\f$ kg/m^3 \f$]
        */
        void cH2ONaCl::Rho_phase(const double& T, const double& P, const double& X, double& rho, double& dRhodP, double& dRhodT,PhaseType phase)
        {
            // X==0 means the input phase doesn't exist, how to deal with that ? Actually we don't need special process for X==0, because the "correction method" already include the case of X==0
            // if(X==0)
            // {
            //     m_Water->Boiling_p(T);
            //     if (phase == Liquid)rho = m_Water->props_sat_liquid().rho;
            //     else if (phase == Vapor)rho = m_Water->props_sat_vapor().rho;
            //     return;
            // }

            double X_mol;
            Wt2Mol(X,X_mol);
            double T_C = T - 273.15;
            double P_bar = P/1E5;
            double n1, n2;
            n1n2_Tstar_V(P_bar, X_mol, n1, n2);
            double D = D_Tstar_V(T_C, P_bar, X_mol);
            double Tstar = n1 + n2*T_C ; // + D; // deg.C   T=273.16, P=20E5, X=0.1, Tstar = n1 + n2*T_C + D will be negative. do not use D at this moment
            double dRhodP_water, dRhodT_water, kappa_water, beta_water;
            double rho_water = _Rho_water(Tstar + 273.15, P, dRhodP_water, dRhodT_water, kappa_water, beta_water, phase);  //input T in unit K
            double molarVolume2MassDensity = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/m_constants_Water.molar_mass;
            rho = molarVolume2MassDensity * rho_water;
            dRhodP = molarVolume2MassDensity * dRhodP_water;
            dRhodT = molarVolume2MassDensity * dRhodT_water;

            // If calculate density of vapor, we dont't need to make extrapolation

            // printf("n1: %.6E, n2: %.6E, D: %.6E, Tstar: %.6E, rho_water: %.6E\n",n1,n2,D,Tstar,rho_water);
            // printf("vol: %.6E, vol_l: %.6E, mass_sol_l: %.6E, coeff: %.6E\n",1.0/rho_water, m_constants_Water.molar_mass/rho_water,(m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol),(m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/m_constants_Water.molar_mass);
            // Extrapolation for liquid phase in low-p low-T and low-p high-T region
            if(phase==Liquid)
            {
                if(T>=873.15 && P<P_Peak_VLH && X_mol>0.1)
                {
                    rho = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_highT(T_C,P_bar,X_mol);
                    double dT = 0.1;
                    double rho_dT = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_highT(T_C + dT,P_bar,X_mol);
                    dRhodT = (rho_dT - rho)/dT;
                    double dP_bar = 0.1;
                    double rho_dP = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_highT(T_C,P_bar + dP_bar,X_mol);
                    dRhodP = (rho_dP - rho)/(dP_bar*1E5);
                    // std::cout<<rho<<std::endl;
                }else if(P<15E5 && rho_water < m_constants_Water.rhomass_critical) // see also Line 142 of calc_rho.m
                {
                    rho = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_lowPlowT(Tstar, P_bar, X_mol);
                    double dT = 0.1;
                    double rho_dT = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_lowPlowT(Tstar +dT, P_bar, X_mol);
                    dRhodT = (rho_dT - rho)/dT;
                    double dP_bar = 0.1;
                    double rho_dP = (m_constants_Water.molar_mass*(1-X_mol) + NaCl::M*X_mol)/extrapolation_V_lowPlowT(Tstar, P_bar + dP_bar, X_mol);
                    dRhodP = (rho_dP - rho)/(dP_bar*1E5);
                }
            }
        }

        /**
         * @brief Calculate \f$q_1, q_2 \f$ for the \f$ T^*_h \f$ in Eq. 22 of \cite Driesner2007Part2 .
         * 
         * * How to calculate \f$ q_{10}, q_{12} \f$ by given \f$ q_{11} \f$ formula in Tab. 5 of \cite Driesner2007Part2.
         *   \f{align}
         *   \left. \begin{matrix}
         *   q_1(X=1) &=& q_{10} = q_{1, X_{NaCl}=1} \\ 
         *   q_1(X=0) &=& q_{10} + q_{11} + q_{12} = 0
         *   \end{matrix}  \right\} \Longrightarrow \left\{ 
         *       \begin{matrix}
         *       q_{10} = q_{1, X_{NaCl}=1} \\
         *       q_{12} = - q_{10} - q_{11}
         *       \end{matrix}
         *       \right.
         *   \f}
         * 
         * * How to calculate \f$ q_{20}, q_{23} \f$ by given \f$ q_{21}, q_{22} \f$ formula in Tab. 5 of \cite Driesner2007Part2.
        * 
        *     \f{align}
        *    \left. \begin{matrix}
        *    q_2(X=1) &=& q_{20} + q_{21}\sqrt{1+q_{22}} + q_{23} = q_{2, X_{NaCl}=1} \\ 
        *    q_2(X=0) &=& q_{20} + q_{21}\sqrt{q_{22}} = 1
        *    \end{matrix}  \right\} \Longrightarrow \left\{ 
        *        \begin{matrix}
        *        q_{20} = 1- q_{21}\sqrt{q_{22}} \\
        *        q_{23} = q_{2, X_{NaCl}=1} - q_{21}\sqrt{1+q_{22}} - q_{20}
        *        \end{matrix}
        *        \right.
        *    \f}
        * 
        * @param P [bar]
        * @param X [mol/mol]
        * @param q1 
        * @param q2 
        */
        void cH2ONaCl::q1q2_Tstar_H(const double& P_bar, const double& X_mol, double& q1, double& q2)
        {
            double one_minus_X = 1-X_mol;
            double P_sqrt = sqrt(P_bar), PP = P_bar*P_bar;

            double q11 = -32.1724 + 0.0621255*P_bar; //Tab.5
            double q21 = -1.69513 - 4.52781E-4*P_bar - 6.04279E-8 * PP;
            double q22 = 0.0612567 + 1.88082E-5 * P_bar;
            double q1_XNaCl_1 = 47.9048 - 9.36994E-3*P_bar + 6.51059E-6 * PP;           //eq. 25, X_NaCl=1
            double q2_XNaCl_1 = 0.241022 + 3.45087E-5 * P_bar - 4.28356E-9 * PP;        //eq. 26, X_NaCl=1
            double q10 = q1_XNaCl_1;                                                    //eq. 23 and eq.25, X_NaCl=1
            double q12 = -q10 - q11;
            double q20 = 1 - q21*sqrt(q22);
            double q23 = q2_XNaCl_1 - q21*sqrt(1+q22) - q20;

            q1 = q10 + q11*one_minus_X + q12*one_minus_X*one_minus_X;
            q2 = q20 + q21*sqrt(X_mol + q22) + q23*X_mol;
        }

        double cH2ONaCl::extrapolation_H_highT(const double& T_C, const double& P_bar, const double& X_mol)
        {
            double q1, q2, Tstar;
            double dP = 1; //bar
            // extrapolation in high temperature but low pressure region
            PhaseType phase = Liquid;
            // 1.1 calculate rho at peak of VLH zone
            q1q2_Tstar_H(P_Peak_VLH_bar, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
            double H_water_P390 = _H_water(Tstar+273.15, P_Peak_VLH, phase);
            double H_water_P390_vapor = _H_water(Tstar+273.15, P_Peak_VLH, phase);
            // printf("Tstar390: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);
            // 1.2 dHdP at PmaxVLH
            q1q2_Tstar_H(P_Peak_VLH_bar + dP, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
            double H_water_P390_dP = _H_water(Tstar+273.15, (P_Peak_VLH_bar+dP)*1E5, phase);
            double H_water_P390_dP_vapor = _H_water(Tstar+273.15, (P_Peak_VLH_bar+dP)*1E5, phase);
            double dHdP = (H_water_P390_dP - H_water_P390)/dP;
            // printf("Tstar400: %.8E, dVdP: %.8E\n",Tstar,dVdP);
            // 1.3 rho at P=1000 bar
            q1q2_Tstar_H(1000, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
            double H_water_1000bar = _H_water(Tstar+273.15, 1000E5, phase);
            double H_water_1000bar_vapor = _H_water(Tstar+273.15, 1000E5, phase);
            // printf("Tstar1000: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);

            // 1.4 delta H
            double delta_H = H_water_1000bar - H_water_P390;
            double delta_P = 1000 - P_Peak_VLH_bar;
            double P_Peak_VLH_bar_plus1000 = P_Peak_VLH_bar + 1000;
            // 2.1 calculate o4
            double o4 = (delta_H - dHdP*delta_P)/(log(2000) - log(P_Peak_VLH_bar_plus1000) - delta_P/P_Peak_VLH_bar_plus1000);
            double o5 = dHdP - o4/P_Peak_VLH_bar_plus1000;
            double o3 = H_water_P390 - o4*log(P_Peak_VLH_bar_plus1000) - o5*P_Peak_VLH_bar;
            // printf("rho390: %.4E, rho400: %.4E, rho1000: %.4E\n",rho_water_P390,rho_water_P390_dP,rho_water_1000bar);
            // std::cout<<"T: "<<T_C<<", P: "<<P_bar<<", X: "<<X_mol<<", o3: "<<o3<<", o4: "<<o4<<", o5: "<<o5<<", v_water390: "<<V_water_P390<<std::endl;
            // printf("v390: %.4E, v400: %.4E, v1000: %.4E\n",V_water_P390, V_water_P390_dP,V_water_1000bar);
            return o3 + o4*log(P_bar + 1000) + o5*P_bar; //J/kg
        }

        /**
         * @brief Calculate extrapolation enthalpy and specific heat at high temperature region.
         *
         * \cite Driesner2007Part2 did not provide such extrapolation formula, but it can be extended in the similar way of density, see eq.18 of \cite Driesner2007Part2. Valid T-P-X range: T >600deg.C, P<#P_Peak_VLH, X>0.2 mole fraction.
         *
         * \f{equation}
         * H = o_3 + o_4ln^{P+1000} + o_5P
         * \f}
         *
         * How to calculate \f$ o_3, o_4, o_5\f$ in eq. 18 ? First, need to calculate three values: \f$ H_{solution}(P=P_{max}[VLH]), \left( \frac{dH_{solution}}{d P} \right)_{P=P_{max}[VLH]}, H_{solution}(P=1000bar) \f$, where \f$ P_{max}[VLH] \f$ in unit bar.
         *
         * * 1. \f$ {\color{blue}{\Delta H}} = H_{solution}(P=1000bar) - H_{solution}(P=P_{max}[VLH]) = \color{red}{o_4}(ln^{1000+1000} - ln^{P_{max}[VLH] + 1000}) + \color{red}{o_5}(1000 - P_{max}[VLH]) \f$
         * * 2. \f$ {\color{blue}{dHdP}} = \left(\frac{dH_{solution}}{d P} \right)_{P=P_{max}[VLH]} =  \frac{\color{red}{o_4}}{P_{max}[VLH] + 1000} + \color{red}{o_5} \f$
         *
         * Combine the above two equations, we can get \f$ o_4 \f$ and \f$ o_5 \f$ :
         *
         * \f{equation}
         *      o_4 = \frac{ {\color{blue}{\Delta H}} - {\color{blue}{dHdP}(1000 - P_{max}[VLH])}}{ln^{1000+1000} - ln^{P_{max}[VLH] + 1000} - \frac{1000 - P_{max}[VLH]}{P_{max}[VLH] + 1000}}
         * \f}
         *
         * \f{equation}
         *      o_5 = {\color{blue}{dHdP}} - \frac{\color{red}{o_4}}{P_{max}[VLH] + 1000}
         * \f}
         *
         * \f{equation}
         *      o_3 = H_{solution}(P=P_{max}[VLH]) - o_4ln^{P_{max}[VLH] + 1000} - o_5P_{max}[VLH]
         * \f}
         *
         * Note that \f$o_3, o_4\f$ and \f$o_5\f$ are function of temperature. Therefore, the isobaric specific heat capacity can be expressed
         * \f{equation}
         * C_p = \left( \frac{\partial H}{\partial T} \right)_P = \frac{\partial o_3}{\partial T} + ln^{P+1000}\frac{\partial o_4}{\partial T} + P \frac{\partial o_5}{\partial T}
         * \f}
         *
         * @param T_C
         * @param P_bar
         * @param X_mol
         * @param H Specific enthalpy
         * @param Cp Isobaric specific heat capacity
         */
        void cH2ONaCl::extrapolation_H_Cp_highT(const double& T_C, const double& P_bar, const double& X_mol, double& H, double& Cp)
        {
            double q1, q2, Tstar;
            double dP = 1; //bar
            // extrapolation in high temperature but low pressure region
            PhaseType phase = Liquid;
            // 1.1 calculate rho at peak of VLH zone
            q1q2_Tstar_H(P_Peak_VLH_bar, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
            double H_water_P390, Cp_water_P390;
            _H_Cp_water(Tstar+273.15, P_Peak_VLH, H_water_P390, Cp_water_P390, phase);
            Cp_water_P390 = Cp_water_P390*q2;
            // printf("Tstar390: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);
            // 1.2 dHdP at PmaxVLH
            q1q2_Tstar_H(P_Peak_VLH_bar + dP, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
//            Cp = Cp_water * q2;
            double H_water_P390_dP,Cp_water_P390_dP;
            _H_Cp_water(Tstar+273.15, (P_Peak_VLH_bar+dP)*1E5, H_water_P390_dP, Cp_water_P390_dP, phase);
            Cp_water_P390_dP = Cp_water_P390_dP*q2;
            double H_water_P390_dP_vapor, Cp_water_P390_dP_vapor;
            _H_Cp_water(Tstar+273.15, (P_Peak_VLH_bar+dP)*1E5, H_water_P390_dP_vapor, Cp_water_P390_dP_vapor, phase);
            Cp_water_P390_dP_vapor = Cp_water_P390_dP_vapor*q2;
            double dHdP = (H_water_P390_dP - H_water_P390)/dP;
            double dCpdP = (Cp_water_P390_dP - Cp_water_P390)/dP;
            // printf("Tstar400: %.8E, dVdP: %.8E\n",Tstar,dVdP);
            // 1.3 rho at P=1000 bar
            q1q2_Tstar_H(1000, X_mol, q1, q2);
            Tstar = q1 + q2*T_C; // deg.C
            double H_water_1000bar, Cp_water_1000bar;
            _H_Cp_water(Tstar+273.15, 1000E5, H_water_1000bar, Cp_water_1000bar, phase);
            Cp_water_1000bar = Cp_water_1000bar*q2;
            double H_water_1000bar_vapor, Cp_water_1000bar_vapor;
            _H_Cp_water(Tstar+273.15, 1000E5, H_water_1000bar_vapor, Cp_water_1000bar_vapor, phase);
            Cp_water_1000bar_vapor = Cp_water_1000bar_vapor*q2;
            // printf("Tstar1000: %.8E, P_Peak_VLH: %.8E\n",Tstar,P_Peak_VLH_bar);

            // 1.4 delta H
            double delta_H = H_water_1000bar - H_water_P390;
            double delta_Cp = Cp_water_1000bar - Cp_water_P390;
            double delta_P = 1000 - P_Peak_VLH_bar;
            double P_Peak_VLH_bar_plus1000 = P_Peak_VLH_bar + 1000;
            double log_P_Peak_VLH_bar_plus1000 = log(P_Peak_VLH_bar_plus1000);
            // 2.1 calculate o4
            double denominator_o4 = (log(2000) - log_P_Peak_VLH_bar_plus1000 - delta_P/P_Peak_VLH_bar_plus1000);
            double o4 = (delta_H - dHdP*delta_P)/denominator_o4;
            double do4dT = (delta_Cp - dCpdP*delta_P)/denominator_o4;
            double o5 = dHdP - o4/P_Peak_VLH_bar_plus1000;
            double do5dT = dCpdP - do4dT/P_Peak_VLH_bar_plus1000;
            double o3 = H_water_P390 - o4*log_P_Peak_VLH_bar_plus1000 - o5*P_Peak_VLH_bar;
            double do3dT = Cp_water_P390 - do4dT*log_P_Peak_VLH_bar_plus1000 - do5dT*P_Peak_VLH_bar;
            // printf("rho390: %.4E, rho400: %.4E, rho1000: %.4E\n",rho_water_P390,rho_water_P390_dP,rho_water_1000bar);
            // std::cout<<"T: "<<T_C<<", P: "<<P_bar<<", X: "<<X_mol<<", o3: "<<o3<<", o4: "<<o4<<", o5: "<<o5<<", v_water390: "<<V_water_P390<<std::endl;
            // printf("v390: %.4E, v400: %.4E, v1000: %.4E\n",V_water_P390, V_water_P390_dP,V_water_1000bar);
            H = o3 + o4*log(P_bar + 1000) + o5*P_bar; //J/kg
            Cp = do3dT + do4dT*log(P_bar + 1000) + do5dT * P_bar;
        }

        /**
         * @brief Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl.
         *
         * For high-temperature(T>800 deg.C) and low pressure (P<P_vlh), implement the similar extrapolation method with that in \like Rho_phase \endlink.
         *
         *
         * \warning The low-T low-P(P<1bar) extrapolation doesn't work for VLH surface, actually we don't need that because this part of VLH surface already below the minimum valid pressure (1bar).
         *
        * @param T [K]
        * @param P [Pa]
        * @param X [kg/kg]
        * @param H [J/kg]
        * @param phase [Liquid, Vapor]
        */
        void cH2ONaCl::H_phase(const double& T, const double& P, const double& X, double& H, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(X==0)
            {
                m_Water->Boiling_T(P, props);
                if (phase == Liquid)H = props.H_l;//m_Water->props_sat_liquid().H;
                else if (phase == Vapor)H = props.H_v;//m_Water->props_sat_vapor().H;
                return;
            }
            double X_mol;
            Wt2Mol(X, X_mol);
            double T_C = T - 273.15;
            double P_bar = P/1E5;
            double q1, q2;
            q1q2_Tstar_H(P_bar, X_mol, q1, q2);
            double Tstar = q1 + q2*T_C;
            H = _H_water(Tstar + 273.15, P,phase);
            // ========== extrapolation for high temperature low pressure region, only for liquid phase =======
            // this is not mentioned in Driesner(2007), the coefficients from the matlab code calc_enthalpy.m (L.113-176)
            if (phase==Liquid)
            {
                if (P<=P_Peak_VLH && T>=873.15 && phase==Liquid)
                {
                    H = extrapolation_H_highT(T_C,P_bar,X_mol);
                }
                // low pressure low T extrapolation doesn't work
                // else if (phase==Liquid && T<373 && P< m_Water->Boiling_p(Tstar+273.15)) //extrapolation for low temperature region: the condition used in the matlab code is (h_l > 2.086e6 || std::isnan(h_l))  &&  P_star_l < P_crit  &&  T_in < 375
                else if((H>2.086E6 || isnan(H)) && (P<m_constants_Water.p_critical) && T <=m_constants_Water.T_critical)
                {
                    m_Water->Boiling_T(P, props);
                    double hl_crit = props.H_l;//m_Water->props_sat_liquid().H;
                    double T_boiling_h2o = props.T;//m_Water->props_sat_liquid().T;
                    m_Water->UpdateState_TPX(props, T_boiling_h2o -1, P);
                    double h_l_minus_one=props.Rho;//m_Water->hmass();
                    double dh_ldT = (hl_crit - h_l_minus_one);
                    double o1 = dh_ldT ;
                    double o0 = hl_crit - o1 * (T_boiling_h2o - 273.15) ;
                    H = o0 + o1 *  Tstar;
                }
            }
        }
        /**
         * @brief Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl.
         *
         * For high-temperature(T>800 deg.C) and low pressure (P<P_vlh), implement the similar extrapolation method with that in \like Rho_phase \endlink.
         *
         * \f{align}
         * \begin{cases}
         * H_{solution} (T,P,X) &= H_{H_2O}(T^*_H, P)\\
         * H_{extrapol}^{lowT} (T,P,X) &= o_0 + o_1 T^*_V + o_2 T^{*3}_{H},& ~~  (P<P_{H_2O}^{crit}, T<=T_{H_2O}^{crit})\\
         * H_{extrapol}^{highT} (T,P,X) &= o_3 + o_4 ln^{P + 1000} + o_5 P,& ~~ (P<P_{vlh}^{peak}, T>600 ^{\circ}C)
         * \end{cases}
         * \f}
         *
         * @param T
         * @param P
         * @param X
         * @param H
         * @param Cp
         * @param phase
         */
        void cH2ONaCl::H_phase(const double& T, const double& P, const double& X, double& H, double& Cp, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(X==0)
            {
                m_Water->Boiling_T(P, props);
                if (phase == Liquid)H = props.H_l;//m_Water->props_sat_liquid().H;
                else if (phase == Vapor)H = props.H_v;//m_Water->props_sat_vapor().H;
                return;
            }
            double X_mol;
            Wt2Mol(X, X_mol);
            double T_C = T - 273.15;
            double P_bar = P/1E5;
            double q1, q2;
            q1q2_Tstar_H(P_bar, X_mol, q1, q2);
            double Tstar = q1 + q2*T_C;
            double Cp_water = 0;
            _H_Cp_water(Tstar + 273.15, P, H, Cp_water,phase);
            Cp = Cp_water * q2;
            // ========== extrapolation for high temperature low pressure region, only for liquid phase =======
            // this is not mentioned in Driesner(2007), the coefficients from the matlab code calc_enthalpy.m (L.113-176)
            if (phase==Liquid)
            {
                if (P<=P_Peak_VLH && T>=873.15)
                {
                    extrapolation_H_Cp_highT(T_C,P_bar,X_mol, H, Cp);
                }
                    // low pressure low T extrapolation doesn't work
                    // else if (phase==Liquid && T<373 && P< m_Water->Boiling_p(Tstar+273.15)) //extrapolation for low temperature region: the condition used in the matlab code is (h_l > 2.086e6 || std::isnan(h_l))  &&  P_star_l < P_crit  &&  T_in < 375
                else if((H>2.086E6 || isnan(H)) && (P<m_constants_Water.p_critical) && T <=m_constants_Water.T_critical)
                {
                    m_Water->Boiling_T(P, props);
                    double hl_crit = props.H_l;//m_Water->props_sat_liquid().H;
                    double T_boiling_h2o = props.T;//m_Water->props_sat_liquid().T;
                    m_Water->UpdateState_TPX(props, T_boiling_h2o -1, P);
                    double h_l_minus_one=props.Rho;//m_Water->hmass();
                    double dh_ldT = (hl_crit - h_l_minus_one);
                    double o1 = dh_ldT ;
                    double o0 = hl_crit - o1 * (T_boiling_h2o - 273.15) ;
                    H = o0 + o1 *  Tstar;
                    Cp = o1;
                }
            }
        }
        /**
         * Calculate dynamic viscosity of H2O-NaCl for given T,P,X. See \cite klyukin2017revised equation (3,4,5,6).
         *
         * \warning For low temperature case, the T_star will be negative, e.g. T = 1 + 273.15, P = 200E5, X=0.1; see also Fig. 7a of \cite klyukin2017revised. Use minimum T ?
         *
         * @param T [K]
         * @param P [Pa]
         * @param X [kg/kg]
         * @param Mu [Pa s]
         * @param phase [Liquid|Vapor]
         */
        void cH2ONaCl::Mu_phase(const double& T, const double& P, const double& X, double& Mu, PhaseType phase)
        {
            ThermodynamicProperties props;
            if(X==0)
            {
                m_Water->Boiling_T(P, props);
                if (phase == Liquid)Mu = props.Mu_l;
                else if (phase == Vapor)Mu = props.Mu_v;
                return;
            }
            double T_C = T - 273.15;
            double e1 = m_coeffs_mu.a1 * pow(X,m_coeffs_mu.a2);
            double e2 = 1.0 - m_coeffs_mu.b1 * pow(T_C,m_coeffs_mu.b2) - m_coeffs_mu.b3 * pow(X,m_coeffs_mu.a2) * pow(T_C,m_coeffs_mu.b2);
            double T_star = e1 + e2 * T_C;
            Mu = _Mu_water(std::max(T_star + 273.15, 273.16), P, phase); // For low temperature and high salinity region, the T_star is negative, just use the minimum T instead. (Added 0.01 to avoid rounding errors leading to negative T)
        }

        void cH2ONaCl::P_X_Critical(std::vector<double> T, std::vector<double> &P_crit, std::vector<double> &X_crit) {
            P_crit.clear(); X_crit.clear();
            P_crit.resize(T.size());
            X_crit.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)
            {
                P_X_Critical(T[i], P_crit[i], X_crit[i]);
            }
        }

        void cH2ONaCl::P_Critical(std::vector<double> T, std::vector<double> &P_crit) {
            P_crit.clear();
            P_crit.resize(T.size());
            for(size_t i =0; i<T.size();i++)P_Critical(T[i], P_crit[i]);
        }

        void cH2ONaCl::T_Critical(std::vector<double> P, std::vector<double> &T_crit) {
            T_crit.clear();
            T_crit.resize(P.size());
            for(size_t i =0; i<P.size();i++)T_Critical(P[i], T_crit[i]);
        }

        void cH2ONaCl::T_X_Critical(std::vector<double> P, std::vector<double> &T_crit, std::vector<double> &X_crit) {
            T_crit.clear(); X_crit.clear();
            T_crit.resize(P.size());
            X_crit.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)T_X_Critical(P[i],T_crit[i],X_crit[i]);
        }

        void cH2ONaCl::X_HaliteLiquidus(std::vector<double> T, std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){ res[i] = X_HaliteLiquidus(T[i],P[i]);}
        }

        void cH2ONaCl::P_VLH(std::vector<double> T, std::vector<double> &P) {
            P.clear();
            P.resize(T.size());
            for(size_t i =0; i<T.size();i++)P_VLH(T[i], P[i]);
        }

        void cH2ONaCl::X_VH(std::vector<double> T, std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){ res[i] = X_VH(T[i],P[i]); }
        }

        void cH2ONaCl::X_VLH(std::vector<double> T, std::vector<double> P, std::vector<double> &X_L,
                             std::vector<double> &X_V) {
            X_L.clear();                X_V.clear();
            X_L.resize(T.size());       X_V.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){ X_VLH(T[i],P[i], X_L[i], X_V[i]);}
        }

        void cH2ONaCl::X_VL(std::vector<double> T, std::vector<double> P, std::vector<double>& X_L, std::vector<double>& X_V)
        {
            X_L.clear();                X_V.clear();
            X_L.resize(T.size());       X_V.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){ X_VL(T[i],P[i], X_L[i], X_V[i]);}
        }

        void cH2ONaCl::XL_VL(std::vector<double> T, std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){res[i] = XL_VL(T[i],P[i]);}
        }

        void cH2ONaCl::XV_VL(std::vector<double> T, std::vector<double> P, std::vector<double> &res) {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++){res[i] = XV_VL(T[i],P[i]);}
        }

        std::vector<int> cH2ONaCl::getPhaseRegion_TPX(std::vector<double> T, std::vector<double> P, std::vector<double> X, std::vector<double>& X_V, std::vector<double>& X_L)
        {
            std::vector<int> phase_region(T.size());
            X_L.clear(); X_V.clear();
            X_L.resize(T.size());
            X_V.resize(T.size());
            PhaseRegion phase_region_;
            for (size_t i = 0; i < T.size(); i++)
            {
                getPhaseRegion_TPX(T[i], P[i], X[i], phase_region_, X_V[i], X_L[i]);
                phase_region[i] = phase_region_;
            }
            return phase_region;
        }

        std::vector<double> cH2ONaCl::Rho_phase(std::vector<double> T, std::vector<double> P, std::vector<double> X, PhaseType phase)
        {
            std::vector<double> rho(T.size());
            double dRhodP, dRhodT;
            for (size_t i = 0; i < T.size(); i++)
            {
                 Rho_phase(T[i], P[i], X[i], rho[i], dRhodP, dRhodT, phase);
            }
            return rho;
        }

        std::vector<double> cH2ONaCl::H_phase(std::vector<double> T, std::vector<double> P, std::vector<double> X, PhaseType phase)
        {
            std::vector<double> H(T.size());
            for (size_t i = 0; i < T.size(); i++)
            {
                H_phase(T[i], P[i], X[i], H[i], phase);
            }
            return H;
        }

        /**
         * Create phase boundary of halite liquidus in a deformed linear "structured" mesh.
         * @param fmt_out
         * @param p_max
         * @param dT
         * @param dp
         * @return
         */
        DeformLinearMesh cH2ONaCl::PhaseBoundary_HaliteLiquidus_DeformLinear(double p_max, double dT, double dp)
        {
            double T_min = Tmin(), p_min = pmin();
            int np = int((p_max - p_min)/dp), nT = int((m_NaCl->Melting_T(p_max) -T_min)/dT);
            std::vector<double> p_hm = linspace(p_min, p_max, np); // Use CoolProp's numeric function directly
            // 1. calculate halite melting temperature as the maximum temperature bound
            std::vector<double> T_hm;
            m_NaCl->Melting_T(p_hm, T_hm);
            // 2. calculate pressure of VLH
            std::vector<double> T_vlh = linspace(T_min,m_NaCl->Melting_T(p_min),nT);
            std::vector<double> p_vlh;
            P_VLH(T_vlh, p_vlh);
            //make sure the p_vlh is valid in the low temperature region: T < T_MIN_VLH = 325.93697348 K
            for (int i = 0; i < nT; ++i) {
                p_vlh[i] = std::max(p_vlh[i], p_min);
            }
            // 3. construct T-p mesh of the halite liquidus surface
            DeformLinearMesh pb_mesh;
            pb_mesh.T.resize(np);
            pb_mesh.p.resize(np);
            pb_mesh.X.resize(np);
            // T
            for (int ip = 0; ip < np; ++ip) {
                pb_mesh.X[ip].resize(nT);
                pb_mesh.p[ip].resize(nT);
                pb_mesh.T[ip] = linspace(T_min, T_hm[ip], nT);
            }
            // p
            for (int iT = 0; iT < nT; ++iT) {
                std::vector<double>p_tmp = linspace(p_vlh[iT], p_max, np);
                for (int ip = 0; ip < np; ++ip) {
                    pb_mesh.p[ip][iT] = std::max(p_tmp[ip],pmin());
                }
            }
            // X
            for (int ip = 0; ip < np; ++ip) {
                for (int iT = 0; iT < nT; ++iT) {
                    pb_mesh.X[ip][iT] = X_HaliteLiquidus(pb_mesh.T[ip][iT], pb_mesh.p[ip][iT]);
                }
            }
            return pb_mesh;
        }

        /**
         * Calculate and return the phase boundary triangle mesh.
         *
         * @param fmt_out [txt|vtk]
         * @return
         */
        TriMesh cH2ONaCl::PhaseBoundary_HaliteLiquidus(std::string fmt_out,double p_max, double dT, double dp) {
            double T_min = Tmin(), p_min = pmin();
            int np = int((p_max - p_min)/dp), nT = int((m_NaCl->Melting_T(p_max) -T_min)/dT);
            std::vector<double> p_hm = linspace(p_min, p_max, np); // Use CoolProp's numeric function directly
            // 1. calculate halite melting temperature as the maximum temperature bound
            std::vector<double> T_hm;
            m_NaCl->Melting_T(p_hm, T_hm);
            // 2. calculate pressure of VLH
            std::vector<double> T_vlh = linspace(T_min,m_NaCl->Melting_T(p_min),nT);
            std::vector<double> p_vlh;
            P_VLH(T_vlh, p_vlh);
            //make sure the p_vlh is valid in the low temperature region: T < T_MIN_VLH = 325.93697348 K
            for (int i = 0; i < nT; ++i) {
                p_vlh[i] = std::max(p_vlh[i], p_min);
            }
            // 3. construct T-p polygon of the halite liquidus surface
            std::vector<double> x_poly, y_poly;
            for (int i = 0; i < T_vlh.size(); ++i) {
                x_poly.push_back(T_vlh[i]);
                y_poly.push_back(p_vlh[i]);
            }
            for (int i = 0; i < T_hm.size(); ++i) {
                x_poly.push_back(T_hm[i]);
                y_poly.push_back(p_hm[i]);
            }
            x_poly.push_back(T_min);
            y_poly.push_back(p_max);
            // 4. pass polygon data to Triangle
            double pointInMesh[2]={T_min+1, p_max-1E5};
            TriMesh triMesh_pb = Triangulation(x_poly, y_poly,pointInMesh[0], pointInMesh[1], dT, dp);
            // 5. calculate composition of (T,p) points in trimesh
            X_HaliteLiquidus(triMesh_pb.x, triMesh_pb.y, triMesh_pb.z);
            //6. write to file if fmt_out is valid
            if(fmt_out=="txt")
            {
                writeTriMesh2Txt(triMesh_pb);
            }
            return triMesh_pb;
        }

        TriMesh
        cH2ONaCl::Triangulation(const std::vector<double> &x_poly, const std::vector<double> &y_poly, double xIn,
                                double yIn, double dx, double dy) {
            TriMesh mesh;
            double pointInMesh[2]={xIn,yIn};
            double dxdy[2] = {dx, dy};
            cxThermal::Triangulation(x_poly, y_poly,pointInMesh, dxdy, mesh);
            return mesh;
        }

        /**
        * Calculate polygons of phase region at constant P.
        * @param P0 [Pa]
        * @return
        */
        PhaseRegion_Slice cH2ONaCl::Slice_constP(const double P0, size_t nPoints)
        {
            PhaseRegion_Slice slice;
            double T_crit, X_crit;
            T_X_Critical(P0, T_crit, X_crit);
            // first check if cross VLH zone
            if (P0<P_Peak_VLH)
            {
                double T1, T2;
                T_VLH_P0(P0,T1, T2);
                if(P0<m_Water->p_critical())
                {
                    double T_boil_water=m_Water->Boiling_T(P0);
                    ThermodynamicProperties prop_boil=m_Water->Boiling_p_props(T_boil_water);
                    T_crit = T_boil_water;
                    double H_boil_l = prop_boil.H_l;
                    double H_boil_v = prop_boil.H_v;
                    slice.points["Two phase water: L"].push_back(Point_slice(T_boil_water, P0, 0, H_boil_l, "Two phase water: L", "o", COLOR_mfc_twoPhaseWater, COLOR_mec_twoPhaseWater));
                    slice.points["Two phase water: V"].push_back(Point_slice(T_boil_water, P0, 0, H_boil_v, "Two phase water: V", "o", COLOR_mfc_twoPhaseWater, COLOR_mec_twoPhaseWater));
                    std::vector<double> T,P,X,H;
                    slice.lines["Two phase water"].push_back(Line_slice({T_boil_water, T_boil_water}, {P0, P0}, {0,0}, {H_boil_l, H_boil_v}, "Two phase water", COLOR_lc_twoPhaseWater));
                }
                else
                {
                    slice.points["Critical point"].push_back(Point_slice(T_crit, P0, X_crit, H_phase(T_crit, P0, X_crit, Vapor), "Critical point", "o", COLOR_mfc_CriticalPoint, COLOR_mec_CriticalPoint));
                }
                // VH
                std::vector<double> T_vh = xThermal::cxThermal::linspace(T1, T2, nPoints);
                std::vector<double> P_vh = xThermal::cxThermal::linspace(P0, P0, nPoints);
                std::vector<double> X_vh;
                X_VH(T_vh, P_vh, X_vh);
                std::vector<double> H_vh = H_phase(T_vh, P_vh, X_vh, Vapor);
                slice.lines["V+H"].push_back(Line_slice(T_vh, P_vh, X_vh, H_vh, "V+H", COLOR_lc_VH));
                // VL, LH
                std::vector<double>T, P, X, H;
                // -- [T2, Tmax]
                T = xThermal::cxThermal::linspace(T2, Tmax(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                    // vapor branch
                XV_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:V"].push_back(Line_slice(T,P,X,H, "V+L:V", COLOR_lc_VL_V));
                    // liquid branch
                XL_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:L"].push_back(Line_slice(T,P,X,H, "V+L:L", COLOR_lc_VL_L));
                    // halite liquidus
                T = xThermal::cxThermal::linspace(T2, m_NaCl->Melting_T(P0), nPoints);
                X_HaliteLiquidus(T,P,X);
                H = H_phase(T,P,X,Liquid);
                slice.lines["L+H"].push_back(Line_slice(T,P,X,H,"L+H",COLOR_darkgray));
                // -- [Tcrit, T1]
                T = xThermal::cxThermal::linspace(T_crit, T1, nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                    // vapor branch
                XV_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                    //append
                slice.lines["V+L:V"].push_back(Line_slice(T,P,X,H, "V+L:V", COLOR_lc_VL_V));
                    // liquid branch
                XL_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:L"].push_back(Line_slice(T,P,X,H, "V+L:L", COLOR_lc_VL_L));
                    // halite liquidus
                T = xThermal::cxThermal::linspace(Tmin(), T1, nPoints);
                X_HaliteLiquidus(T,P,X);
                H = H_phase(T,P,X,Liquid);
                slice.lines["L+H"].push_back(Line_slice(T,P,X,H,"L+H",COLOR_darkgray));
                // Enthalpy boundary at T=Tmax
                T = xThermal::cxThermal::linspace(Tmax(), Tmax(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                // ======= option 1: linear change H =================
                X = xThermal::cxThermal::linspace(slice.lines["V+L:V"][0].X[nPoints-1], slice.lines["V+L:L"][0].X[nPoints-1], nPoints);
                H = xThermal::cxThermal::linspace(slice.lines["V+L:V"][0].H[nPoints-1], slice.lines["V+L:L"][0].H[nPoints-1], nPoints);
                // ======= option 2: calculate bulk enthalpy =========
                // X = xThermal::cxThermal::linspace(slice.lines["V+L:V"][0].X[nPoints-1], slice.lines["V+L:L"][0].X[nPoints-1], nPoints, true);
                // double Hl, Hv, Rhol, Rhov, Rho, Sl, Sv, Xl=X[nPoints-1], Xv=X[0];
                // Hv = H_phase(Tmax(), P0, Xv, Vapor);
                // Hl = H_phase(Tmax(), P0, Xl, Liquid);
                // Rhov = Rho_phase(Tmax(), P0, Xv, Vapor);
                // Rhol = Rho_phase(Tmax(), P0, Xl, Liquid);
                // for (int ii = 0; ii < nPoints; ++ii) {
                //     Sl = Sl_VL(Xl, Xv, Rhol, Rhov, X[ii]);
                //     Sv = 1.0 - Sl;
                //     Rho = Rhol*Sl + Rhov*Sv;
                //     H[ii] = (Hl*Sl*Rhol + Hv*Sv*Rhov)/Rho;
                // }
                //=================================================
                slice.lines["Tmax"].push_back(Line_slice(T,P,X,H,"T=Tmax", COLOR_k, "dashed", 0.5));
                X = xThermal::cxThermal::linspace(slice.lines["V+L:L"][0].X[nPoints-1], 1.0, nPoints);
                H = H_phase(T,P,X,Liquid);
                slice.lines["Tmax"].push_back(Line_slice(T,P,X,H,"T=Tmax", COLOR_k, "dashed", 0.5));
                // VLH
                T = {T2, T1}, P={P0, P0};
                std::vector<double> XL_VLH, XV_VLH;
                X_VLH(T,P,XL_VLH, XV_VLH);
                std::vector<double> HL_VLH = H_phase(T,P, XL_VLH, Liquid);
                std::vector<double> HV_VLH = H_phase(T,P, XV_VLH, Vapor);
                std::vector<double> HH_VLH;
                m_NaCl->H_Solid(T,P, HH_VLH);
                std::vector<COLOR > ecs_VLH={COLOR_red, COLOR_blue};
                std::vector<COLOR > fcs_VLH={COLOR_VLH_highT, COLOR_VLH_lowT};
                for (int i = 0; i < 2; ++i) {
                    // points
                    slice.points["V+L+H: V"].push_back(Point_slice(T[i], P[i], XV_VLH[i], HV_VLH[i], "V+L+H: V", "o", COLOR_lime, COLOR_k));
                    slice.points["V+L+H: L"].push_back(Point_slice(T[i], P[i], XL_VLH[i], HL_VLH[i], "V+L+H: L", "o", COLOR_yellow, COLOR_k));
                    slice.points["V+L+H: H"].push_back(Point_slice(T[i], P[i], 1.0, HH_VLH[i], "V+L+H: H", "o", COLOR_blue, COLOR_k));
                    // edges
                    double Hl, Hv, Hh, Rhol, Rhov, Rhoh, Rho, Sl, Sv, Sh, Xl=XL_VLH[i], Xv=XV_VLH[i];
                    Hv = H_phase(T[i], P0, Xv, Vapor);
                    Hl = H_phase(T[i], P0, Xl, Liquid);
                    Rhov = Rho_phase(Tmax(), P0, Xv, Vapor);
                    Rhol = Rho_phase(Tmax(), P0, Xl, Liquid);
                    Rhoh = m_NaCl->Rho_Solid(T[i], P0);
                    Hh = m_NaCl->H_Solid(T[i], P0);
                    X = xThermal::cxThermal::linspace(XV_VLH[i], 0.1, (int)(nPoints/2), true);
                    std::vector<double> X_l = xThermal::cxThermal::linspace(0.2, XL_VLH[i], nPoints-(int)(nPoints/2));
                    X.insert(X.end(), X_l.begin(), X_l.end());
                    // X = xThermal::cxThermal::linspace(XV_VLH[i], XL_VLH[i], nPoints, false);
                    H=xThermal::cxThermal::linspace(HV_VLH[i], HL_VLH[i], nPoints);
                    for (int ii = 0; ii < nPoints; ++ii) {
                        Sl = Sl_VL(Xl, Xv, Rhol, Rhov, X[ii]);
                        Sv = 1.0 - Sl;
                        Rho = Rhol*Sl + Rhov*Sv;
                        H[ii] = (Hl*Sl*Rhol + Hv*Sv*Rhov)/Rho;
                    }
                    slice.lines["V+L+H: V->L"].push_back(Line_slice(xThermal::cxThermal::linspace(T[i], T[i], nPoints),
                                                              xThermal::cxThermal::linspace(P[i], P[i], nPoints),
                                                              X,
                                                              H,
                                                              "V+L+H: V->L", COLOR_lc_VLH));
                    slice.lines["V+L+H: L->H"].push_back(Line_slice(xThermal::cxThermal::linspace(T[i], T[i], nPoints),
                                                                    xThermal::cxThermal::linspace(P[i], P[i], nPoints),
                                                                    xThermal::cxThermal::linspace(XL_VLH[i], 1.0, nPoints),
                                                                    xThermal::cxThermal::linspace(HL_VLH[i], HH_VLH[i], nPoints),
                                                                    "V+L+H: L->H", COLOR_lc_VLH));
                    // X = xThermal::cxThermal::linspace(1.0, 0.1, (int)(nPoints/2));
                    // X_l = xThermal::cxThermal::linspace(0.1, XV_VLH[i], nPoints-(int)(nPoints/2), true);
                    // X.insert(X.end(), X_l.begin(), X_l.end());
                    X = xThermal::cxThermal::linspace(1.0, XV_VLH[i], nPoints, true);
                    H=xThermal::cxThermal::linspace(HH_VLH[i], HV_VLH[i], nPoints);
                    for (int ii = 0; ii < nPoints; ++ii) {
                        Sh = Saturation_Phase1(1.0, Xv, Rhoh, Rhov, X[ii]);
                        Sv = 1.0 - Sh;
                        Rho = Rhov*Sv + Rhoh*Sh;
                        H[ii] = (Hv*Sv*Rhov + Hh*Sh*Rhoh)/Rho;
                    }
                    slice.lines["V+L+H: H->V"].push_back(Line_slice(xThermal::cxThermal::linspace(T[i], T[i], nPoints),
                                                                    xThermal::cxThermal::linspace(P[i], P[i], nPoints),
                                                                    X,
                                                                    H,
                                                                    "V+L+H: H->V", COLOR_lc_VLH));
                    // polygon
                    char str_name_vlh[40];
                    sprintf(str_name_vlh,"V+L+H: T=%.3f$^{\\circ}$C", T[i]-273.15);
                    Polygon_slice polygon_vlh(slice.lines["V+L+H: V->L"][i].T, slice.lines["V+L+H: V->L"][i].P, slice.lines["V+L+H: V->L"][i].X, slice.lines["V+L+H: V->L"][i].H,string(str_name_vlh),fcs_VLH[i], ecs_VLH[i]);
                    polygon_vlh.T.insert(polygon_vlh.T.end(), slice.lines["V+L+H: L->H"][i].T.begin(), slice.lines["V+L+H: L->H"][i].T.end());
                    polygon_vlh.T.insert(polygon_vlh.T.end(), slice.lines["V+L+H: H->V"][i].T.begin(), slice.lines["V+L+H: H->V"][i].T.end());

                    polygon_vlh.P.insert(polygon_vlh.P.end(), slice.lines["V+L+H: L->H"][i].P.begin(), slice.lines["V+L+H: L->H"][i].P.end());
                    polygon_vlh.P.insert(polygon_vlh.P.end(), slice.lines["V+L+H: H->V"][i].P.begin(), slice.lines["V+L+H: H->V"][i].P.end());

                    polygon_vlh.X.insert(polygon_vlh.X.end(), slice.lines["V+L+H: L->H"][i].X.begin(), slice.lines["V+L+H: L->H"][i].X.end());
                    polygon_vlh.X.insert(polygon_vlh.X.end(), slice.lines["V+L+H: H->V"][i].X.begin(), slice.lines["V+L+H: H->V"][i].X.end());

                    polygon_vlh.H.insert(polygon_vlh.H.end(), slice.lines["V+L+H: L->H"][i].H.begin(), slice.lines["V+L+H: L->H"][i].H.end());
                    polygon_vlh.H.insert(polygon_vlh.H.end(), slice.lines["V+L+H: H->V"][i].H.begin(), slice.lines["V+L+H: H->V"][i].H.end());
                    slice.regions[string(str_name_vlh)].push_back(polygon_vlh);
                }
                // region of VL: part1
                Polygon_slice region1_VL(slice.lines["V+L:V"][0].T,slice.lines["V+L:V"][0].P,slice.lines["V+L:V"][0].X,slice.lines["V+L:V"][0].H,"V+L",COLOR_VL,COLOR_gainsboro);
                region1_VL.T.insert(region1_VL.T.end(), slice.lines["Tmax"][0].T.begin()+1, slice.lines["Tmax"][0].T.end());
                region1_VL.T.insert(region1_VL.T.end(), slice.lines["V+L:L"][0].T.rbegin(), slice.lines["V+L:L"][0].T.rend());
                region1_VL.T.insert(region1_VL.T.end(), slice.lines["V+L+H: V->L"][0].T.rbegin(), slice.lines["V+L+H: V->L"][0].T.rend());

                region1_VL.P.insert(region1_VL.P.end(), slice.lines["Tmax"][0].P.begin()+1, slice.lines["Tmax"][0].P.end());
                region1_VL.P.insert(region1_VL.P.end(), slice.lines["V+L:L"][0].P.rbegin(), slice.lines["V+L:L"][0].P.rend());
                region1_VL.P.insert(region1_VL.P.end(), slice.lines["V+L+H: V->L"][0].P.rbegin(), slice.lines["V+L+H: V->L"][0].P.rend());

                region1_VL.X.insert(region1_VL.X.end(), slice.lines["Tmax"][0].X.begin()+1, slice.lines["Tmax"][0].X.end());
                region1_VL.X.insert(region1_VL.X.end(), slice.lines["V+L:L"][0].X.rbegin(), slice.lines["V+L:L"][0].X.rend());
                region1_VL.X.insert(region1_VL.X.end(), slice.lines["V+L+H: V->L"][0].X.rbegin(), slice.lines["V+L+H: V->L"][0].X.rend());

                region1_VL.H.insert(region1_VL.H.end(), slice.lines["Tmax"][0].H.begin()+1, slice.lines["Tmax"][0].H.end());
                region1_VL.H.insert(region1_VL.H.end(), slice.lines["V+L:L"][0].H.rbegin(), slice.lines["V+L:L"][0].H.rend());
                region1_VL.H.insert(region1_VL.H.end(), slice.lines["V+L+H: V->L"][0].H.rbegin(), slice.lines["V+L+H: V->L"][0].H.rend());
                slice.regions["V+L"].push_back(region1_VL);
                // region of VL: part2
                Polygon_slice region2_VL(slice.lines["V+L:V"][1].T,slice.lines["V+L:V"][1].P,slice.lines["V+L:V"][1].X,slice.lines["V+L:V"][1].H,"V+L",COLOR_VL,COLOR_gainsboro);
                //cout<<"XVmax: "<<slice.lines["V+L:V"][1].X[0]<<endl;
                region2_VL.T.insert(region2_VL.T.end(), slice.lines["V+L+H: V->L"][1].T.begin(), slice.lines["V+L+H: V->L"][1].T.end());
                region2_VL.T.insert(region2_VL.T.end(), slice.lines["V+L:L"][1].T.rbegin(), slice.lines["V+L:L"][1].T.rend());

                region2_VL.P.insert(region2_VL.P.end(), slice.lines["V+L+H: V->L"][1].P.begin(), slice.lines["V+L+H: V->L"][1].P.end());
                region2_VL.P.insert(region2_VL.P.end(), slice.lines["V+L:L"][1].P.rbegin(), slice.lines["V+L:L"][1].P.rend());

                region2_VL.X.insert(region2_VL.X.end(), slice.lines["V+L+H: V->L"][1].X.begin(), slice.lines["V+L+H: V->L"][1].X.end());
                region2_VL.X.insert(region2_VL.X.end(), slice.lines["V+L:L"][1].X.rbegin(), slice.lines["V+L:L"][1].X.rend());

                region2_VL.H.insert(region2_VL.H.end(), slice.lines["V+L+H: V->L"][1].H.begin(), slice.lines["V+L+H: V->L"][1].H.end());
                region2_VL.H.insert(region2_VL.H.end(), slice.lines["V+L:L"][1].H.rbegin(), slice.lines["V+L:L"][1].H.rend());
                slice.regions["V+L"].push_back(region2_VL);

                // region of VH
                Polygon_slice region_VH(slice.lines["V+H"][0].T, slice.lines["V+H"][0].P, slice.lines["V+H"][0].X, slice.lines["V+H"][0].H, "V+H", COLOR_VH, COLOR_lightgray);
                region_VH.T.insert(region_VH.T.end(), slice.lines["V+L+H: H->V"][0].T.rbegin(), slice.lines["V+L+H: H->V"][0].T.rend());
                region_VH.T.insert(region_VH.T.end(), slice.lines["V+L+H: H->V"][1].T.begin(), slice.lines["V+L+H: H->V"][1].T.end());

                region_VH.P.insert(region_VH.P.end(), slice.lines["V+L+H: H->V"][0].P.rbegin(), slice.lines["V+L+H: H->V"][0].P.rend());
                region_VH.P.insert(region_VH.P.end(), slice.lines["V+L+H: H->V"][1].P.begin(), slice.lines["V+L+H: H->V"][1].P.end());

                region_VH.X.insert(region_VH.X.end(), slice.lines["V+L+H: H->V"][0].X.rbegin(), slice.lines["V+L+H: H->V"][0].X.rend());
                region_VH.X.insert(region_VH.X.end(), slice.lines["V+L+H: H->V"][1].X.begin(), slice.lines["V+L+H: H->V"][1].X.end());

                region_VH.H.insert(region_VH.H.end(), slice.lines["V+L+H: H->V"][0].H.rbegin(), slice.lines["V+L+H: H->V"][0].H.rend());
                region_VH.H.insert(region_VH.H.end(), slice.lines["V+L+H: H->V"][1].H.begin(), slice.lines["V+L+H: H->V"][1].H.end());
                slice.regions["V+H"].push_back(region_VH);

                // region of LH: part1
                Polygon_slice region1_LH(slice.lines["L+H"][0].T, slice.lines["L+H"][0].P, slice.lines["L+H"][0].X, slice.lines["L+H"][0].H, "L+H", COLOR_LH, COLOR_lightgray);
                region1_LH.T.insert(region1_LH.T.end(), slice.lines["V+L+H: L->H"][0].T.rbegin(), slice.lines["V+L+H: L->H"][0].T.rend());
                region1_LH.P.insert(region1_LH.P.end(), slice.lines["V+L+H: L->H"][0].P.rbegin(), slice.lines["V+L+H: L->H"][0].P.rend());
                region1_LH.X.insert(region1_LH.X.end(), slice.lines["V+L+H: L->H"][0].X.rbegin(), slice.lines["V+L+H: L->H"][0].X.rend());
                region1_LH.H.insert(region1_LH.H.end(), slice.lines["V+L+H: L->H"][0].H.rbegin(), slice.lines["V+L+H: L->H"][0].H.rend());
                slice.regions["L+H"].push_back(region1_LH);

                // region of LH: part2
                    // low boundary of H
                T = xThermal::cxThermal::linspace(Tmin(), Tmin(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                X = xThermal::cxThermal::linspace(1E-20, slice.lines["L+H"][1].X[0], nPoints, true);
                H = H_phase(T,P,X,Liquid);
                slice.lines["Tmin"].push_back(Line_slice(T,P,X,H,"T=Tmin",COLOR_k, "dashed", 0.5));
                H = xThermal::cxThermal::linspace(slice.lines["L+H"][1].H[0], m_NaCl->H_Solid(Tmin(), P0), nPoints);
                X = xThermal::cxThermal::linspace(slice.lines["L+H"][1].X[0], 1.0, nPoints);
                slice.lines["Tmin"].push_back(Line_slice(T,P,X,H,"T=Tmin",COLOR_k, "dashed", 0.5));
                Polygon_slice region2_LH(slice.lines["L+H"][1].T, slice.lines["L+H"][1].P, slice.lines["L+H"][1].X, slice.lines["L+H"][1].H, "V+H", COLOR_LH, COLOR_lightgray);
                region2_LH.T.insert(region2_LH.T.end(),slice.lines["V+L+H: L->H"][1].T.begin(),slice.lines["V+L+H: L->H"][1].T.end());
                region2_LH.T.insert(region2_LH.T.end(),slice.lines["Tmin"][1].T.rbegin(),slice.lines["Tmin"][1].T.rend());

                region2_LH.P.insert(region2_LH.P.end(),slice.lines["V+L+H: L->H"][1].P.begin(),slice.lines["V+L+H: L->H"][1].P.end());
                region2_LH.P.insert(region2_LH.P.end(),slice.lines["Tmin"][1].P.rbegin(),slice.lines["Tmin"][1].P.rend());

                region2_LH.X.insert(region2_LH.X.end(),slice.lines["V+L+H: L->H"][1].X.begin(),slice.lines["V+L+H: L->H"][1].X.end());
                region2_LH.X.insert(region2_LH.X.end(),slice.lines["Tmin"][1].X.rbegin(),slice.lines["Tmin"][1].X.rend());

                region2_LH.H.insert(region2_LH.H.end(),slice.lines["V+L+H: L->H"][1].H.begin(),slice.lines["V+L+H: L->H"][1].H.end());
                region2_LH.H.insert(region2_LH.H.end(),slice.lines["Tmin"][1].H.rbegin(),slice.lines["Tmin"][1].H.rend());
                slice.regions["L+H"].push_back(region2_LH);

                // region of T>Tmax
                Polygon_slice region_Tmax(slice.lines["Tmax"][0].T, slice.lines["Tmax"][0].P, slice.lines["Tmax"][0].X, slice.lines["Tmax"][0].H, "T>1000$^{\\circ}$C", COLOR_lightgray, COLOR_lightgray);
                region_Tmax.T.insert(region_Tmax.T.end(), slice.lines["Tmax"][1].T.begin(), slice.lines["Tmax"][1].T.end());
                region_Tmax.P.insert(region_Tmax.P.end(), slice.lines["Tmax"][1].P.begin(), slice.lines["Tmax"][1].P.end());
                region_Tmax.X.insert(region_Tmax.X.end(), slice.lines["Tmax"][1].X.begin(), slice.lines["Tmax"][1].X.end());
                region_Tmax.H.insert(region_Tmax.H.end(), slice.lines["Tmax"][1].H.begin(), slice.lines["Tmax"][1].H.end());
                region_Tmax.T.insert(region_Tmax.T.end(), 1, slice.lines["Tmax"][0].T[0]);
                region_Tmax.P.insert(region_Tmax.P.end(), 1, slice.lines["Tmax"][0].P[0]);
                region_Tmax.X.insert(region_Tmax.X.end(), 1, slice.lines["Tmax"][1].X[nPoints-1]);
                region_Tmax.H.insert(region_Tmax.H.end(), 1, slice.lines["Tmax"][0].H[0]);
                slice.regions["T>$T_{max}$"].push_back(region_Tmax);
            }else
            {
                slice.points["Critical point"].push_back(Point_slice(T_crit, P0, X_crit, H_phase(T_crit, P0, X_crit, Vapor), "Critical point", "o", COLOR_mfc_CriticalPoint, COLOR_mec_CriticalPoint));
                // VL, LH
                std::vector<double>T, P, X, H;
                // -- [T2, Tmax]
                T = xThermal::cxThermal::linspace(T_crit, Tmax(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                // vapor branch
                XV_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:V"].push_back(Line_slice(T,P,X,H, "V+L:V", COLOR_lc_VL_V));
                // liquid branch
                XL_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:L"].push_back(Line_slice(T,P,X,H, "V+L:L", COLOR_lc_VL_L));
                // halite liquidus
                T = xThermal::cxThermal::linspace(Tmin(), m_NaCl->Melting_T(P0), nPoints);
                X_HaliteLiquidus(T,P,X);
                H = H_phase(T,P,X,Liquid);
                slice.lines["L+H"].push_back(Line_slice(T,P,X,H,"L+H",COLOR_darkgray));
                // Enthalpy boundary at T=Tmax
                T = xThermal::cxThermal::linspace(Tmax(), Tmax(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                X = xThermal::cxThermal::linspace(slice.lines["V+L:V"][0].X[nPoints-1], slice.lines["V+L:L"][0].X[nPoints-1], nPoints);
                H = xThermal::cxThermal::linspace(slice.lines["V+L:V"][0].H[nPoints-1], slice.lines["V+L:L"][0].H[nPoints-1], nPoints);

                slice.lines["Tmax"].push_back(Line_slice(T,P,X,H,"T=Tmax", COLOR_k, "dashed", 0.5));
                X = xThermal::cxThermal::linspace(slice.lines["V+L:L"][0].X[nPoints-1], 1.0, nPoints);
                H = H_phase(T,P,X,Liquid);
                slice.lines["Tmax"].push_back(Line_slice(T,P,X,H,"T=Tmax", COLOR_k, "dashed", 0.5));
                // region of VL
                Polygon_slice region_VL(slice.lines["V+L:V"][0].T,slice.lines["V+L:V"][0].P,slice.lines["V+L:V"][0].X,slice.lines["V+L:V"][0].H,"V+L",COLOR_VL,COLOR_gainsboro);
                region_VL.T.insert(region_VL.T.end(), slice.lines["Tmax"][0].T.begin(), slice.lines["Tmax"][0].T.end());
                region_VL.T.insert(region_VL.T.end(), slice.lines["V+L:L"][0].T.rbegin(), slice.lines["V+L:L"][0].T.rend());

                region_VL.P.insert(region_VL.P.end(), slice.lines["Tmax"][0].P.begin(), slice.lines["Tmax"][0].P.end());
                region_VL.P.insert(region_VL.P.end(), slice.lines["V+L:L"][0].P.rbegin(), slice.lines["V+L:L"][0].P.rend());

                region_VL.X.insert(region_VL.X.end(), slice.lines["Tmax"][0].X.begin(), slice.lines["Tmax"][0].X.end());
                region_VL.X.insert(region_VL.X.end(), slice.lines["V+L:L"][0].X.rbegin(), slice.lines["V+L:L"][0].X.rend());

                region_VL.H.insert(region_VL.H.end(), slice.lines["Tmax"][0].H.begin(), slice.lines["Tmax"][0].H.end());
                region_VL.H.insert(region_VL.H.end(), slice.lines["V+L:L"][0].H.rbegin(), slice.lines["V+L:L"][0].H.rend());
                slice.regions["V+L"].push_back(region_VL);
                // region of LH
                Polygon_slice region_LH(slice.lines["L+H"][0].T, slice.lines["L+H"][0].P, slice.lines["L+H"][0].X, slice.lines["L+H"][0].H, "L+H", COLOR_LH, COLOR_lightgray);
                // low boundary of H
                T = xThermal::cxThermal::linspace(Tmin(), Tmin(), nPoints);
                P = xThermal::cxThermal::linspace(P0, P0, nPoints);
                X = xThermal::cxThermal::linspace(1E-20, slice.lines["L+H"][0].X[0], nPoints, true);
                H = H_phase(T,P,X,Liquid);
                slice.lines["Tmin"].push_back(Line_slice(T,P,X,H,"T=Tmin",COLOR_k, "dashed", 0.5));
                H = xThermal::cxThermal::linspace(slice.lines["L+H"][0].H[0], m_NaCl->H_Solid(Tmin(), P0), nPoints);
                X = xThermal::cxThermal::linspace(slice.lines["L+H"][0].X[0], 1.0, nPoints);
                slice.lines["Tmin"].push_back(Line_slice(T,P,X,H,"T=Tmin",COLOR_k, "dashed", 0.5));

                region_LH.T.insert(region_LH.T.end(),slice.lines["Tmin"][1].T.rbegin(),slice.lines["Tmin"][1].T.rend());
                region_LH.P.insert(region_LH.P.end(),slice.lines["Tmin"][1].P.rbegin(),slice.lines["Tmin"][1].P.rend());
                region_LH.X.insert(region_LH.X.end(),slice.lines["Tmin"][1].X.rbegin(),slice.lines["Tmin"][1].X.rend());
                region_LH.H.insert(region_LH.H.end(),slice.lines["Tmin"][1].H.rbegin(),slice.lines["Tmin"][1].H.rend());
                slice.regions["L+H"].push_back(region_LH);
                // region of T>Tmax
                Polygon_slice region_Tmax(slice.lines["Tmax"][0].T, slice.lines["Tmax"][0].P, slice.lines["Tmax"][0].X, slice.lines["Tmax"][0].H, "T>1000$^{\\circ}$C", COLOR_lightgray, COLOR_lightgray);
                region_Tmax.T.insert(region_Tmax.T.end(), slice.lines["Tmax"][1].T.begin(), slice.lines["Tmax"][1].T.end());
                region_Tmax.P.insert(region_Tmax.P.end(), slice.lines["Tmax"][1].P.begin(), slice.lines["Tmax"][1].P.end());
                region_Tmax.X.insert(region_Tmax.X.end(), slice.lines["Tmax"][1].X.begin(), slice.lines["Tmax"][1].X.end());
                region_Tmax.H.insert(region_Tmax.H.end(), slice.lines["Tmax"][1].H.begin(), slice.lines["Tmax"][1].H.end());
                region_Tmax.T.insert(region_Tmax.T.end(), 1, slice.lines["Tmax"][0].T[0]);
                region_Tmax.P.insert(region_Tmax.P.end(), 1, slice.lines["Tmax"][0].P[0]);
                region_Tmax.X.insert(region_Tmax.X.end(), 1, slice.lines["Tmax"][1].X[nPoints-1]);
                region_Tmax.H.insert(region_Tmax.H.end(), 1, slice.lines["Tmax"][0].H[0]);
                slice.regions["T>$T_{max}$"].push_back(region_Tmax);
            }

            return slice;
        }

        PhaseRegion_Slice cH2ONaCl::Slice_constT(const double T0, size_t nPoints)
        {
            PhaseRegion_Slice slice;
            if (T0>=T_MAX_VLH && T0<=Tmax())
            {
                double P_crit, X_crit;
                P_X_Critical(T0,P_crit, X_crit);
                slice.points["Critical point"].push_back(Point_slice(T0, P_crit, X_crit, H_phase(T0, P_crit, X_crit, Vapor), "Critical point", "o", COLOR_mfc_CriticalPoint, COLOR_mec_CriticalPoint));
                // VL
                std::vector<double>T, P, X, H;
                P = xThermal::cxThermal::linspace(pmin(), P_crit, nPoints);
                T = xThermal::cxThermal::linspace(T0, T0, nPoints);
                // vapor branch
                XV_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:V"].push_back(Line_slice(T,P,X,H, "V+L:V", COLOR_lc_VL_V));
                // liquid branch
                XL_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:L"].push_back(Line_slice(T,P,X,H, "V+L:L", COLOR_lc_VL_L));
                // region of VL
                Polygon_slice region_VL(slice.lines["V+L:V"][0].T,slice.lines["V+L:V"][0].P,slice.lines["V+L:V"][0].X,slice.lines["V+L:V"][0].H,"V+L",COLOR_VL,COLOR_gainsboro);
                region_VL.T.insert(region_VL.T.end(), slice.lines["V+L:L"][0].T.rbegin(), slice.lines["V+L:L"][0].T.rend());
                region_VL.P.insert(region_VL.P.end(), slice.lines["V+L:L"][0].P.rbegin(), slice.lines["V+L:L"][0].P.rend());
                region_VL.X.insert(region_VL.X.end(), slice.lines["V+L:L"][0].X.rbegin(), slice.lines["V+L:L"][0].X.rend());
                region_VL.H.insert(region_VL.H.end(), slice.lines["V+L:L"][0].H.rbegin(), slice.lines["V+L:L"][0].H.rend());
                slice.regions["V+L"].push_back(region_VL);
            }else if(T0>T_MIN_VLH && T0<T_MAX_VLH)
            {
                double P_crit, X_crit;
                P_X_Critical(T0,P_crit, X_crit);
                slice.points["Critical point"].push_back(Point_slice(T0, P_crit, X_crit, H_phase(T0, P_crit, X_crit, Vapor), "Critical point", "o", COLOR_mfc_CriticalPoint, COLOR_mec_CriticalPoint));
                double Pvlh = P_VLH(T0);
                // VL
                std::vector<double>T, P, X, H;
                P = xThermal::cxThermal::linspace(Pvlh, P_crit, nPoints);
                T = xThermal::cxThermal::linspace(T0, T0, nPoints);
                // vapor branch
                XV_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:V"].push_back(Line_slice(T,P,X,H, "V+L:V", COLOR_lc_VL_V));
                // liquid branch
                XL_VL(T,P,X);
                H = H_phase(T,P,X,Vapor);
                slice.lines["V+L:L"].push_back(Line_slice(T,P,X,H, "V+L:L", COLOR_lc_VL_L));
                //VLH
                double Xv_vlh, Xl_vlh;
                X_VLH(T0, Pvlh, Xl_vlh, Xv_vlh);
                double Hv_vlh = H_phase(T0, Pvlh, Xv_vlh, Vapor);
                double Hl_vlh = H_phase(T0, Pvlh, Xl_vlh, Liquid);
                double Hh_vlh = m_NaCl->H_Solid(T0, Pvlh);
                slice.lines["V+L+H: V->L"].push_back(Line_slice({T0,T0},{Pvlh, Pvlh},{Xv_vlh, Xl_vlh},{Hv_vlh, Hl_vlh},"V+L+H: V->L", COLOR_lc_VLH));
                slice.lines["V+L+H: L->H"].push_back(Line_slice({T0, T0}, {Pvlh, Pvlh}, {Xl_vlh, 1.0},{Hl_vlh, Hh_vlh},"V+L+H: L->H", COLOR_lc_VLH));
                // VH
                P = xThermal::cxThermal::linspace(pmin(), Pvlh, nPoints);
                X_VH(T, P, X);
                H = H_phase(T,P,X, Vapor);
                slice.lines["V+H"].push_back(Line_slice(T,P,X,H, "V+H", COLOR_lc_VH));
                // LH
                P = xThermal::cxThermal::linspace(pmin(), pmax(), nPoints);
                X_HaliteLiquidus(T,P,X);
                H = H_phase(T,P,X, Liquid);
                slice.lines["L+H"].push_back(Line_slice(T,P,X,H, "L+H", COLOR_lc_LH));
                //------------------ regions --------------------------
                // region of VL
                Polygon_slice region_VL(slice.lines["V+L:V"][0].T,slice.lines["V+L:V"][0].P,slice.lines["V+L:V"][0].X,slice.lines["V+L:V"][0].H,"V+L",COLOR_VL,COLOR_gainsboro);
                region_VL.T.insert(region_VL.T.end(), slice.lines["V+L:L"][0].T.rbegin(), slice.lines["V+L:L"][0].T.rend());
                region_VL.P.insert(region_VL.P.end(), slice.lines["V+L:L"][0].P.rbegin(), slice.lines["V+L:L"][0].P.rend());
                region_VL.X.insert(region_VL.X.end(), slice.lines["V+L:L"][0].X.rbegin(), slice.lines["V+L:L"][0].X.rend());
                region_VL.H.insert(region_VL.H.end(), slice.lines["V+L:L"][0].H.rbegin(), slice.lines["V+L:L"][0].H.rend());
                slice.regions["V+L"].push_back(region_VL);
                // region of VH
                Polygon_slice region_VH(slice.lines["V+H"][0].T, slice.lines["V+H"][0].P, slice.lines["V+H"][0].X, slice.lines["V+H"][0].H, "V+H", COLOR_VH, COLOR_lightgray);
                region_VH.T.push_back(T0);
                region_VH.P.push_back(Pvlh);
                region_VH.X.push_back(Xl_vlh);
                region_VH.H.push_back(Hl_vlh);
                region_VH.T.push_back(T0);
                region_VH.P.push_back(Pvlh);
                region_VH.X.push_back(1.0);
                region_VH.H.push_back(Hh_vlh);
                region_VH.T.push_back(T0);
                region_VH.P.push_back(pmin());
                region_VH.X.push_back(1.0);
                region_VH.H.push_back(m_NaCl->H_Solid(T0,pmin()));
                slice.regions["V+H"].push_back(region_VH);
                // region of LH
                Polygon_slice region_LH(slice.lines["L+H"][0].T, slice.lines["L+H"][0].P, slice.lines["L+H"][0].X, slice.lines["L+H"][0].H, "L+H", COLOR_LH, COLOR_lightgray);
                region_LH.T.push_back(T0);
                region_LH.P.push_back(slice.lines["L+H"][0].P[nPoints-1]);
                region_LH.X.push_back(1.0);
                region_LH.H.push_back(m_NaCl->H_Solid(T0,slice.lines["L+H"][0].P[nPoints-1]));
                region_LH.T.push_back(T0);
                region_LH.P.push_back(Pvlh);
                region_LH.X.push_back(1.0);
                region_LH.H.push_back(Hh_vlh);

                slice.regions["L+H"].push_back(region_LH);
            } else if(T0>=Tmin() && T0<=T_MIN_VLH)
            {
                WARNING("For case when T in [Tmin, T_MIN_VLH], isothermal section of phase diagram is very simple, there are only two regions separated by halite liquidus.");
                std::vector<double>T,P,X,H;
                P = xThermal::cxThermal::linspace(pmin(), pmax(), nPoints);
                T = xThermal::cxThermal::linspace(T0, T0, nPoints);
                X_HaliteLiquidus(T,P,X);
                // LH
                H = H_phase(T,P,X, Liquid);
                slice.lines["L+H"].push_back(Line_slice(T,P,X,H, "L+H", COLOR_lc_LH));
                // region of LH
                Polygon_slice region_LH(slice.lines["L+H"][0].T, slice.lines["L+H"][0].P, slice.lines["L+H"][0].X, slice.lines["L+H"][0].H, "L+H", COLOR_LH, COLOR_lightgray);
                region_LH.T.push_back(T0);
                region_LH.P.push_back(slice.lines["L+H"][0].P[nPoints-1]);
                region_LH.X.push_back(1.0);
                region_LH.H.push_back(H[nPoints-1]);
                region_LH.T.push_back(T0);
                region_LH.P.push_back(pmin());
                region_LH.X.push_back(1.0);
                region_LH.H.push_back(m_NaCl->H_Solid(T0,pmin()));

                slice.regions["L+H"].push_back(region_LH);
            } else
            {
                char errorinfo[500];
                sprintf(errorinfo,"The input T0 = %f deg.C out of valid range of the xThermal package (%f, %f) K, (%f, %f) deg.C",T0-273.15, Tmin(), Tmax(), Tmin()-273.15, Tmax()-273.15);
                WARNING(string(errorinfo));
                // throw OutOfRangeError(string(errorinfo));
            }
            return slice;
        }

        void cH2ONaCl::q1q2_Tstar_H(std::vector<double> P, std::vector<double> X, std::vector<double> &q1,
                                    std::vector<double> &q2) {
            q1.clear(); q2.clear();
            q1.resize(P.size()); q2.resize(P.size());
            for (int i = 0; i < P.size(); ++i) {
                q1q2_Tstar_H(P[i]/1E5, Wt2Mol(X[i]), q1[i], q2[i]);
            }
        }

        void cH2ONaCl::n1n2_Tstar_V(std::vector<double> P, std::vector<double> X, std::vector<double> &n1,
                                    std::vector<double> &n2) {
            n1.clear(); n2.clear();
            n1.resize(P.size()); n2.resize(P.size());
            for (int i = 0; i < P.size(); ++i) {
                n1n2_Tstar_V(P[i]/1E5, Wt2Mol(X[i]), n1[i], n2[i]);
            }
        }

        std::vector<double> cH2ONaCl::D_Tstar_V(std::vector<double> T, std::vector<double> P, std::vector<double> X)
        {
            std::vector<double> D(T.size());
            for (int i = 0; i < T.size(); ++i) {
                D[i] = D_Tstar_V(T[i], P[i], X[i]);
            }
            return D;
        }

        DeformLinearMesh cH2ONaCl::PhaseBoundary_VL_DeformLinear(PhaseType VaporOrLiquid,int nT, int np)
        {
            double Tmin_VLH_ = Tmin_VLH(); // The min T of VLH surface when p=pmin()
            // 3. construct T-p mesh of the halite liquidus surface
            DeformLinearMesh pb_mesh;
            pb_mesh.T.resize(np);
            pb_mesh.p.resize(np);
            pb_mesh.X.resize(np);
            // T
            for (int ip = 0; ip < np; ++ip) {
                pb_mesh.X[ip].resize(nT);
                pb_mesh.p[ip].resize(nT);
                pb_mesh.T[ip] = linspace(Tmin_VLH_,Tmax(), nT);
            }
            // p
            for (int iT = 0; iT < nT; ++iT) {
                double pmin_, pmax_;
                double T0 = pb_mesh.T[0][iT];
                if (T0<m_constants_Water.T_critical)
                {
                    pmax_ = m_Water->Boiling_p(T0)-1E-5;
                    // std::cout<<"T0: "<<T0-273.15<<"ph2o: "<<pmax_<<", pvlh: "<<P_VLH(T0)<<", dp: "<<pmax_ - P_VLH(T0)<<std::endl;
                } else
                {
                    P_Critical(T0,pmax_);
                }
                if (T0>=T_MAX_VLH)
                {
                    pmin_ = pmin();
                }else
                {
                    pmin_ = P_VLH(T0);
                }
                //refine low-p segment and high-p segment
                int n_refine = int(np*0.4);
                int n_low = n_refine, n_high = n_refine, n_center = np-n_refine*2;
                double ratio_split_high=0.95, ratio_split_low = 0.05;
                double p0_low = pmin_+(pmax_-pmin_)*ratio_split_low;
                double p0_high = pmin_+(pmax_-pmin_)*ratio_split_high;
                std::vector<double>p_low = linspace(pmin_, p0_low, n_low);
                std::vector<double>p_center = linspace(p0_low, p0_high, n_center);
                std::vector<double>p_high = linspace(p0_high, pmax_, n_high);
                // std::vector<double>p_tmp = linspace(pmin_, pmax_, np);
                for (int ip = 0; ip < n_low; ++ip) { pb_mesh.p[ip][iT] = p_low[ip];}
                for (int ip = 0; ip < n_center; ++ip) { pb_mesh.p[ip+n_low][iT] = p_center[ip];}
                for (int ip = 0; ip < n_high; ++ip) { pb_mesh.p[ip+n_low+n_center][iT] = p_high[ip];}
            }
            // X
            switch (VaporOrLiquid) {
                case Vapor:
                    {
                        for (int iT = 0; iT < nT; ++iT) {
                            for (int ip = 0; ip < np; ++ip) {
                                pb_mesh.X[ip][iT] = XV_VL(pb_mesh.T[ip][iT], pb_mesh.p[ip][iT]);
                            }
                        }
                    }
                    break;
                case Liquid:
                    {
                        for (int iT = 0; iT < nT; ++iT) {
                            for (int ip = 0; ip < np; ++ip) {
                                pb_mesh.X[ip][iT] = XL_VL(pb_mesh.T[ip][iT], pb_mesh.p[ip][iT]);
                            }
                        }
                    }
                    break;
            }
            return pb_mesh;
        }

        DeformLinearMesh cH2ONaCl::PhaseBoundary_VLH_DeformLinear(int nT, int nX) {
            int nT_high = int(nT/3);
            int nT_low=nT - nT_high;
            double T_min = T_MIN_VLH, T_max = T_MAX_VLH;
            double T_high=T_min + (T_max - T_min)*0.95;
            int n_log = int(nX*0.4);
            int n_linear = nX - n_log;
            // calculate p_vlh, Xv_vlh, Xl_vlh (may not used)
            std::vector<double> T = linspace(0.0,1.0,nT);
            std::vector<double> T_low_ = linspace(T_min, T_high, nT_low);
            std::vector<double> T_high_ = linspace(T_high, T_max, nT_high);
            for (int i = 0; i < nT_low; ++i) T[i] = T_low_[i];
            for (int i = 0; i < nT_high; ++i) T[i+nT_low] = T_high_[i];
            std::vector<double> P_vlh;
            P_VLH(T,P_vlh);
            std::vector<double> Xl_vlh, Xv_vlh;
            X_VLH(T, P_vlh, Xl_vlh, Xv_vlh);
            //assemble to mesh
            DeformLinearMesh pb_mesh;
            pb_mesh.T.resize(nX);
            pb_mesh.p.resize(nX);
            pb_mesh.X.resize(nX);
            // T
            for (int iX = 0; iX < nX; ++iX) {
                pb_mesh.X[iX].resize(nT);
                pb_mesh.p[iX] = P_vlh;
                pb_mesh.T[iX] = T;
            }
            // X
            double Xcenter = 0.1; //kg/kg
            for (int iT = 0; iT < nT; ++iT) {
                std::vector<double> X_V2L_log10 = linspace(log10(Xv_vlh[iT]),log10(Xcenter),n_log);
                std::vector<double> X_V2L_linear = linspace(Xcenter, 1.0,n_linear);
                for (int iX = 0; iX < n_log; ++iX) pb_mesh.X[iX][iT] = pow(10.0, X_V2L_log10[iX]);
                for (int iX = 0; iX < n_linear; ++iX) pb_mesh.X[iX + n_log][iT] = X_V2L_linear[iX];
            }
            return pb_mesh;
        }

        DeformLinearMesh cH2ONaCl::PhaseBoundary_VH_DeformLinear(int nT, int nP) {
            std::vector<double> T = linspace(T_MIN_VLH, T_MAX_VLH, nT);
            std::vector<double> P_vlh;
            P_VLH(T, P_vlh);
            int np_low = int(nP/3);
            int np_high=nP - np_low;
            DeformLinearMesh pb_mesh;
            pb_mesh.T.resize(nT);
            pb_mesh.p.resize(nT);
            pb_mesh.X.resize(nT);
            // refine and calculate
            std::vector<double> p_low_, p_high_;
            for (int iT = 0; iT < nT; ++iT) {
                pb_mesh.T[iT].resize(nP);
                pb_mesh.p[iT].resize(nP);
                pb_mesh.X[iT].resize(nP);
                double p_low=pmin() + (P_vlh[iT] - pmin())*0.1;
                p_low_ = linspace(pmin(), p_low, np_low);
                p_high_ = linspace(p_low, P_vlh[iT],np_high);
                for (int iP = 0; iP < np_low; ++iP) pb_mesh.p[iT][iP] = p_low_[iP];
                for (int iP = 0; iP < np_high; ++iP) pb_mesh.p[iT][iP + np_low] = p_high_[iP];
                for (int iP = 0; iP < nP; ++iP) {
                    pb_mesh.T[iT][iP] = T[iT];
                    pb_mesh.X[iT][iP] = X_VH(pb_mesh.T[iT][iP], pb_mesh.p[iT][iP]);
                }
            }
            return pb_mesh;
        }

        /**
        * Calculate all phase boundaries in the valid p-T-X range.
        * @param scale_X ["linear" | "log" | "loglinear"]
        * @param ratio_log_to_linear [length of log axis]/[length of linear axis]
        * @param Xcenter [kg/kg], X<Xcenter use log scale; X>=Xcenter use linear scale. It is only used with scale_X="loglinear"
        * @return
        */
        PhaseBoundaries cH2ONaCl::calc_PhaseBoundaries(std::string scale_X, double ratio_log_to_linear, double Xcenter) {
            PhaseBoundaries phaseBoundaries;
            phaseBoundaries.Tmin = Tmin();
            phaseBoundaries.Tmax = Tmax();
            phaseBoundaries.pmin = pmin();
            phaseBoundaries.pmax = pmax();
            phaseBoundaries.Xmin = 1E-16;
            phaseBoundaries.Xmax = 1;
            phaseBoundaries.Xcenter = 0;
            phaseBoundaries.ratio_log_to_linear = 0;
            if(scale_X == "linear")
            {
                phaseBoundaries.scale_X = SCALE_X_linear;
            } else if(scale_X == "log")
            {
                phaseBoundaries.scale_X = SCALE_X_log;
            }else if(scale_X == "loglinear")
            {
                phaseBoundaries.scale_X = SCALE_X_loglinear;
                phaseBoundaries.Xcenter = Xcenter;
                phaseBoundaries.ratio_log_to_linear = ratio_log_to_linear;
            }else
            {
                WARNING("The input scale_X is not recognized: " + scale_X +", the supported option is linear, log, loglinear. Set it to linear");
                phaseBoundaries.scale_X = SCALE_X_linear;
            }
            // Surfaces
            // 1. VL: L
            phaseBoundaries.surfaces.push_back({"VL: liquid branch", "VL_L", PhaseBoundary_VL_DeformLinear(Liquid), "lightblue"});
            // 2. VL: V
            phaseBoundaries.surfaces.push_back({"VL: vapor branch", "VL_V", PhaseBoundary_VL_DeformLinear(Vapor), "orange"});
            // 3. Halite liquidus
            phaseBoundaries.surfaces.push_back({"Halite liquidus", "HL", PhaseBoundary_HaliteLiquidus_DeformLinear(), "green"});
            // 4. VLH
            DeformLinearMesh mesh_vlh = PhaseBoundary_VLH_DeformLinear();
            phaseBoundaries.surfaces.push_back({"VLH", "VLH", mesh_vlh, "gray"});
            // 5. VH
            phaseBoundaries.surfaces.push_back({"VH", "VH", PhaseBoundary_VH_DeformLinear(), "purple"});

            // Lines
            //1. Critical curve
            Line line_critical;
            line_critical.color="red"; line_critical.name = "Critical curve"; line_critical.shortName = "Critical";
            line_critical.T = linspace(m_constants_Water.T_critical+0.1, phaseBoundaries.Tmax, 100);
            P_X_Critical(line_critical.T, line_critical.p, line_critical.X);
            phaseBoundaries.lines.push_back(line_critical);
            //2. VLH: vapor, liquid, halite
            Line line_VLH_V,line_VLH_L,line_VLH_H;
            line_VLH_V.color="green"; line_VLH_V.name = "VLH: vapor"; line_VLH_V.shortName = "VLH_V";
            line_VLH_L.color="blue"; line_VLH_L.name = "VLH: liquid"; line_VLH_L.shortName = "VLH_L";
            line_VLH_H.color="black"; line_VLH_H.name = "VLH: halite"; line_VLH_H.shortName = "VLH_H";
            line_VLH_V.T = mesh_vlh.T[0];  line_VLH_V.p = mesh_vlh.p[0];
            line_VLH_L.T = mesh_vlh.T[0];  line_VLH_L.p = mesh_vlh.p[0];
            line_VLH_H.T = mesh_vlh.T[0];  line_VLH_H.p = mesh_vlh.p[0];
            X_VLH(line_VLH_V.T, line_VLH_V.p, line_VLH_L.X, line_VLH_V.X);
            line_VLH_H.X = linspace(1.0,1.0,line_VLH_H.T.size());
            phaseBoundaries.lines.push_back(line_VLH_V);
            phaseBoundaries.lines.push_back(line_VLH_L);
            phaseBoundaries.lines.push_back(line_VLH_H);
            //add boiling curve of H2O only if the scale_X is linear, because X on the boiling curve of water is 0, it can not display in the log scale diagram
            if (phaseBoundaries.scale_X == SCALE_X_linear)
            {
                Line line_boiling_H2O;
                line_boiling_H2O.color="pink"; line_boiling_H2O.name = "Boiling curve of H2O"; line_boiling_H2O.shortName = "BoilingCurve_H2O";
                line_boiling_H2O.T = linspace(phaseBoundaries.Tmin, m_constants_Water.T_critical, 100);
                line_boiling_H2O.p.resize(line_boiling_H2O.T.size());
                line_boiling_H2O.X.resize(line_boiling_H2O.T.size());
                for (int i = 0; i < line_boiling_H2O.T.size(); ++i) {
                    line_boiling_H2O.p[i] = m_Water->Boiling_p(line_boiling_H2O.T[i]);
                    line_boiling_H2O.X[i] = 0;
                }
                phaseBoundaries.lines.push_back(line_boiling_H2O);
                //Points
                //1. critical point of water
                Point point_crit;
                point_crit.color = "red"; point_crit.name = "Critical point of H2O"; point_crit.shortName = "CriticalPoint_H2O";
                point_crit.T = m_constants_Water.T_critical;
                point_crit.p = m_constants_Water.p_critical;
                point_crit.X = 0;
                phaseBoundaries.points.push_back(point_crit);
            }

            return phaseBoundaries;
        }

        /**
        * @brief Calculate liquid saturation in V+L region.
        *
        * \f{equation}
        * S_l = \frac{\rho_v (X_v - X)}{\rho_v (X_v - X) + \rho_l (X - X_l)}
        * \f}
        *
        * @param X_L [kg/kg]
        * @param X_V [kg/kg]
        * @param Rho_L [kg/m^3]
        * @param Rho_V [kg/m^3]
        * @param X [kg/kg]
        * @return [-]
        */
        double cH2ONaCl::Sl_VL(const double &X_L, const double &X_V, const double &Rho_L, const double &Rho_V,
                               const double &X) {
            return (Rho_V * (X_V - X))/(Rho_V * (X_V - X) + Rho_L * (X - X_L));
        }
        /**
        * @brief Calculate Saturation of phase1 of the mixture of phase1 & phase2.
        *
        * \f{equation}
        * S_{phase1} = \frac{\rho_{phase2} (X_{phase2} - X)}{\rho_{phase2} (X_{phase2} - X) + \rho_{phase1} (X - X_{phase1})}
        * \f}
        *
        * * Calculate \f$ S_l \f$ of two phase V+L: Saturation_Phase1(X_l, X_v, Rho_l, Rho_v, X);
        *
        * \f{equation}
        * S_{l} = \frac{\rho_{v} (X_{v} - X)}{\rho_{v} (X_{v} - X) + \rho_{l} (X - X_{l})}
        * \f}
        *
        * * Calculate \f$ S_h \f$ of two phase V+H: Saturation_Phase1(1.0, X_v, Rho_h, Rho_v, X);
        *
        * * Calculate \f$ S_h \f$ of two phase H+L: Saturation_Phase1(1.0, X_l, Rho_h, Rho_l, X);
        *
        * \f{equation}
        * S_{h} = \frac{\rho_{l} (X_{l} - X)}{\rho_{l} (X_{l} - X) + \rho_{h} (X - X_{h})} = \frac{\rho_{l} (X_{l} - X)}{\rho_{l} (X_{l} - X) + \rho_{h} (X - 1)}
        * \f}
        *
        * @param X_Phase1
        * @param X_Phase2
        * @param Rho_Phase1
        * @param Rho_Phase2
        * @param X
        * @return
        */
        double cH2ONaCl::Saturation_Phase1(const double& X_Phase1, const double& X_Phase2, const double& Rho_Phase1, const double& Rho_Phase2, const double& X)
        {
            return (Rho_Phase2 * (X_Phase2 - X))/(Rho_Phase2*(X_Phase2 - X) + Rho_Phase1*(X-X_Phase1));
        }
        void cH2ONaCl::prop_VL(const double &T, const double &P, const double &X, ThermodynamicProperties &prop) {
            prop.X = X;
            prop.T = T;
            prop.phase = TwoPhase_VL_Water;
            // calculate Xl, Xv
            X_VL(T, P, prop.X_l, prop.X_v);
            // calculate Rho_l, Rho_v
            Rho_phase(T, P, prop.X_l, prop.Rho_l, prop.dRhodP_l, prop.dRhodT_l, Liquid);
            Rho_phase(T, P, prop.X_v, prop.Rho_v, prop.dRhodP_v, prop.dRhodT_v, Vapor);
            prop.Rho_h = m_NaCl->Rho_Solid(T, P);
            // calculate H_l, H_v
            H_phase(T, P, prop.X_l, prop.H_l, Liquid);
            H_phase(T, P, prop.X_v, prop.H_v, Vapor);
            prop.H_h = m_NaCl->H_Solid(T, P);
            // calculate saturation
            prop.S_l = Sl_VL(prop.X_l, prop.X_v, prop.Rho_l, prop.Rho_v, X);
            prop.S_v = 1-prop.S_l;
            prop.S_h = 0;
            // calculate bulk properties
            prop.H = prop.H_l * prop.S_l + prop.H_v * prop.S_v + prop.H_h * prop.S_h;
            prop.Rho = prop.Rho_l * prop.S_l + prop.Rho_v * prop.S_v + prop.Rho_h * prop.S_h;
        }

        /**
        * @brief Calculate isothermal compressibility of two phase.
        *
        * According to the definition of compressibility \f$ \beta = - \frac{1}{V} \left( \frac{\partial V}{\partial p} \right)_T \f$, the compressibility of two phase can be expressed as,
        * \f{equation}
        * \beta = \frac{1}{\rho} \left( \frac{\partial p}{\partial p} \right)_{T,X}
        * \f}
        * where \f$ \rho = S_1 \rho_1 + S_2\rho_2 \f$ is the bulk density. The analytical expression is a bit complicated, so we can calculated it using numerical difference.
        * \f{equation}
        * \beta (T,P,X) = \frac{\rho(T,P + 0.5\delta p, X) - \rho (T, P - 0.5\delta P, X)}{\delta p}
        * \f}
        *
        * @param T
        * @param P
        * @param X
        * @param bulkRho
        * @param compressibility
        * @param dp
        */
        void cH2ONaCl::compressibility_VL(const double& T, const double& p, const double& X, const double& bulkRho, double& compressibility, double dp)
        {
            const double half_dp = dp/2.0;
            double Rho_l, Rho_v, S_l, tmp_dRhodP, tmp_dRhodT, tmp_Cp;
            double Rho_plus, Rho_minus;

            //plus
            const double p_plus = p + half_dp;
            double X_l = XL_VL(T,p_plus), X_v = XV_VL(T,p_plus);
            Rho_phase(T, p_plus, X_l, Rho_l, tmp_dRhodP, tmp_dRhodT, Liquid);
            Rho_phase(T, p_plus, X_v, Rho_v, tmp_dRhodP, tmp_dRhodT, Vapor);
            S_l = Sl_VL(X_l, X_v, Rho_l, Rho_v, X);
            Rho_plus = S_l * Rho_l + (1 - S_l) * Rho_v;

            // minus
            const double p_minus = p - half_dp;
            X_l = XL_VL(T,p_minus), X_v = XV_VL(T,p_minus);
            Rho_phase(T, p_minus, X_l, Rho_l, tmp_dRhodP, tmp_dRhodT, Liquid);
            Rho_phase(T, p_minus, X_v, Rho_v, tmp_dRhodP, tmp_dRhodT, Vapor);
            S_l = Sl_VL(X_l, X_v, Rho_l, Rho_v, X);
            Rho_minus = S_l * Rho_l + (1 - S_l) * Rho_v;

            compressibility = (Rho_plus - Rho_minus)/(dp * bulkRho);
        }

        /**
        * @brief Calculate isothermal compressibility of L+H phase.
        *
        * @param T0
        * @param p0
        * @param X0
        * @param compressibility
        * @param dp
        */
        void cH2ONaCl::compressibility_LH(const double& T, const double& p, const double& X, const double& bulkRho, double& compressibility, double dp)
        {
            const double half_dp = dp/2.0;
            double Rho_vorl, Rho_h, S_h, tmp_dRhodP, tmp_dRhodT, tmp_Cp;
            double Rho_plus, Rho_minus;

            //plus
            const double p_plus = p + half_dp;
            double X_l = X_HaliteLiquidus(T, p_plus);
            Rho_phase(T, p_plus, X_l, Rho_vorl, tmp_dRhodP, tmp_dRhodT, Liquid);
            Rho_h = m_NaCl->Rho_Solid(T,p_plus);
            S_h = Saturation_Phase1(1.0, X_l, Rho_h, Rho_vorl, X);
            Rho_plus = S_h * Rho_h + (1 - S_h) * Rho_vorl;

            // minus
            const double p_minus = p - half_dp;
            X_l = X_HaliteLiquidus(T, p_minus);
            Rho_phase(T, p_minus, X_l, Rho_vorl, tmp_dRhodP, tmp_dRhodT, Liquid);
            Rho_h = m_NaCl->Rho_Solid(T,p_minus);
            S_h = Saturation_Phase1(1.0, X_l, Rho_h, Rho_vorl, X);
            Rho_minus = S_h * Rho_h + (1 - S_h) * Rho_vorl;

            compressibility = (Rho_plus - Rho_minus)/(dp * bulkRho);
        }

        void cH2ONaCl::compressibility_VH(const double& T, const double& p, const double& X, const double& bulkRho, double& compressibility, double dp)
        {
            const double half_dp = dp/2.0;
            double Rho_vorl, Rho_h, S_h, tmp_dRhodP, tmp_dRhodT, tmp_Cp;
            double Rho_plus, Rho_minus;

            //plus
            const double p_plus = p + half_dp;
            double X_v = X_VH(T, p_plus);
            Rho_phase(T, p_plus, X_v, Rho_vorl, tmp_dRhodP, tmp_dRhodT, Vapor);
            Rho_h = m_NaCl->Rho_Solid(T,p_plus);
            S_h = Saturation_Phase1(1.0, X_v, Rho_h, Rho_vorl, X);
            Rho_plus = S_h * Rho_h + (1 - S_h) * Rho_vorl;

            // minus
            const double p_minus = p - half_dp;
            X_v = X_VH(T, p_minus);
            Rho_phase(T, p_minus, X_v, Rho_vorl, tmp_dRhodP, tmp_dRhodT, Vapor);
            Rho_h = m_NaCl->Rho_Solid(T,p_minus);
            S_h = Saturation_Phase1(1.0, X_v, Rho_h, Rho_vorl, X);
            Rho_minus = S_h * Rho_h + (1 - S_h) * Rho_vorl;

            compressibility = (Rho_plus - Rho_minus)/(dp * bulkRho);
        }

        /**
        * @brief Calculate minimum and maximum specific enthalpy of a specific VLH "triangle" for given bulk salinity X.
        * The "VLH triangle" is boundary (edge) of a three-phase region at constant P slice. The max and min H will be calculated for the latent heat of evaporation, see also \cite grant1979compressibility.
        *
        * @param H_v
        * @param H_l
        * @param H_h
        * @param X_v
        * @param X_l
        */
        void cH2ONaCl::HminHmax_VLH_triangle(const double& H_v, const double& H_l, const double& H_h, const double& X_v, const double& X_l, const double& X, double& Hmin, double& Hmax)
        {
            // Hmin
            double X1 = X_l, X2 = (X < X_l ? X_v : 1.0);
            double H1 = H_l, H2 = (X < X_l ? H_v : H_h);
            double ratio = (H2-H1)/(X2-X1); // linear changes from one endpoint to another
            Hmin = H1 + ratio * (X - X1);

            // Hmax
            X1 = X_v, X2 = 1.0;
            H1 = H_v, H2 = H_h;
            ratio = (H2-H1)/(X2-X1);
            Hmax = H1 + ratio * (X - X1);

            // // for the lower triangle (because there are two VLH-triangles in the H-X space at constant pressure), Hmin and Hmax calculated above are correct, but for the upper triangle, the `Hmin` should be the bigger value and `Hmax` should be the smaller value.
            // // However, from the input parameters, we have no ideas on this triangle is the upper or the lower one, so swap Hmax and Hmin if necessary.
            // if(Hmin>Hmax)
            // {
            //     double tmp_H = Hmin;
            //     Hmin = Hmax;
            //     Hmax = tmp_H;
            // }
        }

        void cH2ONaCl::UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X) {
            props.fluidName = name();
            // props.input_pair = TPX;
            props.T = T;    props.p = p;   props.X = X;
            if (X==0)
            {
                m_Water->UpdateState_TPX(props,T,p);
                return;
            }
            // 1. calculate phase region
            getPhaseRegion_TPX(T, p, X, props.phase, props.X_v, props.X_l);
            // 2. calculate phase properties
            switch (props.phase) {
                case SinglePhase_L:
                {
                    Rho_phase(T,p,props.X_l, props.Rho_l, props.dRhodP_l, props.dRhodT_l, Liquid);
                    props.Rho_v =0; props.Rho_h = 0;
                    H_phase(T,p,props.X_l, props.H_l, props.Cp_l, Liquid);
                    props.H_v = 0; props.H_h = 0; props.Cp_v = 0; props.Cp_h = 0;
                    Mu_phase(T,p,props.X_l, props.Mu_l, Liquid);
                    props.Mu_v = 0;
                    props.S_l = 1.0;
                    props.S_v = 0.0;
                    props.S_h = 0.0;
                    props.dRhodP = props.dRhodP_l;
                    props.dRhodT = props.dRhodT_l;
                    props.IsothermalCompressibility_l = props.dRhodP_l/props.Rho_l;
                    props.IsobaricExpansivity_l = -props.dRhodT_l/props.Rho_l;
                    props.IsothermalCompressibility = props.IsothermalCompressibility_l;
                    props.IsobaricExpansivity = props.IsobaricExpansivity_l;
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v + props.Cp_h*props.S_h * props.Rho_h)/props.Rho;
                    // LHR: there is no bulk viscosity but to keep the consistency of the properties, we set it to the viscosity of the liquid phase
                    props.Mu = props.Mu_l;
                }
                    break;
                case TwoPhase_VL:
                {
                    Rho_phase(T,p,props.X_l, props.Rho_l, props.dRhodP_l, props.dRhodT_l, Liquid);
                    H_phase(T,p,props.X_l, props.H_l, props.Cp_l, Liquid);
                    Mu_phase(T, p, props.X_l, props.Mu_l, Liquid);
                    Rho_phase(T,p,props.X_v, props.Rho_v, props.dRhodP_v, props.dRhodT_v, Vapor);
                    H_phase(T,p,props.X_v, props.H_v, props.Cp_v, Vapor);
                    Mu_phase(T, p, props.X_v, props.Mu_v, Vapor);
                    props.S_l = Sl_VL(props.X_l, props.X_v, props.Rho_l, props.Rho_v,props.X);
                    props.S_v = 1.0 - props.S_l;
                    props.S_h = 0.0;
                    props.Rho_h = 0; props.H_h = 0; props.Cp_h = 0;
                    props.IsothermalCompressibility_l = props.dRhodP_l/props.Rho_l;
                    props.IsobaricExpansivity_l = -props.dRhodT_l/props.Rho_l;
                    props.IsothermalCompressibility_v = props.dRhodP_v/props.Rho_v;
                    props.IsobaricExpansivity_v = -props.dRhodT_v/props.Rho_v;
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v + props.Cp_h*props.S_h * props.Rho_h)/props.Rho;
                    // calculate isothermal compressibility
                    compressibility_VL(T,p, X, props.Rho, props.IsothermalCompressibility);
                }
                    break;
                case TwoPhase_VH:
                {
                    Rho_phase(T,p,props.X_v, props.Rho_v, props.dRhodP_v, props.dRhodT_v, Vapor);
                    H_phase(T,p,props.X_v, props.H_v, props.Cp_v, Vapor);
                    Mu_phase(T,p, props.X_v, props.Mu_v, Vapor);
                    props.Rho_h = m_NaCl->Rho_Solid(T,p);
                    props.H_h = m_NaCl->H_Solid(T,p);
                    props.Cp_h = m_NaCl->Cp_Solid(T,p);
                    props.S_h = Saturation_Phase1(1.0, props.X_v, props.Rho_h, props.Rho_v, props.X);
                    props.S_v = 1.0 - props.S_h;
                    props.S_l = 0;
                    props.Rho_l = 0; props.H_l = 0; props.Cp_l = 0; props.Mu_l = 0;
                    props.IsothermalCompressibility_v = props.dRhodP_v/props.Rho_v;
                    props.IsobaricExpansivity_v = -props.dRhodT_v/props.Rho_v;
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v + props.Cp_h*props.S_h * props.Rho_h)/props.Rho;
                    // calculate isothermal compressibility
                    compressibility_VH(T,p, X, props.Rho, props.IsothermalCompressibility);
                }
                    break;
                case TwoPhase_LH:
                {
                    Rho_phase(T,p,props.X_l, props.Rho_l, props.dRhodP_l, props.dRhodT_l, Liquid);
                    H_phase(T,p,props.X_l, props.H_l, props.Cp_l, Liquid);
                    Mu_phase(T,p,props.X_l, props.Mu_l, Liquid);
                    props.Rho_h = m_NaCl->Rho_Solid(T,p);
                    props.H_h = m_NaCl->H_Solid(T,p);
                    props.S_h = Saturation_Phase1(1.0, props.X_l, props.Rho_h, props.Rho_l, props.X);
                    props.S_l = 1.0 - props.S_h;
                    props.S_v = 0;
                    props.Rho_v = 0; props.H_v = 0; props.Cp_v = 0; props.Mu_v = 0;
                    props.IsothermalCompressibility_l = props.dRhodP_l/props.Rho_l;
                    props.IsobaricExpansivity_l = -props.dRhodT_l/props.Rho_l;
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v + props.Cp_h*props.S_h * props.Rho_h)/props.Rho;
                    // calculate isothermal compressibility
                    // BUG: \TODO need to fix bug of compressibility negative for parameters : double T0=438.76642135933298, p0=774575.32599110622, X0=0.42121320366209464;
                    compressibility_LH(T,p, X, props.Rho, props.IsothermalCompressibility);
                    // Use liquid compressibility at this moment
                    // props.IsothermalCompressibility = props.IsothermalCompressibility_l;
                    // // 近似
                    // double a = - props.Rho_h * (props.X_l - X)/pow(props.Rho_l*(props.X_l - X) + props.Rho_h*(X - 1.0), 2.0);
                    // // Rho_Solid(const double& T, const double& P, double& rho, double& dRhodP, double& dRhodT, double& kappa, double& beta);
                    // double  Rho_h, dRhodP_h, dRhodT_h, compressibility_h, expan_h;
                    // m_NaCl->Rho_Solid(T, p, Rho_h, dRhodP_h, dRhodT_h, compressibility_h, expan_h);
                    // double betaT = (props.IsothermalCompressibility_l * props.Rho_l* (props.S_l + (props.Rho_h - props.Rho_l)* a * (props.X_l - X)) + props.S_h*dRhodP_h)/props.Rho;
                    // // cout<<"L+H 近似: "<<betaT*1E5<<", per. bar"<<endl;
                    // props.IsothermalCompressibility = betaT;
                }
                    break;
                case SinglePhase_V:
                {
                    Rho_phase(T,p,props.X_v, props.Rho_v, props.dRhodP_v, props.dRhodT_v, Vapor);
                    H_phase(T,p,props.X_v, props.H_v, props.Cp_v, Vapor);
                    Mu_phase(T, p, props.X_v, props.Mu_v, Vapor);
                    props.S_h = 0;
                    props.S_v = 1.0;
                    props.S_l = 0;
                    props.Rho_l = 0; props.H_l = 0; props.Cp_l = 0; props.Mu_l = 0;
                    props.Rho_h = 0; props.H_h = 0; props.Cp_h = 0;
                    props.IsothermalCompressibility_v = props.dRhodP_v/props.Rho_v;
                    props.IsobaricExpansivity_v = -props.dRhodT_v/props.Rho_v;
                    props.IsothermalCompressibility = props.IsothermalCompressibility_v;
                    props.IsobaricExpansivity = props.IsobaricExpansivity_v;
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v + props.Cp_h*props.S_h * props.Rho_h)/props.Rho;
                    // there is no bulk viscosity but to keep the consistency of the properties, we set it to the viscosity of the vapor phase
                    props.Mu = props.Mu_v;

                }
                    break;
                default:
                    throw NotImplementedError("Unsupported phase region in function cH2ONaCl::UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X): phase index "
                                              +std::to_string(props.phase) + ", name "+ phase_name(props.phase));
                    break;
            }

            // Dummy bulk viscosity
            // props.Mu = (props.Mu_l * props.S_l  + props.Mu_v*props.S_v );
            // props.dRhodP = (props.dRhodP_l * props.S_l  + props.dRhodP_v*props.S_v );
            // props.dRhodT = (props.dRhodT_l * props.S_l  + props.dRhodT_v*props.S_v );
            // props.IsothermalCompressibility = props.dRhodP/props.Rho;
            // props.IsobaricExpansivity = -props.dRhodT/props.Rho;
        }

        PhaseRegion cH2ONaCl::findPhaseRegion_TPX(const double &T, const double &p, const double &X) {
            PhaseRegion phase;
            double X_l, X_v;
            getPhaseRegion_TPX(T,p,X, phase, X_l, X_v);
            return phase;
        }

        double func_T_HPX(double T, void *params)
        {
            Params_Inversion_PTX* param = (Params_Inversion_PTX*)params;
            cH2ONaCl* sw = param->sw;
            double p   = param->P;
            double X = param->X;
            double H = param->H;
            ThermodynamicProperties props;
            sw->UpdateState_TPX(props, T, p, X);
            if(isnan(props.H))throw ValueError("NAN H value in func_T_HPX(double T, void *params), T="+std::to_string(T)+", p="+std::to_string(p)+", X="+std::to_string(X)+", H = "+std::to_string(H));
            return props.H - H;
        }

        double
        cH2ONaCl::T_HPX(const double &H, const double &p, const double &X, const double &Tmin, const double &Tmax) {
            if(isnan(H))return NAN;
            // calculate maximum H for given P,X
            double Hmin =0, Hmax = 0;
            ThermodynamicProperties props;
            UpdateState_TPX(props, T_MAX, p, X);
            Hmax= props.H;
            if(H>props.H)
            {
                static bool ShowWarning_Hmax = true;
                if(ShowWarning_Hmax)WARNING("Value exceed bound in T_HPX(const double &H, const double &p, const double &X, const double &Tmin, const double &Tmax)\nThe input H greater than Hmax of given p="+std::to_string(p)+"Pa, X="+std::to_string(X)+" kg/kg, Hmax="+std::to_string(Hmax)+"J/kg, input H="+std::to_string(H)+"J/kg. The T will be NAN. This warning information will only display once.");
                ShowWarning_Hmax = false;
                return NAN;
            }
            UpdateState_TPX(props, T_MIN, p, X);
            Hmin = props.H;
            if(H<Hmin)
            {
                static bool ShowWarning_Hmin = true;
                if(ShowWarning_Hmin)WARNING("Value exceed bound in T_HPX(const double &H, const double &p, const double &X, const double &Tmin, const double &Tmax)\nThe input H smaller than Hmin of given p="+std::to_string(p)+"Pa, X="+std::to_string(X)+" kg/kg, Hmin="+std::to_string(Hmin)+"J/kg, input H="+std::to_string(H)+"J/kg. The T will be NAN. This warning information will only display once.");
                ShowWarning_Hmin = false;
                return NAN;
            }
            //
            double T_ = 573; // need a guess T ? but now it works quite well.
            double x_lo = Tmin, x_hi = Tmax;
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;

            gsl_function F;
            Params_Inversion_PTX params;
            params.sw = this;
            params.P = p;
            params.X = X;
            params.H = H;

            F.function = &func_T_HPX;
            F.params = &params;
            // printf ("%5d [%.7f, %.7f] %.7f, p=%.4f bar, H=%.4f J/kg, Hmin = %.4f, Hmax = %.4f, X=%.4f\n",iter, x_lo, x_hi,T_,p/1E5, H,Hmin, Hmax, X);
            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                T_ = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/1E4);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
                // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            if(status != GSL_SUCCESS)
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                ERROR("Fatal error in double cH2ONaCl::T_VLH_P0");
            }
            gsl_root_fsolver_free (s);
            return T_;
        }

        /**
         * Get phase region for given H,P,X.
         * @param H
         * @param p
         * @param X
         * @param phase_region
         */
        void cH2ONaCl::getPhaseRegion_HPX(const double &H, const double &p, const double &X, PhaseRegion &phase_region, double& T, double S_lvh[3])
        {
            //1. first check three phase region: VLH, if in the VLH, return T and Xl, Xv, Xh. Otherwise, calculate inversion of T and then call getPhaseRegion_TPX
            double dRhodP, dRhodT;
            if (p<P_Peak_VLH)
            {
                // calculate two possible T for give p
                double T1T2[2];
                T_VLH_P0(p, T1T2[0], T1T2[1]);
                // calculate saturation for two possible T, if in VLH region, there must be a valid combination and a invalid combination.
                double X_lvh[3]={0,0,1.0};      // in order of L,V,H
                double Rho_lvh[3]={0,0,0};      // in order of L,V,H
                double H_lvh[3] = {0, 0, 0};
                // double S_lvh[3]={0,0,0};        // in order of L,V,H
                for (int i = 0; i < 2; ++i) {
                    X_VLH(T1T2[i], p, X_lvh[0], X_lvh[1]);
                    Rho_phase(T1T2[i], p, X_lvh[0], Rho_lvh[0], dRhodP, dRhodT, Liquid);
                    Rho_phase(T1T2[i], p, X_lvh[1], Rho_lvh[1], dRhodP, dRhodT, Vapor);
                    Rho_lvh[2]= m_NaCl->Rho_Solid(T1T2[i], p); // note that in VLH region, halite is in solid phase.
                    H_phase(T1T2[i], p, X_lvh[0], H_lvh[0], Liquid);
                    H_phase(T1T2[i], p, X_lvh[1], H_lvh[1], Vapor);
                    H_lvh[2]= m_NaCl->H_Solid(T1T2[i], p);
                    if(SolveSaturation_ThreePhase_HX(H, X, X_lvh, Rho_lvh, H_lvh, S_lvh))
                    {
                        T=T1T2[i];
                        phase_region = ThreePhase_VLH;
                        return;
                    }
                }
                // if not in VLH region, calculate inversion of T in a general way
                T = T_HPX(H, p, X);
                double Xv, Xl;
                getPhaseRegion_TPX(T, p, X, phase_region, Xv, Xl);
            } else if(p == P_Peak_VLH)
            {
                // calculate two possible T for give p
                double T1T2 = T_Peak_VLH;
                // calculate saturation for two possible T, there must be a valid combination and a invalid combination.
                double X_lvh[3]={0,0,1.0};      // in order of L,V,H
                double Rho_lvh[3]={0,0,0};      // in order of L,V,H
                double H_lvh[3] = {0, 0, 0};
                // double S_lvh[3]={0,0,0};        // in order of L,V,H
                X_VLH(T1T2, p, X_lvh[0], X_lvh[1]);
                Rho_phase(T1T2, p, X_lvh[0], Rho_lvh[0], dRhodP, dRhodT, Liquid);
                Rho_phase(T1T2, p, X_lvh[1], Rho_lvh[1], dRhodP, dRhodT, Vapor);
                Rho_lvh[2]= m_NaCl->Rho_Solid(T1T2, p); // note that in VLH region, halite is in solid phase.
                H_phase(T1T2, p, X_lvh[0], H_lvh[0], Liquid);
                H_phase(T1T2, p, X_lvh[1], H_lvh[1], Vapor);
                H_lvh[2]= m_NaCl->H_Solid(T1T2, p);
                if(SolveSaturation_ThreePhase_HX(H, X, X_lvh, Rho_lvh, H_lvh, S_lvh))
                {
                    T=T1T2;
                    phase_region = ThreePhase_VLH;
                    return;
                }
                // if not in VLH region, calculate inversion of T in a general way
                T = T_HPX(H, p, X);
                double Xv, Xl;
                getPhaseRegion_TPX(T, p, X, phase_region, Xv, Xl);
            }
            else
            {
                T = T_HPX(H, p, X);
                double Xv, Xl;
                getPhaseRegion_TPX(T, p, X, phase_region, Xv, Xl);
            }
        }

        /**
        * Calculate thermodynamic properties and T for given H,P,X in VLH zone.
        * @param props
        * @param H [J/kg]
        * @param p [Pa]
        * @param X [kg/kg]
        * @param T1_T2 [K] Temperature of V+L+H surface for given P
        * @return
        */
        bool cH2ONaCl::UpdateState_HPX_vlh(ThermodynamicProperties& props, const double& H, const double& p, const double& X, const std::vector<double>& T1_T2)
        {
            props.fluidName = name();
            // double T1_T2[2]={0,0};
            // calculate two possible T for give p
            // T_VLH_P0(p, T1_T2[0], T1_T2[1]);
            // calculate saturation for two possible T, if in VLH region, there must be a valid combination and a invalid combination.
            double X_lvh[3]={0,0,1.0};      // in order of L,V,H
            double Rho_lvh[3]={0,0,0};      // in order of L,V,H
            double H_lvh[3] = {0, 0, 0};
            double Cp_lvh[3] = {0, 0,0};
            double S_lvh[3]={0,0,0};        // in order of L,V,H
            for (int i = 0; i < T1_T2.size(); ++i) {
                X_VLH(T1_T2[i], p, X_lvh[0], X_lvh[1]);
                Rho_phase(T1_T2[i], p, X_lvh[0], Rho_lvh[0], props.dRhodP_l, props.dRhodT_l, Liquid);
                Rho_phase(T1_T2[i], p, X_lvh[1], Rho_lvh[1], props.dRhodP_v, props.dRhodT_v, Vapor);
                Rho_lvh[2]= m_NaCl->Rho_Solid(T1_T2[i], p); // note that in VLH region, halite is in solid phase.
                H_phase(T1_T2[i], p, X_lvh[0], H_lvh[0], Cp_lvh[0], Liquid);
                H_phase(T1_T2[i], p, X_lvh[1], H_lvh[1], Cp_lvh[1], Vapor);
                H_lvh[2]= m_NaCl->H_Solid(T1_T2[i], p);
                Cp_lvh[2] = m_NaCl->Cp_Solid(T1_T2[i], p);
                if(SolveSaturation_ThreePhase_HX(H, X, X_lvh, Rho_lvh, H_lvh, S_lvh))
                {
                    props.T=T1_T2[i];
                    props.phase = ThreePhase_VLH;
                    props.X_l = X_lvh[0];
                    props.X_v = X_lvh[1];
                    props.H_l = H_lvh[0];
                    props.H_v = H_lvh[1];
                    props.H_h = H_lvh[2];
                    props.Rho_l = Rho_lvh[0];
                    props.Rho_v = Rho_lvh[1];
                    props.Rho_h = Rho_lvh[2];
                    props.S_l = S_lvh[0];
                    props.S_v = S_lvh[1];
                    props.S_h = S_lvh[2];
                    props.Cp = 0;    //There is no Cp defined in VLH zone.
                    props.Cp_l = Cp_lvh[0];
                    props.Cp_v = Cp_lvh[1];
                    props.Cp_h = Cp_lvh[2];
                    Mu_phase(props.T, p, props.X_l, props.Mu_l, Liquid);
                    Mu_phase(props.T, p, props.X_v, props.Mu_v, Vapor);
                    //3. calculate bulk properties
                    props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v + props.Rho_h*props.S_h;
                    props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v + props.H_h*props.S_h * props.Rho_h)/props.Rho;
                    props.Mu = (props.Mu_l * props.S_l  + props.Mu_v*props.S_v );
                    return true;
                }
            }

            return false;
        }

        /**
        * Calculate T and thermodynamic properties in V+L zone in the low salinity (<0.01 kg/kg) and low pressure (<P_crit_h2o) zone.
        *
        * There are some issues in low salinity region and close VL boundary for T inversion, therefore we need to do some special process to accurately invert T.
        * Note that below T_crit_h2o, the boundary of VL is an isothermal boundary, i.e., for given P, the T is not changed with X, and the T is the boiling T of H2O for given P
        *
        * @param props
        * @param H [J/kg]
        * @param p [Pa]
        * @param X [kg/kg]
        * @return
        */
        bool cH2ONaCl::UpdateState_HPX_vl_lowXlowP(ThermodynamicProperties& props, const double& H, const double& p, const double& X)
        {
            if (X>0.01 || p>m_Water->p_critical())return false;
            double H_vl_min = H_phase(m_Water->Boiling_T(p), p, X, Liquid);
            if (H<H_vl_min)return false;
            double Tl = T_VL(p, X, Liquid);
            double T_boil = m_Water->Boiling_T(p);
            // double Tv = T_VL(p, X, Vapor);
            double H_l, H_v, Cp_l, Cp_v, Rho_v;
            H_phase(Tl, p, X, H_l, Cp_l, Liquid);
            // H_phase(Tv, p, X, H_v, Cp_v, Vapor);
            ThermodynamicProperties props_boil;
            m_Water->Boiling_T(p, props_boil); // Because vapor salinity is very low, it is difficulty invert a T, use pure water property instead.
            H_v = props_boil.H_v;
            Cp_v = props_boil.Cp_v;
            Rho_v = props_boil.Rho_v;
            // check H,P,X in VL zone or not
            if (H>H_v || H<H_l)return false;

            // then H,P,X in VL zone, solve it
            // props.T= Tl + (H - H_l)/(H_v - H_l) * (Tv - Tl); //
            props.phase = TwoPhase_VL;
            props.X_l = XL_VL(Tl, p);
            props.X_v = XV_VL(Tl, p);
            props.H_l = H_l;
            props.H_v = H_v;
            props.H_h = 0; //because in the VL zone, so we don't need to calculate property of halite! set it to 0
            props.Rho_l = Rho_phase(Tl, p, X, Liquid);
            props.Rho_v = Rho_v; //Rho_phase(Tv, p, X, Vapor);
            props.Rho_h = 0;
            props.S_l = (props.Rho_v * ( props.H_v - H))/(H * (props.Rho_l - props.Rho_v) - (props.Rho_l * props.H_l - props.Rho_v * props.H_v));
            props.S_v = 1.0 - props.S_l;
            props.S_h = 0;
            props.Cp = 0;    //There is no Cp defined in VLH zone.
            props.Cp_l = Cp_l;
            props.Cp_v = Cp_v;
            props.Cp_h = 0;
            Mu_phase(props.T, p, props.X_l, props.Mu_l, Liquid);
            Mu_phase(props.T, p, props.X_v, props.Mu_v, Vapor);
            //3. calculate bulk properties
            props.Rho = props.Rho_l * props.S_l + props.Rho_v*props.S_v;
            props.H = (props.H_l * props.S_l * props.Rho_l + props.H_v*props.S_v * props.Rho_v)/props.Rho;
            props.Cp = (props.Cp_l * props.S_l * props.Rho_l + props.Cp_v*props.S_v * props.Rho_v)/props.Rho; //\todo Bulk Cp calculation method need to confirm!
            props.Mu = (props.Mu_l * props.S_l  + props.Mu_v*props.S_v );
            return true;
        }

        void cH2ONaCl::UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X)
        {
            props.fluidName = name();
            // ------ maybe safe check is not needed ----
            ThermodynamicProperties prop;
            UpdateState_TPX(prop, T_MAX, p, X);
            if(H>prop.H)return;
            //-----------------------------------------

            props.H = H;
            props.X = X;
            props.p = p;
            // double T1_T2[2]={0,0};
            //1. first check three phase region: VLH, if in the VLH, return T and Xl, Xv, Xh. Otherwise, calculate inversion of T and then call getPhaseRegion_TPX
            if (p<P_Peak_VLH)
            {
                std::vector<double> T1_T2(2);
                // calculate two possible T for give p
                T_VLH_P0(p, T1_T2[0], T1_T2[1]);
                // try to solve VLH zone
                if(UpdateState_HPX_vlh(props, H, p, X, T1_T2))return;

                // if not in VLH region, calculate inversion of T in a general way
                // do general inversion
                double T_inv = T_HPX(H, p, X);
                UpdateState_TPX(props, T_inv, p, X);
                // if the bulk H calculated from the inverted T is far away from the input H (e.g., error > 1.0), it should be close to the V+L boundary , wee need to make a special process.
                if (fabs(props.H - H)>1.0)UpdateState_HPX_vl_lowXlowP(props, H, p, X);
                props.T = T_inv;

            } else if(p == P_Peak_VLH)
            {
                // try to solve VLH zone
                std::vector<double> T1_T2 = {T_Peak_VLH};
                if(UpdateState_HPX_vlh(props, H, p, X, T1_T2))return;

                // if not in VLH region, calculate inversion of T in a general way
                // do general inversion
                double T_inv = T_HPX(H, p, X);
                UpdateState_TPX(props, T_inv, p, X);
                // if the bulk H calculated from the inverted T is far away from the input H, it should be close to the V+L boundary , wee need to make a special process.
                if (fabs(props.H - H)>1.0)UpdateState_HPX_vl_lowXlowP(props, H, p, X);
                props.T = T_inv;
            }
            else
            {
                props.T = T_HPX(H, p, X);
                UpdateState_TPX(props, props.T, p, X);
            }
        }

        /**
         * Solve saturation \f$S_l, S_v, S_h \f$ for given \f$ X_l, X_v, X_h, \rho_l, \rho_v, \rho_h \f$ and bulk \f$ H, X \f$. See equation (6-9) of \cite vehling2021brine.
         *
         * * Definition of saturation, bulk density and bulk specific enthalpy.
         *
         * \f{align}
         * \rho = S_l \rho_l + S_v \rho_v + S_h \rho_h \ \ \ \text{(1)}\\
         * H = \frac{S_l \rho_l H_l + S_v \rho_v H_v + S_h \rho_h H_h}{\rho} \ \ \ \text{(2)}\\
         * X = \frac{S_l \rho_L X_l + S_v \rho_v X_v + S_h \rho_h X_h}{\rho} \ \ \ \text{(3)}\\
         * S_l + S_v + S_h = 1 \ \ \ \text{(4)}
         * \f}
         *
         * Substitute equation (1) and (4) into equation (2) and (3),
         * \f{align}
         * (H\rho_l - H_l\rho_l){\color{blue}{S_l}} + (H\rho_v - H_v\rho_v){\color{blue}{S_v}} + (H\rho_h - H_h\rho_h)(1-{\color{blue}{S_l - S_v}}) = 0  \ \ \ \text{(5)} \ \ \ \rightarrow  A{\color{blue}{S_l}} + B{\color{blue}{S_v}} + C(1-{\color{blue}{S_l - S_v}}) = 0\\
         * (X\rho_l - X_l\rho_l){\color{blue}{S_l}} + (X\rho_v - X_v\rho_v){\color{blue}{S_v}} + (X\rho_h - X_h\rho_h)(1-{\color{blue}{S_l - S_v}}) = 0  \ \ \ \text{(6)} \ \ \ \rightarrow  D{\color{blue}{S_l}} + E{\color{blue}{S_v}} + F(1-{\color{blue}{S_l - S_v}}) = 0
         * \f}
         *
         * Rearrange equation (5) and (6),
         *
         * \f{align}
         *  (H\rho_l - H_l\rho_l - H\rho_h + H_h\rho_h){\color{blue}{S_l}} + (H\rho_v - H_v\rho_v - H\rho_h + H_h\rho_h){\color{blue}{S_v}} = H_h\rho_h - H\rho_h  \ \ \ \text{(7)} \ \ \ \rightarrow AA{\color{blue}{S_l}} + BB{\color{blue}{S_v}} = CC \\
         *  (X\rho_l - X_l\rho_l - X\rho_h + X_h\rho_h){\color{blue}{S_l}} + (X\rho_v - X_v\rho_v - X\rho_h + X_h\rho_h){\color{blue}{S_v}} = X_h\rho_h - X\rho_h  \ \ \ \text{(8)} \ \ \ \rightarrow DD{\color{blue}{S_l}} + EE{\color{blue}{S_v}} = FF
         * \f}
         *
         * Therefore,
         * \f{align}
         * {\color{blue}{S_v}} &= \frac{CC-\frac{AA}{DD}FF}{BB - \frac{AA}{DD}EE} \\
         * {\color{blue}{S_l}} &= \frac{FF - EE {\color{blue}{S_v}}}{DD} \\
         * {\color{blue}{S_h}} &= 1-{\color{blue}{S_v}}-{\color{blue}{S_l}}
         * \f}
         *
         * @param H [J/kg]
         * @param X [kg/kg]
         * @param X_lvh [kg/kg] [X_l, X_v, X_h]
         * @param Rho_lvh [kg/m^3] [Rho_l, Rho_v, Rho_h]
         * @param S_lvh [-] [S_l, S_v, S_h]
         * @return Validity
         */
        bool cH2ONaCl::SolveSaturation_ThreePhase_HX(const double &H, const double &X, const double *X_lvh, const double *Rho_lvh, const double *H_lvh, double *S_lvh)
        {
            // equation (5,6)
            double A = (H - H_lvh[0])*Rho_lvh[0];
            double B = (H - H_lvh[1])*Rho_lvh[1];
            double C = (H - H_lvh[2])*Rho_lvh[2];
            double D = (X - X_lvh[0])*Rho_lvh[0];
            double E = (X - X_lvh[1])*Rho_lvh[1];
            double F = (X - X_lvh[2])*Rho_lvh[2];
            // equation (7,8)
            double AA = A - C;
            double BB = B - C;
            double CC = -C;
            double DD = D - F;
            double EE = E - F;
            double FF = -F;
            // result
            double AAbyDD = AA/DD;
            S_lvh[1] = (CC - AAbyDD * FF)/(BB - AAbyDD * EE);
            S_lvh[0] = (FF - EE * S_lvh[1])/DD;
            S_lvh[2] = 1.0 - S_lvh[1] - S_lvh[0];
            // check validity
            if((S_lvh[0]>=0 && S_lvh[0]<=1) && (S_lvh[1]>=0 && S_lvh[1]<=1) && (S_lvh[2]>=0 && S_lvh[2]<=1)) {
                return true;
            }
            return false;
        }

        // ThermodynamicProperties cH2ONaCl::UpdateState_TPX(const double &T, const double &p, const double &X) {
        //     ThermodynamicProperties props;
        //     UpdateState_TPX(props, T, p, X);
        //     return props;
        // }
        //
        // ThermodynamicProperties cH2ONaCl::UpdateState_HPX(const double &H, const double &p, const double &X) {
        //     ThermodynamicProperties props;
        //     UpdateState_HPX(props, H, p, X);
        //     return props;
        // }

    };

};