/**
 * @file NaCl.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of NaCl.
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

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

        /**
         * @brief Update state by given T,p as input pair.
         * 
         * \warning Here we only consider liquid/solid region, because for geological conditions, the vapor region will never be reached.
         * The minimum valid pressure \link pmin \endlink is limited to 1E5.
         * 
         * @param T [K]
         * @param p [Pa]
         */
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
        
        /**
         * @brief Calculate coefficient \f$ q_1, q_2 \f$ of \f$ T^*_h \f$ (Eq.22) at condition of \f$ X_{NaCl} =1 \f$, see equation 25,26 of \cite Driesner2007Part2.
         * The \f$ q_{1, X_{NaCl}} \f$ will be used to calculate \f$ T^*_h \f$ for enthalpy of liquid (melting) NaCl.
         * 
         * \note The unit of P in the original formula (Eq. 25, 26) is bar. So change the second and third coefficient.
         * 
         * @param P [Pa]
         * @param q1 
         * @param q2 
         */
        void cNaCl::get_q1q2_Tstar_H(const double& P, double& q1, double& q2)
        {
            double P_sqr = P * P;

            q1 = 47.9048 - 9.36994E-8 * P + 6.51059E-16 * P_sqr;
            q2 = 0.241022 + 3.45087E-10 * P - 4.28356E-19 * P_sqr;
        }

        /**
         * @brief Calculate \f$ T^*_h \f$ of liquid NaCl. See equation 22 of \cite Driesner2007Part2.
         * 
         * \note The unit of T in the original formula is \f$ ^{\circ}C \f$.
         * 
         * @param q1 see \link get_q1q2_Tstar_H \endlink
         * @param q2 see \link get_q1q2_Tstar_H \endlink
         * @param T [K]
         * @return double [K]
         */
        double cNaCl::Tstar_H(const double& q1, const double& q2, const double& T)
        {
            return q1 + q2*(T-273.15) + 273.15;
        }

        /**
         * @brief The T-P relations of the boiling curve of NaCl. See equation 3, table 3 and figure 4 of \cite Driesner2007Part1.
         * 
         * \note The unit of result is Pa. Unit of #T_Triple and #P_Triple is [K] and [Pa], respectively.
         * 
         * \image html NaCl/Melting_Sublimation_Boiling_NaCl.svg "The boiling, sublimation and melting curve of NaCl." width=20%.
         * 
         * @param T [K]
         * @return double [Pa]
         */
        double cNaCl::Boiling_p(const double& T)
        {
            return pow(10.0, (log10_P_Triple + 9418.12 * (1.0 / T_Triple - 1.0 /T)));
        }

        /**
         * @brief Boiling temperature as a function of pressure. See also \link Boiling_p \endlink
         * 
         * @param P [Pa]
         * @return double [K] 
         */
        double cNaCl::Boiling_T(const double& P)
        {
            return 1.0/(inv_T_Triple - (log10(P) - log10_P_Triple)/9418.12);
        }

        /**
         * @brief Melting curve. See equation 1, table 3 and figure 3 of \cite Driesner2007Part1.
         * 
         * \note Because the unit of P_Triple and P is [Pa], so the parameter \f$ a=2.4726E-7 \f$.
         * 
         * \image html NaCl/Melting_Sublimation_Boiling_NaCl.svg "The boiling, sublimation and melting curve of NaCl." width=20%.
         * 
         * @param P [Pa]
         * @return double [K]
         */
        double cNaCl::Melting_T(const double& P)
        { 
            return T_Triple + 2.4726e-7*(P - P_Triple);
        }

        /**
         * @brief Melting temperature, the same as \link Melting_T \endlink except temperature unit.
         * 
         * @param P [Pa]
         * @return double [deg.C]
         */
        double cNaCl::Melting_T_C(const double& P)
        {
            return T_Triple_C + 2.4726e-7*(P - P_Triple);
        }

        /**
         * @brief Melting pressure as a function of temperature. See also \link Melting_T \endlink
         * 
         * @param T [K]
         * @return double [Pa] 
         */
        double cNaCl::Melting_p(const double& T)
        {
            return P_Triple + 1.0/2.4726e-7 * (T - T_Triple);
        }

        /**
         * @brief Sublimation curve. See equation 2, table 3 and figure 4 of \cite Driesner2007Part1. 
         * 
         * \image html NaCl/Melting_Sublimation_Boiling_NaCl.svg "The boiling, sublimation and melting curve of NaCl." width=20%.
         * 
         * @param T [K]
         * @return double [Pa] 
         */
        double cNaCl::Sublimation_p(const double& T)
        {
            return pow(10, (log10_P_Triple + 11806.1 * (1.0 / T_Triple - 1.0 / T)));
        };

        /**
         * @brief Sublimation temperature as a function of pressure. See also \link Sublimation_P \endlink
         * 
         * @param P [Pa]
         * @return double [K] 
         */
        double cNaCl::Sublimation_T(const double& P)
        {
            return 1.0/(inv_T_Triple - (log10(P) - log10_P_Triple)/11806.1);
        }

        /**
         * @brief Calculate vapor pressure of NaCl. See also \link Sublimation_P \endlink and \link Boiling_p \endlink.
         * 
         * \image html NaCl/Melting_Sublimation_Boiling_NaCl.svg "The boiling, sublimation and melting curve of NaCl." width=20%.
         * 
         * @param T [K]
         * @return double [Pa] 
         */
        double cNaCl::P_Vapor(const double& T)
        {
            return T<T_Triple ? Sublimation_p(T) : Boiling_p(T);
        }
        
        /**
         * @brief Calculate vapor temperature of NaCl. See also \link Sublimation_T \endlink and \link Boiling_T \endlink.
         * 
         * @param P [Pa]
         * @return double [K] 
         */
        double cNaCl::T_Vapor(const double& P)
        {
            return P<P_Triple ? Sublimation_T(P) : Boiling_T(P); 
        }

        /**
         * @brief Halite density as function of (T, P). See equation (1-3) of \cite Driesner2007Part2.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [kg/m3]
         */
        double cNaCl::Rho_Solid(const double& T, const double& P)
        {
            double T_C = T - 273.15;

            return m_coeff_rho.l[0] + m_coeff_rho.l[1]*T_C +m_coeff_rho.l[2]*T_C*T_C + (m_coeff_rho.l[3] + m_coeff_rho.l[4]*exp(T_C/m_coeff_rho.l[5]))*P;
        }

        /**
        * Calculate halite density \f$ \rho \f$ and derivatives \f$ \left( \frac{\partial \rho}{\partial P} \right)_T\f$, \f$ \left( \frac{\partial \rho}{\partial T} \right)_P \f$
        *
        * * \f$ \rho = l_0 + l_1T + l_2 T^2 + (l_3 + l_4 e^{T/l_5}) P\f$: Halite density as function of (T, P). See equation (1-3) of \cite Driesner2007Part2.
        *
        * * \f$ \left( \frac{\partial \rho}{\partial P} \right)_T = l = l_3 + l_4 e^{T/l_5}\f$
        *
        * * \f$ \left( \frac{\partial \rho}{\partial T} \right)_P = l_1+ 2l_2T+ \frac{l_4}{l_5}e^{T/l_5}P \f$
        *
        * @param T [K]
        * @param P [Pa]
        * @param dRhodP [kg/m3/Pa]
        * @param dRhodT [kg/m3/K]
        * @param kappa [1/Pa] Isothermal Compressibility
        * @param beta [1/K] Isobaric Expansivity
        * @return
        */
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
        
        /**
         * @brief Density of liquid (molten) NaCl as function of (T, P). See equation (4-6) of \cite Driesner2007Part2.
         * 
         * @param T 
         * @param P 
         * @return double 
         */
        double cNaCl::Rho_Liquid(const double& T, const double& P)
        {
            double T_C = T - 273.15;

            return m_coeff_rho.m[0]/((m_coeff_rho.m[1] + m_coeff_rho.m[2]*T_C + m_coeff_rho.m[3]*T_C*T_C) * (1 - 0.1*log(1 + 10*P*(m_coeff_rho.m[4] + m_coeff_rho.m[5]*T_C))));
        }

        /**
        * Calculate melting NaCl density \f$ \rho \f$ and derivatives \f$ \left( \frac{\partial \rho}{\partial P} \right)_T\f$, \f$ \left( \frac{\partial \rho}{\partial T} \right)_P \f$
        *
        * * \f$ \rho = \frac{\rho^0}{1 - 0.1ln(1 + 10P \kappa)}\f$: Density of liquid (molten) NaCl as function of (T, P). See equation (4-6) of \cite Driesner2007Part2.
        *
        * * \f$ v = \frac{1}{\rho} = \frac{1}{\rho^0} \left[ 1 - 0.1ln(1 + 10P \kappa) \right] \f$
        *
        * * \f$ \frac{\partial v}{ \partial P} = - \frac{\kappa}{\rho^0 (1 + 10\kappa P)} \f$
        *
        * * \f$ \left( \frac{\partial \rho}{\partial P} \right)_T = -\frac{1}{v^2} \frac{\partial v}{\partial P} = \frac{\kappa \rho^2}{\rho^0 (1 + 10\kappa P)}\f$
        *
        * * \f$ v = \frac{1}{\rho} = \frac{1}{\rho^0} \left[ 1 - 0.1ln(1 + 10P \kappa) \right] =  \frac{1 - 0.1ln(1 + 10P \kappa)}{m_0} (m_1 + m_2T + m_3 T^2) \f$
        *
        * * \f$ \left( \frac{\partial \rho}{\partial T} \right)_P = -\frac{1}{v^2} \frac{\partial v}{\partial T} =  -\rho^2 \frac{1 - 0.1ln(1 + 10P \kappa)}{m_0} (m_2 + 2m_3T)\f$
        *
        * @param T [K]
        * @param P [Pa]
        * @param dRhodP [kg/m3/Pa]
        * @param dRhodT [kg/m3/K]
        * @param kappa [1/Pa] Isothermal Compressibility
        * @param beta [1/K] Isobaric Expansivity
        * @return
        */
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

        /**
         * @brief Calculate specific enthalpy of liquid NaCl. See equation 21 of \cite Driesner2007Part2.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [J/kg]
         */
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

        /**
         * @brief Specific heat capacity of molten NaCl (liquid phase), which is calculated from Eq. 27 of \cite Driesner2007Part2 when \f$ X_{NaCl} =1 \f$
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [J/kg]
         */
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

        /**
         * @brief Specific heat capacity of halite (solid NaCl). See equation 30 of \cite Driesner2007Part2.
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [J/kg/K]
         */
        double cNaCl::Cp_Solid(const double& T, const double& P)
        {
            double T_C = T - 273.15;
            double T_minus_Ttriple = T - T_Triple;
            m_coeff_H.r[3] = m_coeff_H.r3[0] + m_coeff_H.r3[1]*T_C + m_coeff_H.r3[2]*T_C*T_C;

            return m_coeff_H.r[0] + 2.0*m_coeff_H.r[1]*T_minus_Ttriple + 3.0*m_coeff_H.r[2]*T_minus_Ttriple*T_minus_Ttriple + m_coeff_H.r[3]*P + m_coeff_H.r[4]*P*P;
        }

        /**
         * @brief Calculate enthalpy of fusion \f$ \Delta H_{fus} = \frac{dP}{dT}T\Delta V_{fus} \f$ at pressure = P. See Eq. 28 of \cite Driesner2007Part2.
         * 
         * \f$ \Delta V_{fus} = \frac{1}{\rho_l(T, P)} - \frac{1}{\rho_s(T_{melting, P})} \f$ [\f$ \frac{m^3}{kg} \f$]
         * 
         * \image html NaCl/DeltaH_fus.svg "Enthalpy of fusion along melting curve." width=20%
         * 
         * \note The value of \f$ \Delta H_{fus} \f$ at the triple point is 29.78 which is constent with that in section 3.2 of \cite Driesner2007Part2.
         * 
         * @param P [Pa]
         * @return double [J/kg] 
         */
        double cNaCl::DeltaH_fus(const double& T, const double& P)
        {
            return slope_Clapeyron * T * (1.0/Rho_Liquid(T, P) - 1.0/Rho_Solid(T, P));
        }

        /**
         * @brief Calculate the (T, P) function parts of halite enthalpy formula integrated from halite \f$ C_p \f$ and \f$(\partial H/\partial P)_T \f$ (Eq. 30, 29 of \cite Driesner2007Part2). See also doc of \link H_Solid \endlink.
         * 
         * 
         * @param T [K]
         * @param P [Pa]
         * @return double [J/kg] 
         */
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
        /**
         * @brief Calculate enthalpy component \f$ H_1\f$ of integration along T. See doc of \link H_Solid \endlink
         * 
         * \f{align}
         * H_1(T, P) = \left(r_0 - 2r_1 T_{triple, NaCl} + 3r_2T_{triple, NaCl}^2 + \color{red}{r_{3a}P} + r_4P^2 \right)\color{blue}{T} + \left(r_1 - 3r_2T_{triple, NaCl} + \frac{1}{2}\color{red}{r_{3b}P} \right)\color{blue}{T^2} + \left(r_2 + \frac{1}{3}\color{red}{r_{3c}}P \right)\color{blue}{T^3}
         * \f}
         * 
         * @param T_C [deg.C]
         * @param P [Pa]
         * @return double [J/kg] 
         */
        double cNaCl::_H1(const double& T_C, const double& P)
        {
            // R1*T + R2*T^2 + R3*T^3
            // R1 = R1[0] + R1[1]*P + R1[2]*P^2
             double R1[3]={m_coeff_H.r[0] - 2.0*m_coeff_H.r[1]*T_Triple_C + 3.0*m_coeff_H.r[2]*T_Triple_sqr_C,  m_coeff_H.r3[0], m_coeff_H.r[4]};
            return (m_coeff_H.r[0] - 2.0*m_coeff_H.r[1]*T_Triple_C + 3.0*m_coeff_H.r[2]*T_Triple_sqr_C + m_coeff_H.r3[0]*P + m_coeff_H.r[4]*P*P)*T_C 
                + (m_coeff_H.r[1] - 3.0*m_coeff_H.r[2]*T_Triple_C + 0.5*m_coeff_H.r3[1]*P)*T_C*T_C 
                + (m_coeff_H.r[2] + 1.0/3.0*m_coeff_H.r3[2]*P)*pow(T_C, 3.0);
        }

        /**
         * This function is only called from construction function since these coefficients are constants which are also listed in the paper directly.
         * @param coeff_haliteH_ref
         */
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
        /**
         * This function is only used in paper writing stage.
         */
        void cNaCl::_printCoeffs_haliteH()
        {
            // print R
            printf("-- Coefficients of R for halite enthalpy\n");
            printf("R0: %.8E\n", m_H_h_ref.R0);
            printf("R1: %.8E\n", m_H_h_ref.R1);
            printf("R2: %.8E\n", m_H_h_ref.R2);
            printf("R3: %.8E\n", m_H_h_ref.R3);
        }
        /**
         * @brief Calculate enthalpy component \f$ H_2\f$ of integration along P. See doc of \link H_Solid \endlink
         * 
         * \f{align}
         * H_2^*(T, P) &= \frac{ln\rho_h}{l} \\
         * H_2^{**}(T, P) &=  -\frac{ln\rho_h}{l^2} \frac{l_4}{l_5}e^{T/l_5} + \frac{1}{l\rho_h} \left( l_1 + 2l_2T + P\frac{l_4}{l_5}e^{T/l_5} \right) = \frac{1}{l}\left[ -H_2^*(T, P) \color{red}{\frac{l_4}{l_5}e^{T/l_5}} + \frac{1}{\rho_h}\left( l_1 + 2l_2T + P\color{red}{\frac{l_4}{l_5}e^{T/l_5}} \right) \right]
         * \f}
         * 
         * See equation 1-3 and table 3 of \cite Driesner2007Part2 for \f$ \rho_h \f$ and \f$ l \f$.
         * 
         * \note The unit of \f$ T\f$ in the above equations is \f$ ^{\circ}C \f$. Unlike the pressure \f$ P \f$, it is not possible to simply change the coefficient to make the unit of T to [K], so just keep it as \f$ ^{\circ} C\f$.
         * 
         * @param T_C [deg.C]
         * @param P [Pa]
         * @return double [J/kg] 
         */
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
        /**
         * @brief Calculate specific enthalpy of halite (solid NaCl).
         * 
         * 1. Basic formula: 
         * \f{align}
         * H_{halite}(T, P) = \int_{(T_0,P_0)}^{(T,P)} dH = \int_{(T_0,P_0)}^{(T,P)} \left( \frac{\partial H}{\partial T}\right)_P dT + \left( \frac{\partial H}{\partial P} \right)_T dP = \int_{T_0}^{T} \left( \frac{\partial H}{\partial T}\right)_{P=P_0} dT + \int_{P_0}^{P} \left( \frac{\partial H}{\partial P} \right)_T dP
         * \f}
         * 
         * 2. Integration along T (See Eq. 30 of \cite Driesner2007Part2 for \f$C_p \f$ of halite)
         * \f{align}
         * \int_{T_0}^{T} \left( \frac{\partial H}{\partial T}\right)_{P=P_0} dT &= \int_{T_0}^{T} C_p(T, P_0) dT =
         * \left[
         * \underbrace{\left(r_0 - 2r_1 T_{triple, NaCl} + 3r_2T_{triple, NaCl}^2 + \color{red}{r_{3a}P} + r_4P^2 \right)}_{R_1}\color{blue}{T}
         * + \underbrace{\left(r_1 - 3r_2T_{triple, NaCl} + \frac{1}{2}\color{red}{r_{3b}P} \right)}_{R_2}\color{blue}{T^2}
         * + \underbrace{\left(r_2 + \frac{1}{3}\color{red}{r_{3c}}P \right)}_{R_3}\color{blue}{T^3}
         * \right] \bigg|_{T_0}^{T} = H_1(T, P_0) - H_1(T_0, P_0) \\
         * & = \sum_{i=1}^{3}R_i T^i - R_0, ~~(R_0=\sum_{i=1}^{3}R_i T_0^i = H_1(T_0, P_0))
         * \f}
         * 
         * 3. Integration along P (See Eq. 29 of \cite Driesner2007Part2 for \f$ \left(\frac{\partial H}{\partial P}\right)_T \f$ of halite, and Eq. 1-3 for density \f$ \rho_h \f$ of halite)
         * \f{align}
         * \int_{P_0}^P \left( \frac{\partial H}{\partial P} \right)_T dP = \int_{P_0}^P \left[ V - {\color{orange}{T}}\left(\frac{\partial V}{\partial T} \right)_P \right] dP = \color{green}{\int_{P_0}^P V dP} -\color{orange}{T}\frac{\partial }{\partial T} \left( \color{green}{\int_{P_0}^P V dP} \right) = H_2(T, P) - H_2(T, P_0)
         * \f}
         * 
         * \warning The unit of \f$ \color{orange}{T}\f$ in front of \f$ \left( \frac{\partial V}{\partial T} \right)_P \f$ is [K]
         * 
         * Where
         * \f{align}
         * H_2(T, P) = H_2^*(T, P) - \color{orange}{T} H_2^{**}(T, P)
         * \f}
         * 
         * \f{align}
         * \color{green}{\int_{P_0}^P V dP} = \int_{P_0}^P \frac{1}{lP + \rho^0_h} dP = \frac{1}{l}ln\rho_h\bigg|_{P_0}^P = H_2^*(T, P) - H_2^*(T, P_0)
         * \f}
         * 
         * \f{align}
         * \frac{\partial }{\partial T} \left( \color{green}{\int_{P_0}^P V dP} \right) = \frac{\partial H_2(T, P)}{\partial T} - \frac{\partial H_2(T, P_0)}{\partial T} = H_2^{**}(T, P) - H_2^{**}(T, P_0)
         * \f}
         * 
         * \f{align}
         * H_2^{**}(T, P) = \frac{\partial H_2(T, P)}{\partial T} = \frac{\partial }{\partial T}\left(\frac{ln\rho_h}{l} \right) = -\frac{ln\rho_h}{l^2}\frac{\partial l}{\partial T} + \frac{1}{l\rho_h}\frac{\partial \rho_h}{\partial T} = -\frac{ln\rho_h}{l^2} \frac{l_4}{l_5}e^{T/l_5} + \frac{1}{l\rho_h} \left( l_1 + 2l_2T + P\frac{l_4}{l_5}e^{T/l_5} \right)
         * \f}
         * 
         * Therefore, 
         * \f{align}
         * \boxed{ H_{halite}(T, P) = H_1(T, P) - H_1(T_0, P) + \left[ H_2(T, P) - H_2(T, P_0)\right] + C }, \\
         * \text{where } C=T_{halite}(T_0, P_0) \text{ is the reference value, we select a reference point and value as } \boxed{T_{halite}(T=100^{\circ}, P=100bar)=9.415867359\times10^4 ~ J/kg}
         * \f}
         * 
         * \note The unit of \f$ T\f$ in the above equations is \f$ ^{\circ}C \f$. Unlike the pressure \f$ P \f$, it is not possible to simply change the coefficient to make the unit of T to [K], so just keep it as \f$ ^{\circ} C\f$.
         * 
         * **Some original ideas as following**
         * 
         * \note **How to calculate halite (solid phase) enthalpy at the triple point ?**
         * 1. The enthalpy of molten NaCl (liquid phase) \f$ H_l \f$ can be calculated from Eq. 21 of \cite Driesner2007Part2 by subsituting \f$ X_{NaCl} = 1 \f$.
         * 
         * 2. The enthalpy difference \f$ \Delta H_{fus} = H_l - H_h \f$ of phase changing (solid phase to liquid phase) can be calculated from Eq. 28 of \cite Driesner2007Part2. 
         * 
         *      (1). Clapeyron slope of the NaCl melting curve is \f$ \frac{dP}{dT} = \frac{1}{a} = 4044325.810887325~ Pa/K \f$; (**notice the unit!**) 
         * 
         *      (2). Triple temperature \f$ T_{Triple} = 1073.85~ K \f$; 
         * 
         *      (3). Volume of melting at triple point \f$ \Delta V_{fus} = \frac{1}{\rho_l} - \frac{1}{\rho_h} = \frac{1}{1561.224722497394} - \frac{1}{1912.0183915214473} = 0.00011751525885297855  ~ kg/m^3 = 6.867944273144626\times 10^{-6} ~ m^3/mol\f$ (see \link Rho_Liquid \endlink and \link Rho_Solid \endlink, and #M). 
         * 
         *      Therefore, \f$ \Delta H_{fus} = \frac{dP}{dT}\times T_{Triple}\times \Delta V_{fus} = 4044325.810887325 \times 1073.85 \times 6.867944273144626\times 10^{-6} = 29827.476978550334 kJ/mol \f$ which is constent with value \f$ 29.8 ~ kJ/mol \f$ in \cite Driesner2007Part2. 
         * 
         * 3. Enthalpy of halite at triple point temperature \f$ H_h(T_{Triple}, P) = H_l(T_{Triple}, P) - \Delta H_{fus}(T_{Triple}, P) \f$.
         * 
         * \todo Confirm that \f$ \Delta H_{fus} (T_{melting}, P) = H_l(T_{melting}, P) - H_s(T_{melting}, P) \f$ ?
         * 
         * \note The constant \f$ C(P) \f$ in the above equation
         * 
         * We easily calculate the constant C(P) once get \f$ H_h \f$, so \f$ C(P) = H_h(T_{Triple}, P) - H_{halite}(T_{Triple}, P) \f$.
         * 
         * \warning Just like the description in 2.2.4 of \cite Driesner2007Part2 : "Pure NaCl vapor as a separate phase exists only at pressures substantially lower than those encountered in any geological system", no worries that the (T, P) point occurs in pure vapor region of NaCl (See figure at \link Boiling_p(const double& T) \endlink), i.e. NaCl always present as liquid or solid state in the geological T-P conditions.
         * 
         * @param T 
         * @param P 
         * @return double 
         */
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