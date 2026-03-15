/**
 * @file NaCl.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief EOS of H2O in formula of NaCl which is based on \cite prost
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef NaCl_xThermal_H
#define NaCl_xThermal_H

#include "IAPS84.h"
#include "IAPWS95_CoolProp.h"
#include "IAPWS95.h"

namespace xThermal
{
    namespace NaCl
    {
        /**
         * @defgroup CONSTANTS_NACL Constants of pure NaCl.
         * 
         * see last pargraph of section 3.2 in \cite Driesner2007Part1.
         * @{
         */
        /** molar mass  [kg/mol] */
        double const M = 0.058443; 
        /** Temperature at triple point of NaCl [\f$ ^{\circ} C\f$] */
        double const T_Triple_C = 800.7; 
        /** Temperature at triple point of NaCl, [K] */ 
        double const T_Triple = T_Triple_C + 273.15; 
        /** Pressure at triple point of NaCl, [Pa] */
        double const P_Triple = 50; 
        /** log_{10}^{P_{Triple}} */
        double const log10_P_Triple = log10(P_Triple);
        /** \f$ \frac{1}{T_{Triple}} \f$ */
        double const inv_T_Triple = 1.0/T_Triple;
        double const T_Triple_sqr = T_Triple*T_Triple;
        double const T_Triple_sqr_C = T_Triple_C*T_Triple_C;
        /** Clapeyron slope \f$ \frac{dP}{dT} = \frac{1}{a} ~ Pa/K \f$ of the NaCl melting curve (see Tab. 3 of \cite Driesner2007Part1), which is used in \f$ \Delta H_{fus} \f$ calculation (Eq. 28 of \cite Driesner2007Part2) \note Notice the unit! */
        double const slope_Clapeyron = 1.0/2.4726e-7;
        /** Reference point and enthalpy of halite. The JANAF table can be found at https://janaf.nist.gov. While this reference point data comes from Falko Vehling(personal communication: he said this table comes from an attached table to a paper) */
        struct H_halite_ref
        {
            const double H = 9.415867359e4; /**< Enthalpy of halite at reference point (T0, P0). \todo need to confirm the reference point of halite enthalpy. */
            const double T_C = 100;   /**< Temperature [deg.C] at reference point */
            const double P = 100E5;   /**< Pressure [Pa] at reference point */
            double R0;    /**< H1 = R0 + R1*T + R2*T^2 + R3*T^3 */
            double R1; /**< Coefficient R1 as function of P in equation of H1: R1 = R1[0] + R1[1]*P0 + R1[2]*P0^2 */
            double R2; /**< Coefficient R2 as function of P in equation of H1: R2 = R2[0] + R2[1]*P0 */
            double R3; /**< Coefficient R2 as function of P in equation of H1: R2 = R1[0] + R2[1]*P0 */
        };
        /** @} */

        /**
         * @brief Parameters for halite and liquid NaCl volumetric properties. See table 3 of \cite Driesner2007Part2.
         * 
         * \note In order to make the input pressure with unit [Pa], change \f$l_3, l_4 \f$ to be 0.005727E-5 and 0.002715E-5, respectively. Because the pressure unit in equation (1-3) is [bar]. In addition, change \f$ m_4, m_5 \f$ to be -1.5259E-10, 5.5058E-13 for the same reason. (See equation 5,6 of \cite Driesner2007Part2)
         */
        struct Coeff_Rho
        {
            const double l[6] = {2170.4, -0.24599, -9.5797E-05, 0.005727E-5, 0.002715E-5, 733.4};
            const double m[6] = {58443, 23.772, 0.018639, -1.9687E-06, -1.5259E-10, 5.5058E-13};
        };

        /**
         * @brief Parameters for enthalpies. See table 5 of \cite Driesner2007Part2.
         * 
         * \note r[3] is a function of T, need to updated inside function \link Cp \endlink.
         * 
         * \note Again, in order to make input P with unit [Pa], multiply 1E-5 and 1E-10 to \f$ r_3, r_4 \f$, respectively, see equation 30 of \cite Driesner2007Part2.
         * 
         * \warning Why the r4 is written as summation of some tiny values in the Tab. 5 ?
         * 
         * \todo Need to confirm r4 in Tab.5 is correct
         * 
         */
        struct Coeff_H
        {
            double r[5] = {1148.81, 0.275774, 8.8103E-05, 
                            0,
                            5.29063E-18 - 9.63084E-21 + 6.50745E-23};
            const double r3[3] = {-0.0017099E-5, - 3.82734E-11, - 8.65455E-14};
        };
        
        /**
         * @brief Class of NaCl EOS.
         * 
         * \image html Melting_Sublimation_Boiling_NaCl.svg The boiling, sublimation and melting curve of NaCl. width=50%
         * 
         */

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
            double Tmin(){return m_constants.Tmin; };                  /**< Get the minimum temperature in K */
            double Tmax(){return m_constants.Tmax;};                          /**< Get the maximum temperature in K.Same as NaCl-H2O */
            double pmin(){return m_constants.pmin; };                             /**< Get the minimum pressure in Pa. Same as NaCl-H2O */
            double pmax(){return m_constants.pmax;};                  /**< Get the maximum pressure in Pa */
            double Ttriple(){return m_constants.Ttriple;};               /**< Get the triple point temperature in K */
            double T_critical(){throw NotImplementedError("T_critical is not available for NaCl");};            /**< Return the critical temperature in K */
            double p_critical(){throw NotImplementedError("T_critical is not available for NaCl");};            /**< Return the critical pressure in Pa */
            double rhomass_critical(){throw NotImplementedError("T_critical is not available for NaCl");};      /**< Return the critical mass density in \f$ kg/m^3 \f$ */
            double molar_mass(){return m_constants.molar_mass;};       /**< Return the molar mass in kg/mol */

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