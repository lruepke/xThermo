/**
 * @file IAPWS95_CoolProp.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief EOS of H2O in formula of IAPWS95_CoolProp which is based on \cite CoolProp
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifdef USE_COOLPROP

#ifndef IAPWS95_CoolProp_xThermal_H
#define IAPWS95_CoolProp_xThermal_H

#include "thermo.h"

namespace xThermal
{
    // IAPWS95_CoolProp
    #define Name_Backend_IAPWS95_CoolProp "IAPWS95_CoolProp"
    namespace COOLPROP
    {

        class cIAPWS95_CoolProp : public cxThermal
        {
        private:
            CONSTENTS_Thermo m_constants;
        public:
            cIAPWS95_CoolProp(/* args */);
            ~cIAPWS95_CoolProp();

        private:
            void initialize_data();
            PhaseRegion phase_CoolProp2xThermal(const int& phase_);
        public:
            std::string name(){return Name_Backend_IAPWS95_CoolProp;};
            // thermodynamic constants of the model
            double Tmin(){return m_constants.Tmin; };                  /**< Get the minimum temperature in K */
            double Tmax(){return m_constants.Tmax;};                  /**< Get the maximum temperature in K */
            double pmin(){return m_constants.pmin; };                    /**< Get the minimum pressure in Pa */
            double pmax(){return m_constants.pmax;};                  /**< Get the maximum pressure in Pa */
            double Ttriple(){return m_constants.Ttriple;};               /**< Get the triple point temperature in K */
            double T_critical(){return m_constants.T_critical;};            /**< Return the critical temperature in K */
            double p_critical(){return m_constants.p_critical;};            /**< Return the critical pressure in Pa */
            double rhomass_critical(){return m_constants.rhomass_critical;};      /**< Return the critical mass density in \f$ kg/m^3 \f$ */
            double molar_mass(){return m_constants.molar_mass;};       /**< Return the molar mass in kg/mol */

            // update thermodynamic state and properties for given [T,p], or [T,p,X], X is the salinity of H2ONaCl. For pure fluid, X=0 as a default.
            PhaseRegion findPhaseRegion_TPX(const double& T, const double& p, const double& X=0);
            void UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);
            void UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X=0);
            ThermodynamicProperties UpdateState_TPX(const double& T, const double& p, const double& X=0); //for Python API
            double Boiling_p(const double& T);
            double Boiling_T(const double& p);
            double Boiling_p(const double& T, double& rho_l, double& rho_v);
            double Boiling_T(const double& p, double& rho_l, double& rho_v);
            double Boiling_p(const double& T, ThermodynamicProperties& props);
            double Boiling_T(const double& p, ThermodynamicProperties& props);
        };
    
    };

};
#endif

#endif