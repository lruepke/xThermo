/**
 * @file IAPS84.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief xThermo backend "IAPS84": a C++ wrapper around PROST 4.1, which implements the IAPS-84 equation of state
 *        for water (Haar, Gallagher & Kell 1984; \cite prost).
 *
 * The equation of state, its derivatives, the saturation routines and the viscosity correlation are computed entirely
 * by PROST 4.1 (O. Bauer and O. Engel, TU Hamburg-Harburg, 1998; ThirdParties/PROST), which is distributed under the
 * GNU Library General Public License v2.0 (ThirdParties/PROST/COPYING.LIB-2.0). This wrapper contains only the
 * xThermo interface: state set-up, unit and derivative conversion, phase classification and result mapping. The
 * changes xThermo made to PROST itself are listed in ThirdParties/PROST/XTHERMO_MODIFICATIONS.txt.
 * @version 0.1
 * @date 2022-03-27
 *
 * @copyright Copyright (c) 2022. This wrapper is part of xThermo (MIT licence, see LICENSE). PROST is not covered by
 *            xThermo's licence (see THIRD_PARTY_NOTICES.md).
 *
 */

#ifndef IAPS84_xThermal_H
#define IAPS84_xThermal_H

#include "thermo.h"

namespace xThermal
{
    #define Name_Backend_IAPS84 "IAPS84"

    namespace PROST
    {
        class xTHERMO_VAR cIAPS84 : public cxThermal
        {
        private:
            CONSTENTS_Thermo m_constants;
        public:
            cIAPS84(/* args */);
            cIAPS84(const cIAPS84& water);
            ~cIAPS84();
        private:
            void initialize_data();
            xThermal::PhaseRegion phase_PROST2xThermal(const double& T, const double& p, const void * prop);
        public:
            std::string name(){return Name_Backend_IAPS84;};
            // thermodynamic constants of the model
            double Tmin(){return m_constants.Tmin; };                  /**< Get the minimum temperature in K */
            double Tmax(){return m_constants.Tmax;};                  /**< Get the maximum temperature in K */
            double pmin(){return m_constants.pmin; };                  /**< Get the minimum pressure in Pa */
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
            // void UpdateState_TPX(ThermodynamicPropertiesArray& stateArray, const size_t& N, const double* T, const double* p, const double* X=NULL); //for Matlab API
            ThermodynamicProperties UpdateState_TPX(const double &T, const double &p, const double &X=0); //for Python API
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