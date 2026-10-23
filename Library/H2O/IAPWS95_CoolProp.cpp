/**
 * @file IAPWS95_CoolProp.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of IAPWS95_CoolProp.
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifdef USE_COOLPROP

#include "IAPWS95_CoolProp.h"
#include "CoolProp.h"
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"
// ----- LUT related head filess --------
#include "LookUpTableForestI.H"
#include "interpolationI.H"
#include "AMR_LUT_RefineFuncI.H"
//---------------------------------------
namespace xThermal
{
    namespace COOLPROP
    {
        void cIAPWS95_CoolProp::initialize_data()
        {
            shared_ptr<CoolProp::AbstractState> water_CoolProp(CoolProp::AbstractState::factory("HEOS","Water"));
            m_constants.Tmin = water_CoolProp->Tmin();
            m_constants.Tmax = water_CoolProp->Tmax();
            m_constants.pmin = 1E4;
            m_constants.pmax = water_CoolProp->pmax();
            m_constants.Ttriple = water_CoolProp->Ttriple();
            m_constants.T_critical = water_CoolProp->T_critical();
            m_constants.p_critical = water_CoolProp->p_critical();
            m_constants.rhomass_critical = water_CoolProp->rhomass_critical();
            m_constants.molar_mass = water_CoolProp->molar_mass();
        }
        cIAPWS95_CoolProp::cIAPWS95_CoolProp(/* args */)
        {
            initialize_data();
        }

        cIAPWS95_CoolProp::~cIAPWS95_CoolProp()
        {
        }

        xThermal::PhaseRegion cIAPWS95_CoolProp::phase_CoolProp2xThermal(const int& phase_)
        {
            using namespace CoolProp;
            switch (phase_)
            {
            case iphase_liquid: return SinglePhase_L;
            case iphase_supercritical: return Supercritical;
            case iphase_supercritical_gas: return Supercritical_vapor;
            case iphase_supercritical_liquid: return Supercritical_liquid;
            case iphase_critical_point: return Critical;
            case iphase_gas: return SinglePhase_V;
            case iphase_twophase: return TwoPhase_VL_Water;
            case iphase_unknown: return Unknown;
            case iphase_not_imposed: return NotImposed;
            default:
                break;
            }
            return Unknown;
        }

        double cIAPWS95_CoolProp::Boiling_p(const double& T)
        {
            if (T <= T_critical())
            {
                shared_ptr<CoolProp::AbstractState> m_Water_tmp(CoolProp::AbstractState::factory("HEOS","Water"));
                m_Water_tmp->update(CoolProp::QT_INPUTS, 0, T);

                return  m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iP);
            }
            return NAN;
        }

        double cIAPWS95_CoolProp::Boiling_T(const double& p)
        {
            if (p <= p_critical())
            {
                shared_ptr<CoolProp::AbstractState> m_Water_tmp(CoolProp::AbstractState::factory("HEOS","Water"));
                m_Water_tmp->update(CoolProp::PQ_INPUTS, p, 0);

                return  m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iT);
            }
            return NAN;
        }

        double cIAPWS95_CoolProp::Boiling_T(const double& p, double& rho_l, double& rho_v)
        {
            if (p <= p_critical())
            {
                shared_ptr<CoolProp::AbstractState> m_Water_tmp(CoolProp::AbstractState::factory("HEOS","Water"));
                m_Water_tmp->update(CoolProp::PQ_INPUTS, p, 0);
                rho_l = m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iDmass);
                rho_v = m_Water_tmp->saturated_vapor_keyed_output(CoolProp::iDmass);
                return  m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iT);
            }
            return NAN;
        }

        double cIAPWS95_CoolProp::Boiling_p(const double& T, double& rho_l, double& rho_v)
        {
            if (T <= T_critical())
            {
                shared_ptr<CoolProp::AbstractState> m_Water_tmp(CoolProp::AbstractState::factory("HEOS","Water"));
                m_Water_tmp->update(CoolProp::QT_INPUTS, 0, T);
                rho_l = m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iDmass);
                rho_v = m_Water_tmp->saturated_vapor_keyed_output(CoolProp::iDmass);
                return  m_Water_tmp->saturated_liquid_keyed_output(CoolProp::iP);
            }
            return NAN;
        }

        void cIAPWS95_CoolProp::UpdateState_TPX(ThermodynamicProperties &props, const double &T, const double &p, const double& X) {
            props.fluidName = name();
            // Because the update function of CoolProp can not accept PT_INPUTS for two phase, so regard all state close to boiling curve as liquid
            PhaseRegion phase = findPhaseRegion_TPX(T,p, X);
            switch (phase) {
                case TwoPhase_VL_Water: //use liquid properties
                {
                    Boiling_p(T,props);
                    props.T = T;
                    props.p = p;
                    props.Rho = props.Rho_l;
                    props.Rho_l = props.Rho_l;
                    props.Rho_v = props.Rho_v;
                    props.H = props.H_l;
                    props.Cp = 0; //Cp in two phase is set to 0, the same as PROST
                    // props.IsothermalCompressibility = props.IsothermalCompressibility_l;
                    // props.IsobaricExpansivity = props.IsobaricExpansivity_l;
                    props.Mu = props.Mu_l;
                    props.phase = TwoPhase_VL_Water; // or maybe use SinglePhase_L ?
                }
                    break;
                case Supercritical:
                case Supercritical_liquid:
                case Supercritical_vapor:
                case SinglePhase_L:
                case SinglePhase_V:
                case Critical:
                {
                    shared_ptr<CoolProp::AbstractState> water_CoolProp(CoolProp::AbstractState::factory("HEOS","Water"));
                    water_CoolProp->update(CoolProp::PT_INPUTS, p, T);
                    props.T = T;
                    props.p = p;
                    props.Rho = water_CoolProp->rhomass();
                    props.H = water_CoolProp->hmass();
                    props.Cp = water_CoolProp->cpmass();
                    props.Mu = water_CoolProp->viscosity();
                    props.dRhodT = water_CoolProp->first_partial_deriv(CoolProp::iDmass, CoolProp::iT, CoolProp::iP);
                    props.dRhodP = water_CoolProp->first_partial_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iT);
                    props.IsothermalCompressibility = water_CoolProp->isothermal_compressibility();
                    props.IsobaricExpansivity = water_CoolProp->isobaric_expansion_coefficient();
                    props.phase = phase;
                }
                    break;
                default:
                    throw NotImplementedError("The phase is unsupported in UpdateState_TPX function for H2O in UpdateState_TPX function of "+name()
                    +", T="+std::to_string(T)+", p="+std::to_string(p)+", phase: "+ phase_name(phase));
                    break;
            }

        }

        void cIAPWS95_CoolProp::UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X)
        {
            props.fluidName = name();

            shared_ptr<CoolProp::AbstractState> water_CoolProp(CoolProp::AbstractState::factory("HEOS","Water"));
            water_CoolProp->update(CoolProp::HmassP_INPUTS, H, p);
            PhaseRegion phase = phase_CoolProp2xThermal(water_CoolProp->phase());
            props.T = water_CoolProp->T();
            props.p = p;
            props.Rho = water_CoolProp->rhomass();
            props.H = water_CoolProp->hmass();
            props.Cp = water_CoolProp->cpmass(); //Cp in two phase is set to 0, the same as PROST
            props.dRhodT = water_CoolProp->first_partial_deriv(CoolProp::iDmass, CoolProp::iT, CoolProp::iP);
            props.dRhodP = water_CoolProp->first_partial_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iT);
            props.IsothermalCompressibility = water_CoolProp->isothermal_compressibility();
            props.IsobaricExpansivity = water_CoolProp->isobaric_expansion_coefficient();
            props.Mu = water_CoolProp->viscosity();
            props.phase = phase;
            switch (phase) {
                case TwoPhase_VL_Water: //use liquid properties
                {
                    props.Cp = 0; //Cp in two phase is set to 0, the same as PROST
                    props.dRhodT = 0;
                    props.dRhodP = 0;
                    props.IsothermalCompressibility = 0;
                    props.IsobaricExpansivity = 0;
                    // props.Cp_l = water_CoolProp->saturated_liquid_keyed_output(CoolProp::iCpmass);
                    // props.Cp_v = water_CoolProp->saturated_vapor_keyed_output(CoolProp::iCpmass);
                    props.H_l = water_CoolProp->saturated_liquid_keyed_output(CoolProp::iHmass);
                    props.H_v = water_CoolProp->saturated_vapor_keyed_output(CoolProp::iHmass);
                    props.Rho_l = water_CoolProp->saturated_liquid_keyed_output(CoolProp::iDmass);
                    props.Rho_v = water_CoolProp->saturated_vapor_keyed_output(CoolProp::iDmass);
                    props.Mu_l = water_CoolProp->saturated_liquid_keyed_output(CoolProp::iviscosity);
                    props.Mu_v = water_CoolProp->saturated_vapor_keyed_output(CoolProp::iviscosity);
                    props.S_v = water_CoolProp->Q(); // Confirmed: Zhikui 2022.08.01. vapor quality (mol/mol) = kg/kg, note that the mole fraction is the same as mass fraction for single component fluid
                    props.S_l = 1.0 - props.S_v;
                }
                    break;
                case SinglePhase_V:
                    props.S_v = 1.0;
                    props.S_l = 0;
                    break;
                case Supercritical:
                case Supercritical_liquid:
                case Supercritical_vapor:
                case SinglePhase_L:
                case Critical:
                    props.S_v = 1;
                    props.S_l = 1.0;
                    break;
                default:
                    throw NotImplementedError("The phase is unsupported in UpdateState_HPX function of "+name()
                                              +", T="+std::to_string(props.T)+", p="+std::to_string(p)+", phase: "+ phase_name(phase));
                    break;
            }
        }

        double cIAPWS95_CoolProp::Boiling_p(const double &T, ThermodynamicProperties &props) {
            if (T <= T_critical())
            {
                props.fluidName = name();
                shared_ptr<CoolProp::AbstractState> water(CoolProp::AbstractState::factory("HEOS","Water"));
                water->update(CoolProp::QT_INPUTS, 0, T);
                props.Rho_l = water->saturated_liquid_keyed_output(CoolProp::iDmass);
                props.Rho_v = water->saturated_vapor_keyed_output(CoolProp::iDmass);
                props.H_l = water->saturated_liquid_keyed_output(CoolProp::iHmass);
                props.H_v = water->saturated_vapor_keyed_output(CoolProp::iHmass);
                props.Cp_l = water->saturated_liquid_keyed_output(CoolProp::iCpmass);
                props.Cp_v = water->saturated_vapor_keyed_output(CoolProp::iCpmass);
                props.IsothermalCompressibility_l = water->saturated_liquid_keyed_output(CoolProp::iisothermal_compressibility);
                props.IsothermalCompressibility_v = water->saturated_vapor_keyed_output(CoolProp::iisothermal_compressibility);
                props.IsobaricExpansivity_l = water->saturated_liquid_keyed_output(CoolProp::iisobaric_expansion_coefficient);
                props.IsobaricExpansivity_v = water->saturated_vapor_keyed_output(CoolProp::iisobaric_expansion_coefficient);
                props.Mu_l = water->saturated_liquid_keyed_output(CoolProp::iviscosity);
                props.Mu_v = water->saturated_vapor_keyed_output(CoolProp::iviscosity);
                props.p = water->saturated_liquid_keyed_output(CoolProp::iP);
                props.T = T;
                return  props.p;
            }
            return NAN;
        }

        double cIAPWS95_CoolProp::Boiling_T(const double &p, ThermodynamicProperties &props) {
            if (p <= p_critical())
            {
                shared_ptr<CoolProp::AbstractState> water(CoolProp::AbstractState::factory("HEOS","Water"));
                water->update(CoolProp::PQ_INPUTS, p, 0);

                props.Rho_l = water->saturated_liquid_keyed_output(CoolProp::iDmass);
                props.Rho_v = water->saturated_vapor_keyed_output(CoolProp::iDmass);
                props.H_l = water->saturated_liquid_keyed_output(CoolProp::iHmass);
                props.H_v = water->saturated_vapor_keyed_output(CoolProp::iHmass);
                props.Cp_l = water->saturated_liquid_keyed_output(CoolProp::iCpmass);
                props.Cp_v = water->saturated_vapor_keyed_output(CoolProp::iCpmass);
                props.Mu_l = water->saturated_liquid_keyed_output(CoolProp::iviscosity);
                props.Mu_v = water->saturated_vapor_keyed_output(CoolProp::iviscosity);
                props.p = water->saturated_liquid_keyed_output(CoolProp::iP);
                props.T = water->saturated_liquid_keyed_output(CoolProp::iT);;
                return  props.T;
            }
            return NAN;
        }
        /**
        * @brief Calculate phase for given T, p
        * @param T [K]
        * @param p [Pa]
        * @param X [default 0]
        * @return
        */
        PhaseRegion cIAPWS95_CoolProp::findPhaseRegion_TPX(const double &T, const double &p, const double &X) {
            CoolProp::phases phase_;
            if (T<m_constants.T_critical && p < m_constants.p_critical) //sub-critical region
            {
                shared_ptr<CoolProp::AbstractState> water(CoolProp::AbstractState::factory("HEOS","Water"));
                water->update(CoolProp::QT_INPUTS, 0, T);
                double p_boil_l = water->saturated_liquid_keyed_output(CoolProp::iP);
                double p_boil_v = water->saturated_vapor_keyed_output(CoolProp::iP);
                // see line 1983 in HelmholtzEOSMixtureBackend.cpp of CoolProp, if do not check this, CoolProp will throw error for T,p close to boiling curve. Keep consistent with CoolProp
                if(p > p_boil_l*(1e-6 + 1)) // two phase
                {
                    phase_ = CoolProp::iphase_liquid;
                } else if(p<p_boil_v*(1 - 1E-6))
                {
                    phase_ = CoolProp::iphase_gas;
                }
                else
                {
                    phase_ = CoolProp::iphase_twophase;
                }
            } else  //critical region
            {
                // see source code of PhaseSI in CoolProp.cpp (Line ~ 1081)
                std::size_t Phase_int = CoolProp::PropsSI("Phase","T", T, "P", p, "Water");
                phase_= static_cast<CoolProp::phases>(Phase_int);
            }
            return phase_CoolProp2xThermal(phase_);
        }

        ThermodynamicProperties cIAPWS95_CoolProp::UpdateState_TPX(const double &T, const double &p, const double &X) {
            ThermodynamicProperties props;
            UpdateState_TPX(props, T, p);
            return props;
        }
    };

};

#endif