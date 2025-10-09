/**
 * @file IAPS84.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of IAPS84.
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "IAPS84.h"
#include "steam4.h"
#include "iaps.h"
#include "LookUpTableForest.h"
// ----- LUT related head filess --------
#include "LookUpTableForestI.H"
#include "interpolationI.H"
#include "AMR_LUT_RefineFuncI.H"

//---------------------------------------

namespace xThermal
{
    namespace PROST
    {
        void cIAPS84::initialize_data()
        {
            m_constants.Tmin = 260;
            m_constants.Tmax = 2500;
            m_constants.pmin = PMIN;
            m_constants.pmax = 3000E6;
            m_constants.Ttriple = 273.16;
            m_constants.T_critical = TCRIT;
            m_constants.p_critical = PCRIT;
            m_constants.rhomass_critical = DCRIT;
            m_constants.molar_mass = 0.018015268; //kg/mol
        }
        cIAPS84::cIAPS84(/* args */)
        {
            initialize_data();
        }
        // copy constructor
        cIAPS84::cIAPS84(const cIAPS84& water)
        {
            initialize_data();
        }
        cIAPS84::~cIAPS84()
        {
        }
        /**
         * Calculate water properties based on IAPS84 and implementation of PROST.
         *
         * \warning The viscosity will be zero when T>900 deg.C, see line 332 of iaps.c in PROST.
         *
         * \warning The T,p validation in the matlab code (see Line 39-43 of water_tp_IAPS84.m) is not complete, see also PROST code (Line 1235 of iaps.c), so for e.g. T = -11 deg.C(T_star of viscosity correction for T=1 deg.C and X=0.7 kg/kg) p=100E5 Pa, the matlab code still give output, but it is in an invalid region.
         *
         * @param props
         * @param T [K]
         * @param P [Pa]
         * @param X [default 0 for H2O]
         */
        void cIAPS84::UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& P, const double& X)
        {
            props.fluidName = name();
            double d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = newProp('t', 'p', 1);
            d = 0.0;
            water_tp(T,P,d,dp,prop);
            // Zhikui: 18 Jul, 2022, If the P,T in two phase region, although it is extremely rare, let's move a bit to a single phase side
            if(prop->phase == TWO)
            {
                Prop* prop_l = newProp('t', 'p', 1);
                Prop* prop_v = newProp('t', 'p', 1);
                sat_p(P, prop_l, prop_v);
                freeProp(prop);
                prop = newProp('t', 'p', 1);
                if(T>=prop_v->T)
                {
                    water_tp(T+0.1,P,d,dp,prop);
                }
                else
                {
                    water_tp(T-0.1,P,d,dp,prop);
                }
                freeProp(prop_l);
                freeProp(prop_v);
            }
            // if (prop->error)
            // {
            //     char errorinfo[200];
            //     sprintf(errorinfo, "Error (from water_tp function of PROST) in UpdateState_TPX of IAPS84.cpp, T=%f K, P=%f Pa. But still try to return result.", T,P);
            //     WARNING(errorinfo);
            // }
            props.T = T;
            props.p = P;
            props.Rho = prop->d;
            props.H = prop->h;
            props.Cp = prop->cp;
            prop->T = T;
            prop->p = P;
            props.Mu = viscos(prop);
            props.phase = phase_PROST2xThermal(T,P, prop);
            // ============== derivatives ==============
            props.dRhodP = 1.0/prop->dp->d_CT;
            props.dRhodT = -prop->dp->T_Cd * props.dRhodP;
            props.IsothermalCompressibility = 1.0/(props.Rho * prop->dp->d_CT);
            props.IsobaricExpansivity = prop->dp->T_Cd/prop->dp->d_CT/props.Rho; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
            //==========================================
            switch (props.phase) {
                case SinglePhase_V:
                    props.Rho_v = props.Rho;
                    props.H_v = props.H;
                    props.Mu_v = props.Mu;
                    props.Cp_v = props.Cp;
                    props.S_v = 1;
                    props.S_l = 0;
                    props.dRhodP_l = props.dRhodP;
                    props.dRhodT_l = props.dRhodT;
                    props.IsothermalCompressibility_l = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_l = props.IsobaricExpansivity;
                    break;
                default:
                    props.Rho_l = props.Rho;
                    props.H_l = props.H;
                    props.Mu_l = props.Mu;
                    props.Cp_l = props.Cp;
                    props.S_l = 1;
                    props.S_v = 0;
                    props.dRhodP_v = props.dRhodP;
                    props.dRhodT_v = props.dRhodT;
                    props.IsothermalCompressibility_v = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_v = props.IsobaricExpansivity;
                    break;
            }
            freeProp(prop);
        }

        void cIAPS84::UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X)
        {
            props.fluidName = name();
            double t, d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = newProp('p', 'h', 1);
            d = 0.0;
            water_ph(p,H,t, d,dp,dh,prop);
            // if (prop->error)
            // {
            //     char errorinfo[200];
            //     sprintf(errorinfo, "Error (from water_tp function of PROST) in UpdateState_HPX of IAPS84.cpp, H=%f J/kg, P=%f Pa. But still try to return result.", H,p);
            //     WARNING(errorinfo);
            // }
            props.T = prop->T;
            props.p = p;
            props.Rho = prop->d;
            props.H = prop->h;
            props.Cp = prop->cp;
            prop->p = p;
            props.Mu = viscos(prop);
            props.phase = phase_PROST2xThermal(prop->T,p, prop);
            // ============== derivatives ==============
            props.dRhodP = 1.0/prop->dp->d_CT;
            props.dRhodT = -prop->dp->T_Cd * props.dRhodP;
            props.IsothermalCompressibility = 1.0/(props.Rho * prop->dp->d_CT);
            props.IsobaricExpansivity = prop->dp->T_Cd/prop->dp->d_CT/props.Rho; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
            //==========================================
            switch (props.phase)
            {
                case TwoPhase_VL_Water:
                {
                    Prop* prop_l = newProp('t', 'p', 1);
                    Prop* prop_v = newProp('t', 'p', 1);
                    sat_t(prop->T, prop_l, prop_v);
                    props.Rho_l = prop_l->d;
                    props.Rho_v = prop_v->d;
                    props.H_l = prop_l->h;
                    props.H_v = prop_v->h;
                    props.Cp_l = prop_l->cp;
                    props.Cp_v = prop_v->cp;
                    props.Mu_l = viscos(prop_l);
                    props.Mu_v = viscos(prop_v);
                    props.S_v = prop->x; // Steam quota: volume fraction
                    props.S_l = 1.0 - prop->x;
                    // ============== derivatives ==============
                    props.dRhodP_l = 1.0/prop_l->dp->d_CT;
                    props.dRhodT_l = -prop_l->dp->T_Cd * props.dRhodP_l;
                    props.IsothermalCompressibility_l = 1.0/(props.Rho_l * prop_l->dp->d_CT);
                    props.IsobaricExpansivity_l = prop_l->dp->T_Cd/prop_l->dp->d_CT/props.Rho_l; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
                    props.dRhodP_v = 1.0/prop_v->dp->d_CT;
                    props.dRhodT_v = -prop_v->dp->T_Cd * props.dRhodP_v;
                    props.IsothermalCompressibility_v = 1.0/(props.Rho_v * prop_v->dp->d_CT);
                    props.IsobaricExpansivity_v = prop_v->dp->T_Cd/prop_v->dp->d_CT/props.Rho_v; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
                    //==========================================
                    freeProp(prop_l);
                    freeProp(prop_v);
                }
                break;
                case SinglePhase_V:
                    props.Rho_v = props.Rho;
                    props.H_v = props.H;
                    props.Mu_v = props.Mu;
                    props.Cp_v = props.Cp;
                    break;
                case Supercritical:
                case Supercritical_liquid:
                case Supercritical_vapor:
                case SinglePhase_L:
                case Critical:
                    props.Rho_l = props.Rho;
                    props.H_l = props.H;
                    props.Mu_l = props.Mu;
                    props.Cp_l = props.Cp;
                    break;
                default:
                    throw NotImplementedError("The phase is unsupported in UpdateState_HPX function of "+name()
                                              +", T="+std::to_string(props.T)+", p="+std::to_string(p)+", h="+std::to_string(H)+", phase: "+ phase_name(props.phase)+", PROST phase: "+std::to_string(prop->phase));
                    break;
            }
            freeProp(prop);
        }

        double cIAPS84::Boiling_p(const double& T)
        {
            Prop* prop_l = newProp('t', 'p', 0);
            Prop* prop_v = newProp('t', 'p', 0);
            sat_t(T, prop_l, prop_v);
            double p_sat = prop_v->p;

            freeProp(prop_l);
            freeProp(prop_v);
            return p_sat;
        }
        
        double cIAPS84::Boiling_T(const double& p)
        {
            Prop* prop_l = newProp('t', 'p', 0);
            Prop* prop_v = newProp('t', 'p', 0);
            sat_p(p, prop_l, prop_v);
            
            double T_sat = prop_v->T;

            freeProp(prop_l);
            freeProp(prop_v);
            return T_sat;
        }
        
        double cIAPS84::Boiling_T(const double& p, double& rho_l, double& rho_v)
        {
            Prop* prop_l = newProp('t', 'p', 0);
            Prop* prop_v = newProp('t', 'p', 0);
            sat_p(p, prop_l, prop_v);
            rho_l = prop_l->d;
            rho_v = prop_v->d;
            double T_boil = prop_l->T;
            freeProp(prop_l);
            freeProp(prop_v);
            return T_boil;
        }

        double cIAPS84::Boiling_p(const double& T, double& rho_l, double& rho_v)
        {
            Prop* prop_l = newProp('t', 'p', 0);
            Prop* prop_v = newProp('t', 'p', 0);
            sat_t(T, prop_l, prop_v);
            rho_l = prop_l->d;
            rho_v = prop_v->d;
            double p_boil = prop_l->p;
            freeProp(prop_l);
            freeProp(prop_v);
            return p_boil;
        }
        /**
         * @brief Get phase
         * Because the phase index in PROST only has option of ONE and TWO.
         * So the detailed phase region is implemented below.
         * @return xThermal::PhaseRegion 
         */
        xThermal::PhaseRegion cIAPS84::phase_PROST2xThermal(const double& T, const double& p, const void * prop)
        {
            if (((Prop*)prop)->phase==ONE)
            {
                if(p==p_critical() && T == T_critical())
                {
                    return Critical;
                }else
                {
                    if (p > p_critical())
                    {
                        if(T < T_critical())
                        return Supercritical_liquid;
                        else return Supercritical;
                    }else
                    {
                        if(T > T_critical())
                        {
                            return SinglePhase_V;
                        }else
                        {
                            // calculate saturated pressure
                            double p_boil = Boiling_p(T);
                            if(p > p_boil)
                            {
                                return SinglePhase_L;
                            }else
                            {
                                return SinglePhase_V;
                            }
                        }
                    }
                }
            }else if (((Prop*)prop)->phase == TWO)
            {
                return TwoPhase_VL_Water;
            }
            return Unknown;
        }

        double cIAPS84::Boiling_p(const double &T, ThermodynamicProperties &props) {
            props.fluidName = name();
            Prop* prop_l = newProp('t', 'p', 1);
            Prop* prop_v = newProp('t', 'p', 1);
            sat_t(T, prop_l, prop_v);
            props.Rho_l = prop_l->d;
            props.Rho_v = prop_v->d;
            props.H_l = prop_l->h;
            props.H_v = prop_v->h;
            props.Cp_l = prop_l->cp;
            props.Cp_v = prop_v->cp;
            props.p = prop_l->p;
            props.T = T;
            props.Mu_l = viscos(prop_l);
            props.Mu_v = viscos(prop_v);
            freeProp(prop_l);
            freeProp(prop_v);
            return props.p;
        }

        double cIAPS84::Boiling_T(const double &p, ThermodynamicProperties &props) {
            props.fluidName = name();
            Prop* prop_l = newProp('t', 'p', 1);
            Prop* prop_v = newProp('t', 'p', 1);
            sat_p(p, prop_l, prop_v);
            props.Rho_l = prop_l->d;
            props.Rho_v = prop_v->d;
            props.H_l = prop_l->h;
            props.H_v = prop_v->h;
            props.Cp_l = prop_l->cp;
            props.Cp_v = prop_v->cp;
            props.Mu_l = viscos(prop_l);
            props.Mu_v = viscos(prop_v);
            props.p = p;
            props.T = prop_l->T;
            freeProp(prop_l);
            freeProp(prop_v);
            return props.T;
        }

        PhaseRegion cIAPS84::findPhaseRegion_TPX(const double &T, const double &p, const double &X) {
            double d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = newProp('t', 'p', 0);
            d = 0.0;
            water_tp(T,p,d,dp,prop);
            PhaseRegion phase = phase_PROST2xThermal(T, p, prop);
            freeProp(prop);

            return phase;
        }
        /**
        * @brief Calculate water properties for given T, p
        * @param stateArray
        * @param T
        * @param p
        * @param N
        */
        // void cIAPS84::UpdateState_TPX(ThermodynamicPropertiesArray &stateArray, const size_t &N, const double *T, const double *p, const double* X) {
        //     ThermodynamicProperties props;
        //     stateArray.N = N;
        //     for (size_t i = 0; i < N; ++i) {
        //         stateArray.T[i] = T[i];
        //         stateArray.p[i] = p[i];
        //         UpdateState_TPX(props, stateArray.T[i], stateArray.p[i]);
        //         stateArray.Rho[i] = props.Rho;
        //         stateArray.H[i] = props.H;
        //         stateArray.Cp[i] = props.Cp;
        //         stateArray.Mu[i] = props.Mu;
        //         stateArray.phase[i]=props.phase;
        //     }
        // }

        ThermodynamicProperties cIAPS84::UpdateState_TPX(const double &T, const double &p, const double &X) {
            ThermodynamicProperties props;
            UpdateState_TPX(props, T, p);
            return props;
        }

    };
};