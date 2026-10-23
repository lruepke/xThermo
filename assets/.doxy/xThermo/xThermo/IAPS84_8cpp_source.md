

# File IAPS84.cpp

[**File List**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**IAPS84.cpp**](IAPS84_8cpp.md)

[Go to the documentation of this file](IAPS84_8cpp.md)


```C++

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
        // LIB-37: PROST's viscos() returns 0.0 when the state lies outside the validity range of its viscosity
        // correlation (ThirdParties/PROST/iaps.c:329-332); a valid viscosity is never 0. Report "not available" as NaN
        // instead of a physically meaningless zero. All other properties of the state are unaffected.
        static double viscosity_or_nan(Prop* prop)
        {
            const double mu = viscos(prop);
            return mu == 0.0 ? NAN : mu;
        }

        // LIB-49: PROST's tp() returns its last Newton iterate with prop->error = 1 when the density iteration has not
        // converged within IMAX steps (ThirdParties/PROST/iaps.c, tp(); since 2022). Near the critical point the iterate
        // can be far off. Such a state is restarted from its last iterate (each restart continues the iteration for up
        // to IMAX steps; PROST itself is unchanged) and returned only if verified: PROST reports convergence, or the
        // density-equivalent residual |p(T, rho) - p| / (rho (dp/drho)_T) of PROST's explicit evaluation water_td is
        // <= 1e-9 (reachable where PROST's relative pressure tolerance is not, e.g. liquid at very low pressure); and
        // rho > 0, (dp/drho)_T > 0 and, below T_c, rho on the side of the phase PROST identified. Otherwise
        // SolutionError. States PROST reports as converged in the first call are not touched; states PROST does not
        // evaluate at all (prop->d == 0, outside valid_tp) are left as they are here.
        static bool density_residual_ok(const double& T, const double& P, Prop* prop)
        {
            if (!(prop->d > 0) || !(prop->dp->d_CT > 0)) return false;
            Prop* q = newProp('t', 'd', 0);
            water_td(T, prop->d, q);
            const double residual = fabs(q->p - P) / (prop->d * prop->dp->d_CT);
            freeProp(q);
            return residual <= 1e-9;
        }
        static bool stable_on_phase_side(const double& T, const double& P, Prop* prop)
        {
            if (!(prop->d > 0) || !(prop->dp->d_CT > 0)) return false;
            if (T >= TCRIT) return true;
            Prop* l = newProp('t', 'p', 0);
            Prop* v = newProp('t', 'p', 0);
            sat_t(T, l, v);
            const bool ok = P > v->p ? prop->d >= l->d * (1.0 - 1e-9) : prop->d <= v->d * (1.0 + 1e-9);
            freeProp(l);
            freeProp(v);
            return ok;
        }
        // LIB-51: in the near-critical band creg.t < T < T_c PROST's saturation (psatc) is an empirical approximation with
        // p_l != p_v, and region_tp reports p_v <= p <= p_l as saturated. Such a state is evaluated with the equation at
        // the requested T (not at T -+ 0.1 K): its mechanically stable roots on the liquid and vapour branches are found
        // by bisection on the fundamental equation (meta_td, phase ONE); of these, the one with the lower Gibbs energy is
        // the equation's stable state. The branches are separated at the spinodals ((dp/drho)_T = 0) around the density
        // of minimum (dp/drho)_T in 0.6-1.4 rho_c (the equation's loop is centred near 280-300 kg/m3, not at PROST's
        // approximate coexisting densities, which close to T_c do not bracket it). Returns a formatted Prop (caller frees)
        // and whether it is the liquid branch; SolutionError if no stable root is found.
        static void eos_at(const double& T, const double& d, double& p, double& dpdrho, double& g)
        {
            Prop* q = newProp('t', 'd', 1);
            q->phase = ONE;
            meta_td(T, d, q);
            p = q->p; dpdrho = q->dp->d_CT; g = q->g;
            freeProp(q);
        }
        static double bisect_density(const double& T, const double& P, double lo, double hi)   // p(lo) <= P <= p(hi)
        {
            double p, a, g;
            for (int i = 0; i < 200 && hi - lo > 1e-13 * hi; i++)
            {
                const double mid = 0.5 * (lo + hi);
                eos_at(T, mid, p, a, g);
                if (p < P) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        }
        static Prop* near_critical_band_state(const double& T, const double& P, bool& liquid)
        {
            Prop* l = newProp('t', 'p', 0);
            Prop* v = newProp('t', 'p', 0);
            sat_t(T, l, v);
            const double dl = l->d, dv = v->d;
            freeProp(l);
            freeProp(v);
            double p, a, g;
            // density of minimum (dp/drho)_T in [0.6, 1.4] rho_c (golden section)
            const double r = 0.5 * (std::sqrt(5.0) - 1.0);
            double lo = 0.6 * DCRIT, hi = 1.4 * DCRIT, c = hi - r * (hi - lo), e = lo + r * (hi - lo), fc, fe;
            eos_at(T, c, p, fc, g);
            eos_at(T, e, p, fe, g);
            for (int i = 0; i < 100 && hi - lo > 1e-9 * hi; i++)
            {
                if (fc < fe) { hi = e; e = c; fe = fc; c = hi - r * (hi - lo); eos_at(T, c, p, fc, g); }
                else { lo = c; c = e; fc = fe; e = lo + r * (hi - lo); eos_at(T, e, p, fe, g); }
            }
            const double d_min = 0.5 * (lo + hi);
            eos_at(T, d_min, p, a, g);
            // branch limits: the spinodals if the isotherm has a loop, otherwise d_min for both
            double sv = d_min, sl = d_min;
            if (a < 0)
            {
                double x = 1e-3, y = d_min;          // vapour spinodal: (dp/drho) > 0 at x, < 0 at y
                for (int i = 0; i < 200 && y - x > 1e-12 * y; i++) { const double m = 0.5 * (x + y); eos_at(T, m, p, a, g); if (a > 0) x = m; else y = m; }
                sv = x;
                x = d_min; y = 2.0 * dl;              // liquid spinodal: (dp/drho) < 0 at x, > 0 at y
                for (int i = 0; i < 200 && y - x > 1e-12 * y; i++) { const double m = 0.5 * (x + y); eos_at(T, m, p, a, g); if (a > 0) y = m; else x = m; }
                sl = y;
            }
            double best_g = 0, best_d = 0;
            bool found = false;
            // liquid branch [sl, 2 dl], vapour branch [1e-3, sv]: p increases with rho on each
            const double branch_lo[2] = {sl, 1e-3}, branch_hi[2] = {2.0 * dl, sv};
            for (int k = 0; k < 2; k++)
            {
                double plo, phi;
                eos_at(T, branch_lo[k], plo, a, g);
                eos_at(T, branch_hi[k], phi, a, g);
                if (!(plo <= P && P <= phi)) continue;
                const double d = bisect_density(T, P, branch_lo[k], branch_hi[k]);
                double pd;
                eos_at(T, d, pd, a, g);
                if (!(a > 0) || !(std::fabs(pd - P) <= 1e-8 * P + d * a * 1e-9)) continue;
                if (!found || g < best_g) { best_g = g; best_d = d; liquid = (k == 0); found = true; }
            }
            if (!found)
                throw SolutionError("cIAPS84: no stable root of the IAPS-84 equation in the near-critical band at T=" + std::to_string(T)
                                    + " K, p=" + std::to_string(P) + " Pa (LIB-51)");
            Prop* prop = newProp('t', 'd', 1);
            prop->phase = ONE;
            meta_td(T, best_d, prop);
            return prop;
        }

        static void water_tp_verified(const double& T, const double& P, const double& dp, Prop* prop)
        {
            water_tp(T, P, 0.0, dp, prop);
            // LIB-45: PROST evaluated nothing (outside its range of evaluation, e.g. the ice regions): no state
            if (prop->error && prop->phase != TWO && !(prop->d > 0))
                throw OutOfRangeError("cIAPS84: PROST does not evaluate (T, p) = (" + std::to_string(T) + " K, " + std::to_string(P)
                                      + " Pa): outside its range of evaluation (LIB-45)");
            if (!prop->error || prop->phase == TWO || !(prop->d > 0)) return;  // converged, saturated, or not evaluated
            for (int restart = 0; prop->error && !density_residual_ok(T, P, prop); restart++)
            {
                if (restart == 9 || !(prop->d > 0)) break;
                water_tp(T, P, prop->d, dp, prop);
            }
            if ((prop->error && !density_residual_ok(T, P, prop)) || !stable_on_phase_side(T, P, prop))
            {
                throw SolutionError("cIAPS84: PROST density iteration not converged (LIB-49) at T=" + std::to_string(T)
                                    + " K, p=" + std::to_string(P) + " Pa");
            }
        }

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
        void cIAPS84::UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& P, const double& X)
        {
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
            props.fluidName = name();
            double d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = newProp('t', 'p', 1);
            d = 0.0;
            water_tp_verified(T,P,dp,prop);
            // LIB-51: near-critical band (creg.t < T < T_c): the equation's stable state at the requested T
            bool band = false, band_liquid = false;
            if(prop->phase == TWO && T > creg.t)
            {
                freeProp(prop);
                prop = near_critical_band_state(T, P, band_liquid);
                band = true;
            }
            // Zhikui: 18 Jul, 2022, If the P,T in two phase region, although it is extremely rare, let's move a bit to a single phase side
            else if(prop->phase == TWO)
            {
                Prop* prop_l = newProp('t', 'p', 1);
                Prop* prop_v = newProp('t', 'p', 1);
                sat_p(P, prop_l, prop_v);
                freeProp(prop);
                prop = newProp('t', 'p', 1);
                if(T>=prop_v->T)
                {
                    water_tp_verified(T+0.1,P,dp,prop);
                }
                else
                {
                    water_tp_verified(T-0.1,P,dp,prop);
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
            props.Mu = viscosity_or_nan(prop);
            // LIB-51: in the band the label follows the branch of the selected root
            props.phase = band ? (band_liquid ? (P > p_critical() ? Supercritical_liquid : SinglePhase_L) : SinglePhase_V)
                               : phase_PROST2xThermal(T,P, prop);
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
                    props.dRhodP_v = props.dRhodP;  // LIB-47: the vapour slot (was written to _l)
                    props.dRhodT_v = props.dRhodT;
                    props.IsothermalCompressibility_v = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_v = props.IsobaricExpansivity;
                    break;
                default:
                    props.Rho_l = props.Rho;
                    props.H_l = props.H;
                    props.Mu_l = props.Mu;
                    props.Cp_l = props.Cp;
                    props.S_l = 1;
                    props.S_v = 0;
                    props.dRhodP_l = props.dRhodP;  // LIB-47: the liquid slot (was written to _v)
                    props.dRhodT_l = props.dRhodT;
                    props.IsothermalCompressibility_l = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_l = props.IsobaricExpansivity;
                    break;
            }
            freeProp(prop);
        }

        void cIAPS84::UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X)
        {
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
            props.fluidName = name();
            double t = 0.0, d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = newProp('p', 'h', 1);
            d = 0.0;
            water_ph(p,H,t, d,dp,dh,prop);
            // LIB-45: no state computed (PROST extrapolates invalid (p, h) from its bounds and then returns a state)
            if (prop->error && !(prop->d > 0))
                throw OutOfRangeError("cIAPS84: PROST computed no state for (p, h) = (" + std::to_string(p) + " Pa, " + std::to_string(H)
                                      + " J/kg) (LIB-45)");
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
            props.Mu = viscosity_or_nan(prop);
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
                    props.Mu_l = viscosity_or_nan(prop_l);
                    props.Mu_v = viscosity_or_nan(prop_v);
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
                    props.S_v = 1;  // LIB-48: saturation of the phase slot, as in UpdateState_TPX
                    props.S_l = 0;
                    props.dRhodP_v = props.dRhodP;  // LIB-47
                    props.dRhodT_v = props.dRhodT;
                    props.IsothermalCompressibility_v = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_v = props.IsobaricExpansivity;
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
                    props.S_l = 1;  // LIB-48
                    props.S_v = 0;
                    props.dRhodP_l = props.dRhodP;  // LIB-47
                    props.dRhodT_l = props.dRhodT;
                    props.IsothermalCompressibility_l = props.IsothermalCompressibility;
                    props.IsobaricExpansivity_l = props.IsobaricExpansivity;
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
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
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
            props.Mu_l = viscosity_or_nan(prop_l);
            props.Mu_v = viscosity_or_nan(prop_v);
            // LIB-52: intrinsic derivatives of each saturated phase at its coexisting density (PROST, phase ONE)
            props.dRhodP_l = 1.0/prop_l->dp->d_CT;
            props.dRhodT_l = -prop_l->dp->T_Cd * props.dRhodP_l;
            props.IsothermalCompressibility_l = 1.0/(props.Rho_l * prop_l->dp->d_CT);
            props.IsobaricExpansivity_l = prop_l->dp->T_Cd/prop_l->dp->d_CT/props.Rho_l;
            props.dRhodP_v = 1.0/prop_v->dp->d_CT;
            props.dRhodT_v = -prop_v->dp->T_Cd * props.dRhodP_v;
            props.IsothermalCompressibility_v = 1.0/(props.Rho_v * prop_v->dp->d_CT);
            props.IsobaricExpansivity_v = prop_v->dp->T_Cd/prop_v->dp->d_CT/props.Rho_v;
            freeProp(prop_l);
            freeProp(prop_v);
            return props.p;
        }

        double cIAPS84::Boiling_T(const double &p, ThermodynamicProperties &props) {
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
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
            props.Mu_l = viscosity_or_nan(prop_l);
            props.Mu_v = viscosity_or_nan(prop_v);
            props.p = p;
            props.T = prop_l->T;
            // LIB-52: intrinsic derivatives of each saturated phase at its coexisting density (PROST, phase ONE)
            props.dRhodP_l = 1.0/prop_l->dp->d_CT;
            props.dRhodT_l = -prop_l->dp->T_Cd * props.dRhodP_l;
            props.IsothermalCompressibility_l = 1.0/(props.Rho_l * prop_l->dp->d_CT);
            props.IsobaricExpansivity_l = prop_l->dp->T_Cd/prop_l->dp->d_CT/props.Rho_l;
            props.dRhodP_v = 1.0/prop_v->dp->d_CT;
            props.dRhodT_v = -prop_v->dp->T_Cd * props.dRhodP_v;
            props.IsothermalCompressibility_v = 1.0/(props.Rho_v * prop_v->dp->d_CT);
            props.IsobaricExpansivity_v = prop_v->dp->T_Cd/prop_v->dp->d_CT/props.Rho_v;
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
```


