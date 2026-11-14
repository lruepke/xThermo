/**
 * @file IAPS84.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of the xThermo backend "IAPS84", a C++ wrapper around PROST 4.1 (see IAPS84.h).
 *
 * All IAPS-84 computations (equation of state, derivatives, saturation, viscosity) are done by PROST 4.1
 * (O. Bauer and O. Engel, TU Hamburg-Harburg, 1998; ThirdParties/PROST; GNU LGPL v2.0). This file contains only the
 * xThermo interface around it.
 * @version 0.1
 * @date 2022-03-27
 *
 * @copyright Copyright (c) 2022. This wrapper is part of xThermo (MIT licence, see LICENSE). PROST is not covered by
 *            xThermo's licence (see THIRD_PARTY_NOTICES.md).
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
        // G2 (review/G2_IAPS84_SAT.md): the isotherm of the equation near T_c (LIB-51 arithmetic, unchanged): density of
        // minimum (dp/drho)_T in [0.6, 1.4] rho_c, and the spinodals if the isotherm has a loop (otherwise d_min for both);
        // the liquid branch is [sl, 2 dl], the vapour branch [1e-3, sv] (dl: PROST's saturated liquid at T, sat_t; above
        // PROST's T_c, where sat_t is not defined, dl = rho_c)
        struct Isotherm { double d_min, sv, sl, d_hi, p_sv, p_sl; bool loop; };
        static Isotherm eq_isotherm(const double& T)
        {
            double dl = DCRIT;
            if (T < TCRIT)
            {
                Prop* l = newProp('t', 'p', 0);
                Prop* v = newProp('t', 'p', 0);
                sat_t(T, l, v);
                dl = l->d;
                freeProp(l);
                freeProp(v);
            }
            double p, a, g;
            const double r = 0.5 * (std::sqrt(5.0) - 1.0);
            double lo = 0.6 * DCRIT, hi = 1.4 * DCRIT, c = hi - r * (hi - lo), e = lo + r * (hi - lo), fc, fe;
            eos_at(T, c, p, fc, g);
            eos_at(T, e, p, fe, g);
            for (int i = 0; i < 100 && hi - lo > 1e-9 * hi; i++)
            {
                if (fc < fe) { hi = e; e = c; fe = fc; c = hi - r * (hi - lo); eos_at(T, c, p, fc, g); }
                else { lo = c; c = e; fc = fe; e = lo + r * (hi - lo); eos_at(T, e, p, fe, g); }
            }
            Isotherm iso;
            iso.d_min = 0.5 * (lo + hi);
            iso.d_hi = 2.0 * dl;
            eos_at(T, iso.d_min, p, a, g);
            iso.loop = a < 0;
            iso.sv = iso.d_min; iso.sl = iso.d_min;
            if (iso.loop)
            {
                double x = 1e-3, y = iso.d_min;          // vapour spinodal: (dp/drho) > 0 at x, < 0 at y
                for (int i = 0; i < 200 && y - x > 1e-12 * y; i++) { const double m = 0.5 * (x + y); eos_at(T, m, p, a, g); if (a > 0) x = m; else y = m; }
                iso.sv = x;
                x = iso.d_min; y = iso.d_hi;               // liquid spinodal: (dp/drho) < 0 at x, > 0 at y
                for (int i = 0; i < 200 && y - x > 1e-12 * y; i++) { const double m = 0.5 * (x + y); eos_at(T, m, p, a, g); if (a > 0) y = m; else x = m; }
                iso.sl = y;
            }
            eos_at(T, iso.sv, iso.p_sv, a, g);
            eos_at(T, iso.sl, iso.p_sl, a, g);
            return iso;
        }
        // the mechanically stable root of the equation at (T, P) on branch k (0 liquid, 1 vapour), if it exists
        static bool eq_branch_root(const double& T, const double& P, const Isotherm& iso, int k, double& d, double& g)
        {
            const double lo = k == 0 ? iso.sl : 1e-3, hi = k == 0 ? iso.d_hi : iso.sv;
            double plo, phi, a, pd;
            eos_at(T, lo, plo, a, g);
            eos_at(T, hi, phi, a, g);
            if (!(plo <= P && P <= phi)) return false;
            d = bisect_density(T, P, lo, hi);
            eos_at(T, d, pd, a, g);
            return (a > 0) && std::fabs(pd - P) <= 1e-8 * P + d * a * 1e-9;
        }
        static Prop* eq_prop(const double& T, const double& d)    // formatted state of the equation at (T, d), phase ONE
        {
            Prop* prop = newProp('t', 'd', 1);
            prop->phase = ONE;
            meta_td(T, d, prop);
            return prop;
        }
        // LIB-51 / G2: the stable state of the equation at (T, P): of the roots on both branches the one with the lower
        // Gibbs energy. Returns a formatted Prop (caller frees) and whether it is the liquid branch.
        static Prop* near_critical_band_state(const double& T, const double& P, bool& liquid, const Isotherm& iso)
        {
            double best_g = 0, best_d = 0;
            bool found = false;
            for (int k = 0; k < 2; k++)
            {
                double d, g;
                if (!eq_branch_root(T, P, iso, k, d, g)) continue;
                if (!found || g < best_g) { best_g = g; best_d = d; liquid = (k == 0); found = true; }
            }
            if (!found)
                throw SolutionError("cIAPS84: no stable root of the IAPS-84 equation in the near-critical band at T=" + std::to_string(T)
                                    + " K, p=" + std::to_string(P) + " Pa (LIB-51)");
            return eq_prop(T, best_d);
        }
        // G2 (endpoint): the critical temperature of the IAPS-84 equation itself, the highest T whose isotherm has a loop
        // (min (dp/drho)_T < 0): 647.1264452 K (critical pressure 22.054008 MPa, density 296.8 kg/m3), 0.45 mK above and
        // 0.9 kPa below PROST's declared critical point (647.126 K, 22.054915 MPa). The equation's coexistence ends there.
        static double eq_critical_T()
        {
            static const double Tc_eq = [] {
                double a = TCRIT - 0.01, b = TCRIT + 0.01;           // loop at a, none at b
                for (int i = 0; i < 100 && b - a > 4 * (std::nextafter(b, HUGE_VAL) - b); i++)
                {
                    const double m = 0.5 * (a + b);
                    if (eq_isotherm(m).loop) a = m; else b = m;
                }
                return a;
            }();
            return Tc_eq;
        }
        // G2: the coexistence of the equation at T (creg.t < T <= T_c of the equation): the pressure where the liquid and vapour roots
        // have equal Gibbs energy (Maxwell), between the spinodal pressures; g_l - g_v decreases with p (v_l < v_v).
        // Regula falsi (Illinois) with bisection safeguard. False if the isotherm has no loop.
        static bool eq_coexistence_T(const double& T, double& p_sat, double& dl, double& dv)
        {
            thread_local double last_T = NAN, last_p, last_dl, last_dv;
            thread_local bool last_ok = false;
            if (T == last_T) { p_sat = last_p; dl = last_dl; dv = last_dv; return last_ok; }
            const Isotherm iso = eq_isotherm(T);
            bool ok = false;
            if (iso.loop)
            {
                double gl, gv, plo_v;
                {   double a, g; eos_at(T, 1e-3, plo_v, a, g); }
                // both roots exist for p_sl <= p <= p_sv; the ends are taken 1e-6 of the width inside, where (dp/drho)_T > 0
                // at the roots (the coexistence pressure lies near the middle; near T_c the width is a few Pa)
                const double lo0 = std::max(iso.p_sl, plo_v), w0 = iso.p_sv - lo0;
                double a = lo0 + 1e-6 * w0, b = iso.p_sv - 1e-6 * w0;
                auto dg = [&](double P, double& d_l, double& d_v, bool& valid) {
                    valid = eq_branch_root(T, P, iso, 0, d_l, gl) && eq_branch_root(T, P, iso, 1, d_v, gv);
                    return gl - gv;
                };
                bool va, vb, vm;
                double dla, dva, dlb, dvb, dlm = NAN, dvm = NAN;
                double fa = dg(a, dla, dva, va), fb = dg(b, dlb, dvb, vb);
                if (va && vb && fa > 0 && fb < 0)
                {
                    int side = 0;
                    double m = a;
                    for (int i = 0; i < 200 && b - a > 1e-15 * b; i++)
                    {
                        m = (a * fb - b * fa) / (fb - fa);
                        if (!(m > a && m < b) || i % 8 == 7) m = 0.5 * (a + b);   // safeguard
                        double fm = dg(m, dlm, dvm, vm);
                        if (!vm && m != 0.5 * (a + b)) { m = 0.5 * (a + b); fm = dg(m, dlm, dvm, vm); }   // retry by bisection
                        if (!vm) break;
                        if (fm == 0) { a = b = m; break; }
                        if (fm > 0) { a = m; fa = fm; if (side == 1) fb *= 0.5; side = 1; }
                        else { b = m; fb = fm; if (side == -1) fa *= 0.5; side = -1; }
                    }
                    p_sat = 0.5 * (a + b);
                    ok = eq_branch_root(T, p_sat, iso, 0, dl, gl) && eq_branch_root(T, p_sat, iso, 1, dv, gv);
                }
                // within ~3e-5 K of the equation's critical temperature g_l - g_v is below the rounding of g (~1e-9 J/kg)
                // and its sign is noise; the loop itself spans < 1e-3 Pa there, so every pressure inside it is the
                // coexistence pressure to that accuracy: the middle of the loop, with the branch roots there
                if (!ok && w0 <= 1e-9 * iso.p_sv)
                {
                    p_sat = 0.5 * (a + b);
                    ok = eq_branch_root(T, p_sat, iso, 0, dl, gl) && eq_branch_root(T, p_sat, iso, 1, dv, gv);
                }
            }
            last_T = T; last_p = p_sat; last_dl = dl; last_dv = dv; last_ok = ok;
            return ok;
        }
        // G2: the coexistence temperature of the equation at P: p_sat(T) of eq_coexistence_T increases with T; Illinois
        // on T within [creg.t - 0.05 K, T_c,eq] (the equation's critical temperature, eq_critical_T). Above the
        // equation's critical pressure (22.054008 MPa, 0.9 kPa below PROST's p_crit) the equation has no coexistence and
        // T-P-X is single-phase: the saturation is degenerate there, T_sat = T_c,eq and both phases the equation's stable
        // state at (T_c,eq, P) (zero latent heat, continuous with the coexistence at its endpoint).
        static bool eq_coexistence_P(const double& P, double& T_sat, double& dl, double& dv)
        {
            thread_local double last_P = NAN, last_T, last_dl, last_dv;
            thread_local bool last_ok = false;
            if (P == last_P) { T_sat = last_T; dl = last_dl; dv = last_dv; return last_ok; }
            bool ok = false;
            // the end of the coexistence curve: the highest T <= T_c,eq at which the coexistence is solvable (the loop
            // vanishes at T_c,eq itself), computed once
            static const double T_end = [] {
                double d_l, d_v, p_s;
                const double Tc = eq_critical_T();
                if (eq_coexistence_T(Tc, p_s, d_l, d_v)) return Tc;
                double ok = TCRIT, bad = Tc;                 // solvable at PROST's T_c (verified), not at T_c,eq
                for (int i = 0; i < 100 && bad - ok > 4 * (std::nextafter(bad, HUGE_VAL) - bad); i++)
                {
                    const double m = 0.5 * (ok + bad);
                    if (eq_coexistence_T(m, p_s, d_l, d_v)) ok = m; else bad = m;
                }
                return ok;
            }();
            double ps, a = creg.t - 0.05, b = T_end, fa, fb;
            if (eq_coexistence_T(b, ps, dl, dv) && ps < P)
            {
                bool liquid;
                Prop* c = near_critical_band_state(b, P, liquid, eq_isotherm(b));
                T_sat = b; dl = dv = c->d;
                freeProp(c);
                ok = true;
            }
            else if (eq_coexistence_T(b, ps, dl, dv) && (fb = ps - P) >= 0)
            {
                if (fb == 0) { T_sat = b; ok = true; }
                else if (eq_coexistence_T(a, ps, dl, dv) && (fa = ps - P) < 0)
                {
                    int side = 0;
                    for (int i = 0; i < 200 && b - a > 4 * (std::nextafter(b, HUGE_VAL) - b); i++)
                    {
                        double m = (a * fb - b * fa) / (fb - fa);
                        if (!(m > a && m < b) || i % 8 == 7) m = 0.5 * (a + b);
                        if (!eq_coexistence_T(m, ps, dl, dv))
                        {
                            // within ~1e-5 K of T_c,eq the Gibbs-energy difference reaches the rounding of g (~1e-9 J/kg):
                            // an unresolvable T is treated as above the coexistence temperature (bisection from above)
                            b = m; fb = std::fabs(fb); side = 0;
                            continue;
                        }
                        const double fm = ps - P;
                        if (fm == 0) { a = b = m; break; }
                        if (fm < 0) { a = m; fa = fm; if (side == 1) fb *= 0.5; side = 1; }
                        else { b = m; fb = fm; if (side == -1) fa *= 0.5; side = -1; }
                    }
                    T_sat = std::fabs(fa) < std::fabs(fb) ? a : b;
                    ok = eq_coexistence_T(T_sat, ps, dl, dv);
                }
            }
            last_P = P; last_T = T_sat; last_dl = dl; last_dv = dv; last_ok = ok;
            return ok;
        }
        // G2: the saturation of the backend above creg uses the equation's coexistence (below creg PROST's psat/tsat,
        // the equation's Maxwell solution)
        static bool use_eq_saturation_p(const double& p) { return p > creg.p * 1.0e6 && p < PCRIT; }
        static bool use_eq_saturation_T(const double& T) { return T > creg.t && T <= eq_critical_T(); }
        // fill the saturated-phase fields of props from PROST Props of the liquid and vapour (LIB-52 derivatives)
        static void fill_saturated(ThermodynamicProperties& props, Prop* prop_l, Prop* prop_v)
        {
            props.Rho_l = prop_l->d;
            props.Rho_v = prop_v->d;
            props.H_l = prop_l->h;
            props.H_v = prop_v->h;
            props.Cp_l = prop_l->cp;
            props.Cp_v = prop_v->cp;
            props.Mu_l = viscosity_or_nan(prop_l);
            props.Mu_v = viscosity_or_nan(prop_v);
            props.dRhodP_l = 1.0/prop_l->dp->d_CT;
            props.dRhodT_l = -prop_l->dp->T_Cd * props.dRhodP_l;
            props.IsothermalCompressibility_l = 1.0/(props.Rho_l * prop_l->dp->d_CT);
            props.IsobaricExpansivity_l = prop_l->dp->T_Cd/prop_l->dp->d_CT/props.Rho_l;
            props.dRhodP_v = 1.0/prop_v->dp->d_CT;
            props.dRhodT_v = -prop_v->dp->T_Cd * props.dRhodP_v;
            props.IsothermalCompressibility_v = 1.0/(props.Rho_v * prop_v->dp->d_CT);
            props.IsobaricExpansivity_v = prop_v->dp->T_Cd/prop_v->dp->d_CT/props.Rho_v;
        }
        // G2 (XT-2): in PROST's saturated band below creg (|p/p_s(T) - 1| <= 1e-6) the root of the equation at the requested
        // T on the side of T_sat (sat_p): bisection on the density in a bracket around the saturated density at T (the root
        // is within ~1e-6 of it; (dp/drho)_T > 0 required). Replaces the former state of T -+ 0.1 K.
        static Prop* saturated_band_state(const double& T, const double& P, bool& liquid)
        {
            Prop* l = newProp('t', 'p', 0);
            Prop* v = newProp('t', 'p', 0);
            sat_p(P, l, v);
            liquid = T < v->T;
            sat_t(T, l, v);
            const double ds = liquid ? l->d : v->d;
            freeProp(l);
            freeProp(v);
            double lo = ds, hi = ds, plo, phi, a, g;
            eos_at(T, ds, plo, a, g);
            phi = plo;
            for (double eps = 1e-12; !(plo <= P && P <= phi) && eps < 0.05; eps *= 4)
            {
                lo = ds * (1.0 - eps); hi = ds * (1.0 + eps);
                double alo, ahi;
                eos_at(T, lo, plo, alo, g);
                eos_at(T, hi, phi, ahi, g);
                if (!(alo > 0 && ahi > 0)) break;
            }
            if (!(plo <= P && P <= phi))
                throw SolutionError("cIAPS84: no root of the IAPS-84 equation on the " + std::string(liquid ? "liquid" : "vapour")
                                    + " branch next to saturation at T=" + std::to_string(T) + " K, p=" + std::to_string(P) + " Pa (G2)");
            const double d = plo == P ? lo : (phi == P ? hi : bisect_density(T, P, lo, hi));
            return eq_prop(T, d);
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
        void cIAPS84::UpdateState_TPX(ThermodynamicProperties& props, const double& T_in, const double& P_in, const double& X_in)
        {
            const double T = T_in, P = P_in, X = X_in;  // copies: the inputs may be fields of props, which is reset next
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
            if (outer.is_outer()) check_state_input("cIAPS84::UpdateState_TPX", "T", T, P);  // C2
            props.fluidName = name();
            double d, dp, ds, dh;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            Prop *prop = nullptr;
            d = 0.0;
            // LIB-51 / G2 (XT-12): near T_c (creg.t < T <= T_c) the equation's stable state at the requested T: of the
            // mechanically stable roots on the liquid and vapour branches the one with the lower Gibbs energy. Its
            // liquid/vapour switch is the equation's coexistence, which the saturation API uses above creg as well. PROST's
            // own density iteration (which near T_c may converge to the metastable root or not at all) is used only where
            // neither branch has a root (pressures beyond the branch brackets)
            bool band = false, band_liquid = false;
            if (use_eq_saturation_T(T))
            {
                const Isotherm iso = eq_isotherm(T);
                try { prop = near_critical_band_state(T, P, band_liquid, iso); band = true; }
                catch (const SolutionError&) { prop = nullptr; }
            }
            if (!band)
            {
                prop = newProp('t', 'p', 1);
                water_tp_verified(T,P,dp,prop);
                if (prop->phase == TWO && T > creg.t)
                {
                    freeProp(prop);
                    throw SolutionError("cIAPS84: no stable root of the IAPS-84 equation in the near-critical band at T=" + std::to_string(T)
                                        + " K, p=" + std::to_string(P) + " Pa (LIB-51)");
                }
                // G2 (XT-2): PROST's saturated band below creg: the equation's root at the requested T on the side of
                // T_sat (formerly the state of T -+ 0.1 K, "Zhikui 18 Jul 2022", reported with the label and T of the input)
                if (prop->phase == TWO)
                {
                    freeProp(prop);
                    prop = saturated_band_state(T, P, band_liquid);
                    band = true;
                }
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
            // G2: at T = T_c the labels of phase_PROST2xThermal (critical point; supercritical above p_c)
            props.phase = band ? ((P == p_critical() && T == T_critical()) ? Critical
                                  : band_liquid ? (P > p_critical() ? (T < T_critical() ? Supercritical_liquid : Supercritical) : SinglePhase_L)
                                                : SinglePhase_V)
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

        void cIAPS84::UpdateState_HPX(ThermodynamicProperties& props, const double& H_in, const double& p_in, const double& X_in)
        {
            const double H = H_in, p = p_in, X = X_in;  // copies: the inputs may be fields of props, which is reset next
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
            if (outer.is_outer()) check_state_input("cIAPS84::UpdateState_HPX", "H", H, p);  // C2
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
            // G2 (XT-12, review/G2_IAPS84_SAT.md): near saturation (PROST two-phase) and above creg the V+L state and the
            // single-phase branch next to it follow the backend's saturation (Boiling_T(p, props): PROST's sat_p below creg,
            // the equation's coexistence above it), so that T, the phase properties and the quality come from one
            // calculation and agree with UpdateState_TPX
            if (p < PCRIT && (prop->phase == TWO || use_eq_saturation_p(p)))
            {
                ThermodynamicProperties sat;
                Boiling_T(p, sat);
                if (H >= sat.H_l && H <= sat.H_v)
                {
                    const double x = (H - sat.H_l) / (sat.H_v - sat.H_l);       // vapour mass fraction (quality)
                    const double rho = 1.0 / ((1.0 - x) / sat.Rho_l + x / sat.Rho_v);
                    Prop* bulk = prop;
                    if (prop->phase != TWO)
                    {
                        // PROST's two-phase conventions for the bulk fields (cp = 0 "infinite", (dp/drho)_T = 0, viscosity)
                        bulk = newProp('t', 'd', 1);
                        bulk->phase = TWO;
                        meta_td(sat.T, rho, bulk);
                    }
                    props.T = sat.T;
                    props.p = p;
                    props.Rho = rho;
                    props.H = sat.H_l + x * (sat.H_v - sat.H_l);
                    props.Cp = bulk->cp;
                    props.Mu = viscosity_or_nan(bulk);
                    props.phase = TwoPhase_VL_Water;
                    props.dRhodP = 1.0/bulk->dp->d_CT;
                    props.dRhodT = -bulk->dp->T_Cd * props.dRhodP;
                    props.IsothermalCompressibility = 1.0/(props.Rho * bulk->dp->d_CT);
                    props.IsobaricExpansivity = bulk->dp->T_Cd/bulk->dp->d_CT/props.Rho;
                    props.Rho_l = sat.Rho_l; props.Rho_v = sat.Rho_v;
                    props.H_l = sat.H_l; props.H_v = sat.H_v;
                    props.Cp_l = sat.Cp_l; props.Cp_v = sat.Cp_v;
                    props.Mu_l = sat.Mu_l; props.Mu_v = sat.Mu_v;
                    props.S_v = x;              // the water models report the vapour mass fraction (docs/limitations.md)
                    props.S_l = 1.0 - x;
                    props.dRhodP_l = sat.dRhodP_l; props.dRhodT_l = sat.dRhodT_l;
                    props.IsothermalCompressibility_l = sat.IsothermalCompressibility_l;
                    props.IsobaricExpansivity_l = sat.IsobaricExpansivity_l;
                    props.dRhodP_v = sat.dRhodP_v; props.dRhodT_v = sat.dRhodT_v;
                    props.IsothermalCompressibility_v = sat.IsothermalCompressibility_v;
                    props.IsobaricExpansivity_v = sat.IsobaricExpansivity_v;
                    if (bulk != prop) freeProp(bulk);
                    freeProp(prop);
                    check_hpx_enthalpy("cIAPS84::UpdateState_HPX", props.H, H, props.T, props.phase);
                    return;
                }
                const bool liquid = H < sat.H_l;
                if (prop->phase == TWO || (liquid ? prop->T > sat.T : prop->T < sat.T))
                {
                    // the equation's state at (p, H) on the branch of the side, from the saturated state of that side
                    freeProp(prop);
                    prop = newProp('p', 'h', 1);
                    prop->phase = ONE;
                    meta_ph(p, H, sat.T, liquid ? sat.Rho_l : sat.Rho_v, dp, dh, prop);
                    if (prop->error || !(prop->d > 0))
                    {
                        freeProp(prop);
                        throw SolutionError("cIAPS84: no single-phase state of the IAPS-84 equation at (p, h) = (" + std::to_string(p)
                                            + " Pa, " + std::to_string(H) + " J/kg) next to saturation (G2)");
                    }
                }
            }
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
            // HPX-CHK (review/HPX_CHK.md): props.H is PROST's enthalpy of the solved state
            check_hpx_enthalpy("cIAPS84::UpdateState_HPX", props.H, H, props.T, props.phase);
        }

        double cIAPS84::Boiling_p(const double& T)
        {
            double p_eq, dl_eq, dv_eq;     // G2: above creg the equation's coexistence
            if (use_eq_saturation_T(T) && eq_coexistence_T(T, p_eq, dl_eq, dv_eq)) return p_eq;
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
            double T_eq, dl_eq, dv_eq;     // G2: above creg the equation's coexistence
            if (use_eq_saturation_p(p) && eq_coexistence_P(p, T_eq, dl_eq, dv_eq)) return T_eq;
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
            double T_eq;                   // G2: above creg the equation's coexistence
            if (use_eq_saturation_p(p) && eq_coexistence_P(p, T_eq, rho_l, rho_v)) return T_eq;
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
            double p_eq;                   // G2: above creg the equation's coexistence
            if (use_eq_saturation_T(T) && eq_coexistence_T(T, p_eq, rho_l, rho_v)) return p_eq;
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
            OuterEvaluation outer(props);  // C7: reset at the outermost entry
            props.fluidName = name();
            double p_eq, dl_eq, dv_eq;     // G2: above creg the equation's coexistence
            if (use_eq_saturation_T(T) && eq_coexistence_T(T, p_eq, dl_eq, dv_eq))
            {
                Prop* l = eq_prop(T, dl_eq);
                Prop* v = eq_prop(T, dv_eq);
                fill_saturated(props, l, v);
                freeProp(l);
                freeProp(v);
                props.p = p_eq;
                props.T = T;
                return props.p;
            }
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
            double T_eq, dl_eq, dv_eq;     // G2: above creg the equation's coexistence
            if (use_eq_saturation_p(p) && eq_coexistence_P(p, T_eq, dl_eq, dv_eq))
            {
                Prop* l = eq_prop(T_eq, dl_eq);
                Prop* v = eq_prop(T_eq, dv_eq);
                fill_saturated(props, l, v);
                freeProp(l);
                freeProp(v);
                props.p = p;
                props.T = T_eq;
                return props.T;
            }
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
            // G2: the label of UpdateState_TPX (the side of the equation's root near saturation, review/G2_IAPS84_SAT.md);
            // formerly PROST's region test alone, which reports saturation within |p/p_s - 1| <= 1e-6
            ThermodynamicProperties props;
            props.fluidName = name();
            UpdateState_TPX(props, T, p, X);
            return props.phase;
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