/**
 * @file IAPWS-IF97.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 
 * @version 0.1
 * @date 2022-04-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "IAPWS-IF97.h"

namespace xThermal
{
    namespace IAPWS_IF97
    {
        void cIAPWS_IF97::initialize_data()
        {
            m_constants.R = 461.526;
            m_constants.Ttriple = 273.16;
            m_constants.Tmin = m_constants.Ttriple + 0.1;
            m_constants.Tmax = 2273.15;
            m_constants.pmin = 1E4;
            m_constants.pmax = 1000E5;
            m_constants.T_critical = 647.096;
            m_constants.p_critical = 22.064E6;
            m_constants.rhomass_critical = 322.0;
            m_constants.molar_mass = 0.018015268;
        }
        cIAPWS_IF97::cIAPWS_IF97(/* args */)
        {
            initialize_data();
        }

        cIAPWS_IF97::~cIAPWS_IF97()
        {
        }

        /**
     * @brief The boundary between regions 2 and 3 (see Fig. 1) is defined by the following simple quadratic pressure-temperature relation, the B23-equation. See equation 5 of \cite IF97.
     *
     * @param P [Pa]
     * @return double
     */
        double cIAPWS_IF97::Boundary_region23_P2T(double P)
        {
            return 0.57254459862746e3 + pow((P/1E6-0.1391883977870e2)/0.10192970039326e-2, 0.5);
        }

        /**
         * @brief Boundary between region 2b and 2c.
         * In order to know whether the T(p,h) equation for subregion 2b or for subregion 2c has to be used for given values of p and h, a special correlation equation for the boundary between subregions 2b and 2c (which approximates s = 5.85 kJ/kg/K #CONST_IF97_S_Region2b2c) is needed; see Fig. 2. This boundary equation, called the B2bc-equation, is a simple quadratic pressure-enthalpy relation which reads Equation 20 of \cite IF97
         *
         * @param P [Pa]
         * @return double [J/kg]
         */
        double cIAPWS_IF97::Boundary_region2b2c_P2H(double P)
        {
            // note that the p_star = 1MPa, h_star = 1 kJ/kg
            return (0.26526571908428e4+pow((P/1E6-0.45257578905948e1)/1.2809002730136e-4, 0.5))*1000;
        }
        /**
         * @brief Boundary between region 2b and 2c. See also #Boundary_region2b2c_P2H.
         *
         * @param H [J/kg]
         * @return double [Pa]
         */
        double cIAPWS_IF97::Boundary_region2b2c_H2P(double H)
        {
            // note that the p_star = 1MPa, h_star = 1 kJ/kg
            double eta = H/1000;
            return (905.84278514723-0.67955786399241*eta+1.2809002730136e-4*eta*eta)*1E6;
        }
        /**
         * @brief Boundary between region3a and 3b. See Equation 1 and Table 2 of \cite IF97-Region3
         * \note This equation is only valid in P range of [#m_constants.p_critical, #CONST_IF97_Pmax_Region1].
         * @param P [Pa]
         * @return double
         */
        double cIAPWS_IF97::Boundary_region3ab_P2H(double P)
        {
            double pi = P/1E6;
            return (0.201464004206875e4 + 3.74696550136983*pi - 0.0219921901054187*pi*pi+0.875131686009950e-4*pi*pi*pi)*1000;
        }
        /**
         * @brief Boundary between region3a and 3b. See Equation 10 and Table 17 of \cite IF97-Region3.
         * \note The valid H range is \f$ [h^{\prime}(#CONST_IF97_Tmin_Region3), h^{\prime\prime}(#CONST_IF97_Tmin_Region3)] \f$. See also pp.17 of \cite IF97-Region3 .
         * @param H [J/kg]
         * @return double [Pa]
         */
        double cIAPWS_IF97::Boundary_region3ab_H2P(double H)
        {
            double eta = H/coeff_boundary3ab.hStar;
            double pi = 0;
            for (int i = 0; i < coeff_boundary3ab.num; i++)
            {
                pi += coeff_boundary3ab.n[i]*pow(eta - 1.02, coeff_boundary3ab.I[i])*pow(eta - 0.608, coeff_boundary3ab.J[i]);
            }
            return pi * coeff_boundary3ab.pStar;
        }
        /**
         * @brief Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::m_constants.pmin, H2O::m_constants.p_critical] See Equation 31 and Table 34 of \cite IF97.
         *
         * @param P [Pa]
         * @return double [K]
         */

        /**
     * @brief Calculate temperature in region 2a by given P and H
     *
     * @param P [Pa]
     * @param H [J/kg]
     * @return double
     */
        double cIAPWS_IF97::Backward_T_PH_region2a(double P, double H)
        {
            // TStar = 1 K;
            double T = 0;
            double pi = P/coeff_backward_T_PH_region2a.pStar;
            double eta = H/coeff_backward_T_PH_region2a.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region2a.num; i++)
            {
                T += coeff_backward_T_PH_region2a.n[i]*pow(pi, coeff_backward_T_PH_region2a.I[i])*pow(eta-2.1, coeff_backward_T_PH_region2a.J[i]);
            }
            return T;
        }
        double cIAPWS_IF97::Backward_T_PH_region2b(double P, double H)
        {
            double pi = P/coeff_backward_T_PH_region2b.pStar;
            double T = 0;
            double eta = H/coeff_backward_T_PH_region2b.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region2b.num; i++)
            {
                T += coeff_backward_T_PH_region2b.n[i]*pow(pi - 2.0, coeff_backward_T_PH_region2b.I[i])*pow(eta - 2.6, coeff_backward_T_PH_region2b.J[i]);
            }
            return T;
        }
        double cIAPWS_IF97::Backward_T_PH_region2c(double P, double H)
        {
            // note that the TStar = 1K
            double pi = P/coeff_backward_T_PH_region2b.pStar;
            double T = 0;
            double eta = H/coeff_backward_T_PH_region2c.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region2c.num; i++)
            {
                T += coeff_backward_T_PH_region2c.n[i]*pow(pi + 25.0, coeff_backward_T_PH_region2c.I[i])*pow(eta - 1.8, coeff_backward_T_PH_region2c.J[i]);
            }
            return T;
        }
        /**
         * @brief The backward equation \f$ T_{3a}(p,h) \f$ for subregion 3a. See equation 2 of \cite IF97-Region3.
         *
         * @param P
         * @param H
         * @return double
         */
        double cIAPWS_IF97::Backward_T_PH_region3a(double P, double H)
        {
            // note that the TStar = 1K
            double pi = P/coeff_backward_T_PH_region3a.pStar;
            double theta = 0;
            double eta = H/coeff_backward_T_PH_region3a.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region3a.num; i++)
            {
                theta += coeff_backward_T_PH_region3a.n[i]*pow(pi + 0.24, coeff_backward_T_PH_region3a.I[i])*pow(eta - 0.615, coeff_backward_T_PH_region3a.J[i]);
            }
            return theta*coeff_backward_T_PH_region3a.TStar;
        }
        /**
         * @brief The backward equation \f$ T_{3b}(p,h) \f$ for subregion 3a. See equation 3 of \cite IF97-Region3.
         *
         * @param P
         * @param H
         * @return double
         */
        double cIAPWS_IF97::Backward_T_PH_region3b(double P, double H)
        {
            // note that the TStar = 1K
            double pi = P/coeff_backward_T_PH_region3b.pStar;
            double theta = 0;
            double eta = H/coeff_backward_T_PH_region3b.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region3b.num; i++)
            {
                theta += coeff_backward_T_PH_region3b.n[i]*pow(pi + 0.298, coeff_backward_T_PH_region3b.I[i])*pow(eta - 0.72, coeff_backward_T_PH_region3b.J[i]);
            }
            return theta*coeff_backward_T_PH_region3b.TStar;
        }
        /**
         * @brief Calculate temperature in region 1 by given P and H. See equation 11 of \cite IF97.
         *
         * @param P [Pa]
         * @param H [J/kg]
         * @return double
         */
        double cIAPWS_IF97::Backward_T_PH_region1(double P, double H)
        {
            // note that the TStar = 1K
            double pi = P/coeff_backward_T_PH_region2b.pStar;
            double T = 0;
            double eta = H/coeff_backward_T_PH_region1.hStar;
            for (int i = 0; i < coeff_backward_T_PH_region1.num; i++)
            {
                T += coeff_backward_T_PH_region1.n[i]*pow(pi, coeff_backward_T_PH_region1.I[i])*pow(eta + 1.0, coeff_backward_T_PH_region1.J[i]);
            }
            return T;
        }

        /**
         * @brief Calculate region index by given pressure and specific enthalpy.
         *
         * @param P [Pa]
         * @param H [J/kg]
         * @return int
         */
        int cIAPWS_IF97::GetRegion_PH(double P, double H)
        {
            if(P>=m_constants.pmin && P<=CONST_IF97_Pmin_Region3) //Region 1, 2, 4, 5
            {
                double T_sat = T_sat_P(P);
                State_Region1 state14, state_Tmin;
                State_Region2 state24, state25;
                State_Region5 state_Tmax;
                getState_Region1(P, m_constants.Tmin, state_Tmin);
                getState_Region1(P, T_sat, state14);
                getState_Region2(P, T_sat, state24);
                getState_Region2(P, CONST_IF97_Tmax_Region2, state25);
                getState_Region5(P, m_constants.Tmax, state_Tmax);
                double h14 = Prop_Region1(state14, Prop_h);
                double h24 = Prop_Region2(state24, Prop_h);
                double h25 = Prop_Region2(state25, Prop_h);
                double hmin = Prop_Region1(state_Tmin, Prop_h);
                double hmax = Prop_Region5(P, m_constants.Tmax, state_Tmax, Prop_h);
                if(H>=hmin && H<=h14)  //Region 1
                {
                    return IF97_REGION1;
                }else if(H>h14 && H<h24) //Region 4: two phase region
                {
                    return IF97_REGION4;
                }else if(H>=h24 && H<=h25) //Region 2
                {
                    if(P<=CONST_IF97_P_Region2a2b) //Region 2a
                    {
                        return IF97_REGION2a;
                    }else if(CONST_IF97_P_Region2a2b<P && P<CONST_IF97_Pmin_Region2c) //Region 2b
                    {
                        return IF97_REGION2b;
                    }else
                    {
                        double h_bd = Boundary_region2b2c_P2H(P);
                        if(H>=h_bd) //Region 2b
                        {
                            return IF97_REGION2b;
                        }else
                        {
                            return IF97_REGION2c; //region 2c
                        }
                    }
                }else if(H>h25 && H<hmax)  //region 5
                {
                    return IF97_REGION5;
                }else
                {
                    return IF97_RegionUndefined;
                    // printf("P = %f Pa, H = %f J/kg\n", P, H);
                    // ERROR("H out of bounds in int cIAPWS_IF97::GetRegion_PH(double P, double H)");
                }
            }else if(P>CONST_IF97_Pmin_Region3 && P<m_constants.p_critical)
            {
                double T_32  = Boundary_region23_P2T(P);
                State_Region1 state_Tmin, state13;
                State_Region2 state32, state25;
                State_Region5 state_Tmax;
                getState_Region1(P, m_constants.Tmin, state_Tmin);
                getState_Region1(P, CONST_IF97_Tmin_Region3, state13);
                getState_Region2(P, T_32, state32);
                getState_Region2(P, CONST_IF97_Tmax_Region2, state25);
                getState_Region5(P, m_constants.Tmax, state_Tmax);
                double hmin = Prop_Region1(state_Tmin, Prop_h);
                double h13  = Prop_Region1(state13, Prop_h);
                double h32  = Prop_Region2( state32, Prop_h);
                double h25  = Prop_Region2( state25, Prop_h);
                double hmax = Prop_Region5(P, m_constants.Tmax, state_Tmax, Prop_h);
                if (H>=hmin && H<=h13) //Region1
                {
                    return IF97_REGION1;
                }else if(H>h13 && H<h32) //Region 3, Region 4
                {
                    if(H<=CONST_IF97_Hl_Tmin_Region3)
                    {
                        return IF97_REGION3a;
                    }else if(H>=CONST_IF97_Hv_Tmin_Region3)
                    {
                        return IF97_REGION3b;
                    }else
                    {
                        double p3ab = Boundary_region3ab_H2P(H); //Use p3ab can only determin in region3 or region 4, but can not determin in 3a, or 3b. We can use H_c to determin region is 3a or 3b.
                        assert(p3ab>=CONST_IF97_Pmin_Region3 && p3ab<m_constants.p_critical);
                        if(P<p3ab)   //Region 4
                        {
                            return IF97_REGION4;
                        }else       //Region 3a or 3b
                        {
                            if(H<=CONST_IF97_H_c)  //Region 3a
                            {
                                return IF97_REGION3a;
                            }else
                            {
                                return IF97_REGION3b;
                            }
                        }
                    }
                }else if(H>=h32 && H<=h25)
                {
                    double h_bd = Boundary_region2b2c_P2H(P);
                    if(H>=h_bd) //Region 2b
                    {
                        return IF97_REGION2b;
                    }else
                    {
                        return IF97_REGION2c; //region 2c
                    }
                }else if(H>h25 && H<hmax) //Region 5
                {
                    return IF97_REGION5;
                }else
                {
                    return IF97_RegionUndefined;
                    // printf("P = %f Pa, H = %f J/kg\n", P, H);
                    // ERROR("H out of bounds in int cIAPWS_IF97::GetRegion_PH(double P, double H)");
                }

            }else if(P>=m_constants.p_critical && P<=CONST_IF97_Pmax_Region5)
            {
                double T_32  = Boundary_region23_P2T(P);
                State_Region1 state_Tmin, state13;
                State_Region2 state32, state25;
                State_Region5 state_Tmax;
                getState_Region1(P, m_constants.Tmin, state_Tmin);
                getState_Region1(P, CONST_IF97_Tmin_Region3, state13);
                getState_Region2(P, T_32, state32);
                getState_Region2(P, CONST_IF97_Tmax_Region2, state25);
                getState_Region5(P, m_constants.Tmax, state_Tmax);
                double hmin = Prop_Region1(state_Tmin, Prop_h);
                double h13  = Prop_Region1(state13, Prop_h);
                double h32  = Prop_Region2(state32, Prop_h);
                double h25  = Prop_Region2( state25, Prop_h);
                double hmax = Prop_Region5(P, m_constants.Tmax, state_Tmax, Prop_h);
                if (H>=hmin && H<=h13)  //Region 1
                {
                    return IF97_REGION1;
                }else if(H>h13 && H<h32) //Region 3
                {
                    double h3ab = Boundary_region3ab_P2H(P);
                    assert(h3ab>h13 && h3ab<h32);
                    if(H<=h3ab)     //Region 3a
                    {
                        return IF97_REGION3a;
                    }else           //Region 3b
                    {
                        return IF97_REGION3b;
                    }
                }else if (H>=h32 && H<=h25) //Region 2
                {
                    double h_bd = Boundary_region2b2c_P2H(P);
                    if(H>=h_bd) //Region 2b
                    {
                        return IF97_REGION2b;
                    }else
                    {
                        return IF97_REGION2c; //region 2c
                    }
                }else if(H>h25 && H<hmax)
                {
                    return IF97_REGION5;
                }else
                {
                    return IF97_RegionUndefined;
                    // printf("P = %f Pa, H = %f J/kg\n", P, H);
                    // ERROR("H out of bounds in int cIAPWS_IF97::GetRegion_PH(double P, double H)");
                }
            }else if(P>CONST_IF97_Pmax_Region5 && P<=CONST_IF97_Pmax_Region1)
            {
                double T_32  = Boundary_region23_P2T(P);
                State_Region1 state_Tmin, state13;
                State_Region2 state32, state25;
                State_Region5 state_Tmax;
                getState_Region1(P, m_constants.Tmin, state_Tmin);
                getState_Region1(P, CONST_IF97_Tmin_Region3, state13);
                getState_Region2(P, T_32, state32);
                getState_Region2(P, CONST_IF97_Tmax_Region2, state25);
                getState_Region5(P, m_constants.Tmax, state_Tmax);
                double hmin = Prop_Region1(state_Tmin, Prop_h);
                double h13  = Prop_Region1(state13, Prop_h);
                double h32  = Prop_Region2(state32, Prop_h);
                double h25  = Prop_Region2( state25, Prop_h);
                double hmax = Prop_Region5(P, m_constants.Tmax, state_Tmax, Prop_h);
                if (H>=hmin && H<=h13)  //Region 1
                {
                    return IF97_REGION1;
                }else if(H>h13 && H<h32) //Region 3
                {
                    double h3ab = Boundary_region3ab_P2H(P);
                    assert(h3ab>h13 && h3ab<h32);
                    if(H<=h3ab)     //Region 3a
                    {
                        return IF97_REGION3a;
                    }else           //Region 3b
                    {
                        return IF97_REGION3b;
                    }
                }else if (H>=h32 && H<=h25) //Region 2
                {
                    double h_bd = Boundary_region2b2c_P2H(P);
                    if(H>=h_bd) //Region 2b
                    {
                        return IF97_REGION2b;
                    }else
                    {
                        return IF97_REGION2c; //region 2c
                    }
                }else if(H>h25 && H<hmax)    //Undefined region
                {
                    return IF97_RegionUndefined;
                }else
                {
                    return IF97_RegionUndefined;
                    // printf("P = %f Pa, H = %f J/kg\n", P, H);
                    // ERROR("H out of bounds in int cIAPWS_IF97::GetRegion_PH(double P, double H)");
                }
            }
            else
            {
                return IF97_RegionUndefined;
                // ERROR("Fatal error in int cIAPWS_IF97::GetRegion_PH(double P, double H): input pressure out of bound.\nP = "+std::to_string(P)+", pressure bound: ["+std::to_string(m_constants.pmin)+", "+std::to_string(CONST_IF97_Pmax_Region1)+"] Pa");
            }

            return IF97_RegionUndefined;
        }
        /**
         * @brief Calculate region index by given pressure and temperature.
         * \note In P-T space, there is no need to introduce subregions, e.g. 2a, 2b, 2c and 3a, 3b.
         * @param P [Pa]
         * @param T [K]
         * @return int
         */
        int cIAPWS_IF97::GetRegion_PT(double P, double T)
        {
            if(P>=m_constants.pmin && P<=CONST_IF97_Pmin_Region3) //Region 1, 2, 5, maybe 4
            {
                double T_sat = T_sat_P(P);
                if (T<=T_sat && T>=m_constants.Tmin)        //Region 1
                {
                    return IF97_REGION1;
                }else if(T>T_sat && T<CONST_IF97_Tmax_Region2)        //Region 2
                {
                    return IF97_REGION2;
                }else if(T>CONST_IF97_Tmax_Region2 && T<=CONST_IF97_Tmax_Region5)
                {
                    return IF97_REGION5;
                }
                else
                {
                    printf("P = %f Pa, T = %f K\n", P, T);
                    ERROR("T out of bound in int cIAPWS_IF97::GetRegion_PT(double P, double T)");
                }
            }else if(P>CONST_IF97_Pmin_Region3 && P<=CONST_IF97_Pmax_Region5)
            {
                double T_b23 = Boundary_region23_P2T(P);
                if(T>=m_constants.Tmin && T<=CONST_IF97_Tmin_Region3)
                {
                    return IF97_REGION1;
                }else if(T>CONST_IF97_Tmin_Region3 && T<T_b23)
                {
                    return IF97_REGION3;
                }else if(T>=T_b23 && T<CONST_IF97_Tmax_Region2)
                {
                    return IF97_REGION2;
                }else if(T>=CONST_IF97_Tmax_Region2 && T<=CONST_IF97_Tmax_Region5)
                {
                    return IF97_REGION5;
                }else
                {
                    printf("P = %f Pa, T = %f K\n", P, T);
                    ERROR("T out of bound in int cIAPWS_IF97::GetRegion_PT(double P, double T)");
                }
            }else if(P>CONST_IF97_Pmax_Region5 && P<=CONST_IF97_Pmax_Region1)
            {
                double T_b23 = Boundary_region23_P2T(P);
                if(T>=m_constants.Tmin && T<=CONST_IF97_Tmin_Region3)
                {
                    return IF97_REGION1;
                }else if(T>CONST_IF97_Tmin_Region3 && T<T_b23)
                {
                    return IF97_REGION3;
                }else if(T>=T_b23 && T<CONST_IF97_Tmax_Region2)
                {
                    return IF97_REGION2;
                }else if(T>=CONST_IF97_Tmax_Region2 && T<=CONST_IF97_Tmax_Region5)
                {
                    return IF97_RegionUndefined;
                }else
                {
                    printf("P = %f Pa, T = %f K\n", P, T);
                    ERROR("T out of bound in int cIAPWS_IF97::GetRegion_PT(double P, double T)");
                }
            }
            else
            {
                ERROR("Fatal error in int cIAPWS_IF97::GetRegion_PT(double P, double T): input pressure out of bound.\nP = "+std::to_string(P)+", pressure bound: ["+std::to_string(m_constants.pmin)+", "+std::to_string(CONST_IF97_Pmax_Region1)+"] Pa");
            }
            return IF97_RegionUndefined;
        }
        /**
         * @brief Calculate Gibbs free energy and its derivatives by given pressure and temperature.
         *
         * @param P [Pa]
         * @param T [K]
         * @param state
         */
        void cIAPWS_IF97::getState_Region1(double P, double T, State_Region1& state)
        {
            state.T = T; state.P = P;
            state.tau = coeff_region1.TStar/T;
            state.pi  = P/coeff_region1.pStar;
            state.RT = m_constants.R*T;
            state.gamma = {0, 0, 0, 0, 0, 0};
            for (int i = 0; i < coeff_region1.num; i++)
            {
                state.gamma.value += coeff_region1.n[i]*pow(7.1 - state.pi, coeff_region1.I[i]) * pow(state.tau - 1.222, coeff_region1.J[i]);
                state.gamma.p -= coeff_region1.n[i]*coeff_region1.I[i]*pow(7.1 - state.pi, coeff_region1.I[i]-1.0)*pow(state.tau-1.222, coeff_region1.J[i]);
                state.gamma.pp += coeff_region1.n[i]*coeff_region1.I[i]*(coeff_region1.I[i]-1.0)*pow(7.1 - state.pi, coeff_region1.I[i]-2.0)*pow(state.tau - 1.222, coeff_region1.J[i]);
                state.gamma.t += coeff_region1.n[i]*pow(7.1 - state.pi, coeff_region1.I[i])*coeff_region1.J[i]*pow(state.tau - 1.222, coeff_region1.J[i]-1.0);
                state.gamma.tt += coeff_region1.n[i]*pow(7.1 - state.pi, coeff_region1.I[i])*coeff_region1.J[i]*(coeff_region1.J[i]-1.0)*pow(state.tau - 1.222, coeff_region1.J[i]-2.0);
                state.gamma.pt -= coeff_region1.n[i]*coeff_region1.I[i]*pow(7.1 - state.pi, coeff_region1.I[i]-1.0)*coeff_region1.J[i]*pow(state.tau - 1.222, coeff_region1.J[i]-1.0);
            }
        }
        /**
         * @brief Calculate dimensionless Gibbs free energy in the region 2. ideal-gas part \f$ \gamma^o \f$ and the residual part \f$ \gamma^r \f$. See equation 15-17 of \cite IF97.
         *
         * @param P
         * @param T
         * @param gammao
         * @param gammar
         */
        void cIAPWS_IF97::getState_Region2(double P, double T, State_Region2& state)
        {
            state.P = P; state.T = T;
            state.pi = P/coeff_region2.pStar;
            state.tau = coeff_region2.TStar/T;
            state.RT = m_constants.R*T;
            // ideal-gass part: gammo
            state.gammao = {log(state.pi), 1.0/state.pi, -1.0/(state.pi*state.pi), 0, 0, 0};
            for (int i = 0; i < coeff_region2.numo; i++)
            {
                state.gammao.value += coeff_region2.no[i]*pow(state.tau, coeff_region2.Jo[i]);
                state.gammao.t += coeff_region2.no[i]*coeff_region2.Jo[i]*pow(state.tau, coeff_region2.Jo[i]-1.0);
                state.gammao.tt += coeff_region2.no[i]*coeff_region2.Jo[i]*(coeff_region2.Jo[i] - 1.0)*pow(state.tau, coeff_region2.Jo[i]-2.0);
            }
            // residual part
            state.gammar = {0, 0, 0, 0, 0, 0};
            for (int i = 0; i < coeff_region2.numr; i++)
            {
                state.gammar.value += coeff_region2.nr[i]*pow(state.pi, coeff_region2.Ir[i])*pow(state.tau - 0.5, coeff_region2.Jr[i]);
                state.gammar.p += coeff_region2.nr[i]*coeff_region2.Ir[i]*pow(state.tau-0.5, coeff_region2.Jr[i]);
                state.gammar.pp += coeff_region2.nr[i]*coeff_region2.Ir[i]*(coeff_region2.Ir[i]-1.0)*pow(state.pi, coeff_region2.Ir[i]-2.0)*pow(state.tau-0.5, coeff_region2.Jr[i]);
                state.gammar.t += coeff_region2.nr[i]*pow(state.pi, coeff_region2.Ir[i])*coeff_region2.Jr[i]*pow(state.tau - 0.5, coeff_region2.Jr[i]-1.0);
                state.gammar.tt += coeff_region2.nr[i]*pow(state.pi, coeff_region2.Ir[i])*coeff_region2.Jr[i]*(coeff_region2.Jr[i]-1.0)*pow(state.tau-0.5, coeff_region2.Jr[i]-2.0);
                state.gammar.pt += coeff_region2.nr[i]*coeff_region2.Ir[i]*pow(state.pi, coeff_region2.Ir[i]-1.0)*coeff_region2.Jr[i]*pow(state.tau-0.5, coeff_region2.Jr[i]-1.0);
            }
            // calculate gamma
            state.gamma = state.gammao + state.gammar;
        }
        /**
         * @brief Calculate dimensionless Gibbs free energy in the region 5. ideal-gas part \f$ \gamma^o \f$ and the residual part \f$ \gamma^r \f$. See equation 32-34 of \cite IF97.
         *
         * @param P
         * @param T
         * @param state
         */
        void cIAPWS_IF97::getState_Region5(double P, double T, State_Region5& state)
        {
            state.T = T; state.P = P;
            state.tau = coeff_region1.TStar/T;
            state.pi  = P/coeff_region1.pStar;
            state.RT = m_constants.R*T;
            // ideal-gass part: gammo
            state.gammao = {log(state.pi), 1.0/state.pi, -1.0/(state.pi*state.pi), 0, 0, 0};
            for (int i = 0; i < coeff_region5.numo; i++)
            {
                state.gammao.value += coeff_region5.no[i]*pow(state.tau, coeff_region5.Jo[i]);
                state.gammao.t += coeff_region5.no[i]*coeff_region5.Jo[i]*pow(state.tau, coeff_region5.Jo[i]-1.0);
                state.gammao.tt += coeff_region5.no[i]*coeff_region5.Jo[i]*(coeff_region5.Jo[i]-1.0)*pow(state.tau, coeff_region5.Jo[i]-2.0);
                // state.gammao.pt = 0;
            }
            // residual part
            state.gammar = {0, 0, 0, 0, 0, 0};
            for (int i = 0; i < coeff_region5.numr; i++)
            {
                state.gammar.value += coeff_region5.nr[i]*pow(state.pi, coeff_region5.Ir[i])*pow(state.tau, coeff_region5.Jr[i]);
                state.gammar.p += coeff_region5.nr[i]*coeff_region5.Ir[i]*pow(state.pi, coeff_region5.Ir[i]-1.0)*pow(state.tau, coeff_region5.Jr[i]);
                state.gammar.pp += coeff_region5.nr[i]*coeff_region5.Ir[i]*(coeff_region5.Ir[i]-1.0)*pow(state.pi, coeff_region5.Ir[i]-2.0)*pow(state.tau, coeff_region5.Jr[i]);
                state.gammar.t += coeff_region5.nr[i]*pow(state.pi, coeff_region5.Ir[i])*coeff_region5.Jr[i]*pow(state.tau, coeff_region5.Jr[i]-1.0);
                state.gammar.tt += coeff_region5.nr[i]*pow(state.pi, coeff_region5.Ir[i])*coeff_region5.Jr[i]*(coeff_region5.Jr[i]-1.0)*pow(state.pi, coeff_region5.Jr[i]-2.0);
                state.gammar.pt += coeff_region5.nr[i]*coeff_region5.Ir[i]*pow(state.pi, coeff_region5.Ir[i]-1.0)*coeff_region5.Jr[i]*pow(state.tau, coeff_region5.Jr[i]-1.0);
            }
            // calculate gamma
            state.gamma = state.gammao + state.gammar;
        }
        /**
         * @brief Calculate thermodynamic properties by given dimensionless Gibbs free energy and its derivatives for the region 1.
         *
         * @param gamma
         * @param which
         * @return double
         */
        double cIAPWS_IF97::Prop_Region1(const State_Region1& state, BasicThermodynamicProperties which)
        {
            switch (which)
            {
                case Prop_rho:
                    return state.P/(state.RT*state.pi*state.gamma.p);
                case Prop_v:
                    return state.RT*state.pi*state.gamma.p/state.P;
                    break;
                case Prop_u:
                    return state.RT*(state.tau*state.gamma.t - state.pi*state.gamma.p);
                case Prop_s:
                    return m_constants.R*(state.tau*state.gamma.t - state.gamma.value);
                case Prop_h:
                    return state.RT*state.tau*state.gamma.t;
                case Prop_cp:
                    return -m_constants.R*state.tau*state.tau*state.gamma.tt;
                case Prop_cv:
                    return m_constants.R*(-state.tau*state.tau*state.gamma.tt + pow(state.gamma.p - state.tau*state.gamma.pt, 2.0)/state.gamma.pp);
                case Prop_w:
                    return state.RT*state.gamma.p*state.gamma.p/(pow(state.gamma.p - state.tau*state.gamma.pt, 2.0)/(state.tau*state.tau*state.gamma.tt) - state.gamma.pp);
                default:
                ERROR("Unsupported basic property in region 1: "+std::to_string(which));
                    break;
            }
            return 0;
        }
        double cIAPWS_IF97::Prop_Region2(State_Region2 state, BasicThermodynamicProperties which)
        {
            switch (which)
            {
                case Prop_rho:
                    return state.P/(state.RT*state.pi*state.gamma.p);
                case Prop_v:
                    return state.RT*state.pi*state.gamma.p/state.P;
                    break;
                case Prop_u:
                    return state.RT*(state.tau*state.gamma.t - state.pi*state.gamma.p);
                case Prop_s:
                    return m_constants.R*(state.tau*state.gamma.t - state.gamma.value);
                case Prop_h:
                    return state.RT*state.tau*state.gamma.t;
                case Prop_cp:
                    return -m_constants.R*state.tau*state.tau*state.gamma.tt;
                case Prop_cv:
                    return m_constants.R*(-state.tau*state.tau * state.gamma.tt) - pow(1+state.pi*state.gammar.p - state.tau*state.pi*state.gammar.pt, 2.0)/(1.0 - state.pi*state.pi*state.gammar.pp);
                case Prop_w:
                    return state.RT * (1.0 + 2.0*state.pi*state.gammar.p + state.pi*state.pi*pow(state.gammar.p, 2.0))/((1.0 - state.pi*state.pi*state.gammar.pp) + pow(1.0+state.pi*state.gammar.p - state.tau*state.pi*state.gammar.pt, 2.0)/(state.tau*state.tau*(state.gammao.tt + state.gammar.tt)));
                default:
                ERROR("Unsupported basic property in region 2: "+std::to_string(which));
                    break;
            }
            return 0;
        }
        double cIAPWS_IF97::Prop_Region5(double P, double T, State_Region5 state, BasicThermodynamicProperties which)
        {
            switch (which)
            {
                case Prop_rho:
                    return state.P/(state.RT*state.pi*state.gamma.p);
                case Prop_v:
                    return state.RT*state.pi*state.gamma.p/state.P;
                    break;
                case Prop_u:
                    return state.RT*(state.tau*state.gamma.t - state.pi*state.gamma.p);
                case Prop_s:
                    return m_constants.R*(state.tau*state.gamma.t - state.gamma.value);
                case Prop_h:
                    return state.RT*state.tau*state.gamma.t;
                case Prop_cp:
                    return -m_constants.R*state.tau*state.tau*state.gamma.tt;
                case Prop_cv:
                    return m_constants.R*(-state.tau*state.tau * state.gamma.tt) - pow(1+state.pi*state.gammar.p - state.tau*state.pi*state.gammar.pt, 2.0)/(1.0 - state.pi*state.pi*state.gammar.pp);
                case Prop_w:
                    return state.RT * (1.0 + 2.0*state.pi*state.gammar.p + state.pi*state.pi*pow(state.gammar.p, 2.0))/((1.0 - state.pi*state.pi*state.gammar.pp) + pow(1.0+state.pi*state.gammar.p - state.tau*state.pi*state.gammar.pt, 2.0)/(state.tau*state.tau*(state.gammao.tt + state.gammar.tt)));
                default:
                ERROR("Unsupported basic property in region 2: "+std::to_string(which));
                    break;
            }
            return 0;
        }
        /**
     * @brief Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::P_MIN, H2O::P_c] See Equation 31 and Table 34 of \cite IF97.
     *
     * @param P [Pa]
     * @return double [K]
     */
        double cIAPWS_IF97::T_sat_P(const double& P)
        {
            //pStar = 1 MPa
            double beta = pow(P/1E6, 0.25);
            double beta_sqr = beta*beta;
            double E = beta_sqr + table34.n[2]*beta + table34.n[5];
            double F = table34.n[0] * beta_sqr + table34.n[3]*beta + table34.n[6];
            double G = table34.n[1] * beta_sqr + table34.n[4]*beta + table34.n[7];
            double D = 2*G/(-F-pow(F*F-4*E*G, 0.5));
            return (table34.n[9]+D-sqrt(pow(table34.n[9]+D, 2.0) - 4*(table34.n[8]+table34.n[9]*D)))/2.0;
        }
        /**
         * @brief Calculate saturated liquid enthalpy and saturated vapor enthalpy for a given P.
         * \note This function is only valid in P range [#P_MIN, #CONST_IF97_Pmin_Region3].
         * @param P
         * @param H_l
         * @param H_v
         */
        void cIAPWS_IF97::H_sat_P(const double& P, double& H_l, double& H_v)
        {
            if(P>=m_constants.pmin && P<=CONST_IF97_Pmin_Region3)
            {
                double T_sat = T_sat_P(P);
                State_Region1 state14;
                getState_Region1(P, T_sat, state14);
                H_l = Prop_Region1(state14, Prop_h);
                State_Region2 state24;
                getState_Region2(P, T_sat, state24);
                H_v = Prop_Region2(state24, Prop_h);
            }else
            {
                H_l = NAN; H_v = NAN;
            }

        }
        /**
         * @brief Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::T_MIN, H2O::T_c] See Equation 30 and Table 34 of \cite IF97.
         *
         * @param T [K]
         * @return double [Pa]
         */
        double cIAPWS_IF97::P_sat_T(const double& T)
        {
            //pStar = 1MPa, TStar = 1K
            double theta = T + table34.n[8]/(T - table34.n[9]);
            double theta_sqr = theta*theta;
            double A = theta_sqr + table34.n[0]*theta + table34.n[1];
            double B = table34.n[2]*theta_sqr + table34.n[3]*theta + table34.n[4];
            double C = table34.n[5]*theta_sqr + table34.n[6]*theta + table34.n[7];
            return pow(2*C/(-B + sqrt(B*B - 4*A*C)), 4.0) * 1E6;
        }
        /**
     * @brief Calculate temperature T by given pressure and specific enthalpy.
     *
     * @param P
     * @param H
     * @return double
     */
        double cIAPWS_IF97::T_PH(const double& P, const double& H)
        {
            int region = GetRegion_PH(P, H);
            switch (region)
            {
                case IF97_REGION1:
                    return Backward_T_PH_region1(P, H);
                case IF97_REGION2a:
                    return Backward_T_PH_region2a(P, H);
                case IF97_REGION2b:
                    return Backward_T_PH_region2b(P, H);
                case IF97_REGION2c:
                    return Backward_T_PH_region2c(P, H);
                case IF97_REGION3a:
                    return Backward_T_PH_region3a(P, H);
                case IF97_REGION3b:
                    return Backward_T_PH_region3b(P, H);
                case IF97_REGION4:
                    return T_sat_P(P);
                case IF97_RegionUndefined:
                    return NAN;
                default:
                ERROR("Region "+std::to_string(region)+" is not supported for backward temperature");
                    break;
            }
        }

        void cIAPWS_IF97::GetRegion_PT(const std::vector<double> P, const std::vector<double> T, std::vector<int>& res)
        {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = GetRegion_PT(P[i], T[i]);
        }
        void cIAPWS_IF97::GetRegion_PH(const std::vector<double> P, const std::vector<double> H, std::vector<int>& res)
        {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = GetRegion_PH(P[i], H[i]);
        }
        void cIAPWS_IF97::P_sat_T(const std::vector<double> T, std::vector<double>& res)
        {
            res.clear();
            res.resize(T.size());
            for (size_t i = 0; i < T.size(); i++)res[i] = P_sat_T(T[i]);
        }
        void cIAPWS_IF97::T_sat_P(const std::vector<double> P, std::vector<double>& res)
        {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = T_sat_P(P[i]);
        }
        void cIAPWS_IF97::H_sat_P(const std::vector<double> P, std::vector<double>& H_l,std::vector<double>& H_v)
        {
            H_l.clear(); H_v.clear();
            H_l.resize(P.size());
            H_v.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)H_sat_P(P[i], H_l[i], H_v[i]);
        }
        void cIAPWS_IF97::T_PH(const std::vector<double> P, const std::vector<double> H, std::vector<double>& res)
        {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = T_PH(P[i], H[i]);
        }
        void cIAPWS_IF97::Boundary_region3ab_H2P(const std::vector<double> H, std::vector<double>& res)
        {
            res.clear();
            res.resize(H.size());
            for (size_t i = 0; i < H.size(); i++)res[i] = Boundary_region3ab_H2P(H[i]);
        }
        void cIAPWS_IF97::Boundary_region3ab_P2H(const std::vector<double> P, std::vector<double>& res)
        {
            res.clear();
            res.resize(P.size());
            for (size_t i = 0; i < P.size(); i++)res[i] = Boundary_region3ab_P2H(P[i]);
        }

    };

};