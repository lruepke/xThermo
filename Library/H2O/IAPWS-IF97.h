/**
 * @file IAPWS-IF97.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of some parts of IAPWS-IF97 EOS.
 * @version 0.1
 * @date 2022-04-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef IAPWS_IF97_xThermal_H
#define IAPWS_IF97_xThermal_H

#include "thermo.h"

namespace xThermal
{
    #define Name_Backend_IAPWS_IF97 "IAPWS-IF97"
    namespace IAPWS_IF97
    {
        /**
         * @defgroup REGIONS_IF97 Define IAPWS-IF97 region index.
         *
         * @{
         */
        #define IF97_REGION1    1
        #define IF97_REGION2    2
        #define IF97_REGION3    3
        #define IF97_REGION4    4
        #define IF97_REGION5    5
        /**  Used for backward equation */
        #define IF97_REGION2a   6
        /**  Used for backward equation */
        #define IF97_REGION2b   7
        /**  Used for backward equation */
        #define IF97_REGION2c   8
        /**  Used for backward equation */
        #define IF97_REGION3a   9
        /**  Used for backward equation */
        #define IF97_REGION3b   10
        #define IF97_RegionUndefined -1
        /** @} */

        /**
         * @defgroup CONSTS_IF97  Constants of IAPWS-IF97.
         *
         * @{
         */

        /**  The boundary between the subregions 2a and 2b is the isobar \f$ p = 4 MPa \f$, see \cite IF97 */
        #define CONST_IF97_P_Region2a2b  4E6
                /**  Equations (20) and (21) give the boundary line between subregions 2b and 2c from the saturation state at T=554.485 K and \f$ p_s \f$ =6.54670 MPa to T=1019.32 K and p=100 MPa. See \cite IF97 */
        #define CONST_IF97_Pmin_Region2c 6.54670E6
                /**  the boundary between the subregions 2b and 2c corresponds to the entropy line s = 5.85 kJ/kg/K. See section 6.3 of \cite IF97 */
        #define CONST_IF97_S_Region2b2c 5.85E3
                /**  Minimum pressure of the region 3. See pp.6 of \cite IF97 */
        #define CONST_IF97_Pmin_Region3 16.5292E6
                /**  Maximum temperature of the region 3. See pp.6 of \cite IF97. */
        #define CONST_IF97_Tmin_Region3 623.15
                /**  Saturated liquid enthalpy at #CONST_IF97_Tmin_Region3 , see section 4.3, pp.17 of \cite IF97-Region3. */
        #define CONST_IF97_Hl_Tmin_Region3 1.670858218E6
                /**  Saturated vapor enthalpy at #CONST_IF97_Tmin_Region3 , see section 4.3, pp.17 of \cite IF97-Region3. */
        #define CONST_IF97_Hv_Tmin_Region3 2.563592004E6
                /**  Maximum pressure [MPa] of the region 1, 2, 3, 4 */
        #define CONST_IF97_Pmax_Region1 100E6
                /**  Maximum temperature [K] of the retion 2 */
        #define CONST_IF97_Tmax_Region2 1073.15
                /**  Maximum pressure [MPa] of the region 5 */
        #define CONST_IF97_Pmax_Region5 50E6
                /**  Maximum temperature [K] of the retion 5 */
        #define CONST_IF97_Tmax_Region5 2273.15
                /**  Specific enthalpy of critical point, calculate from \link IAPWS_IF97::cIAPWS_IF97::Boundary_region3ab_P2H \endlink  (H2O::P_c)  */
        #define CONST_IF97_H_c   2.08754684511650027707E+06
        /** @} */

        /**
         * @brief Dimensionless Gibbs free energy \f$ \gamma = g/(RT) \f$ and its partial derivatives.
         *
         */
        struct GibbsEnergy_dimensionless
        {
            double value;    /**< \f$ \gamma\f$ */
            double p;        /**< \f$ \left[ \frac{\partial \gamma}{\partial \pi} \right]_{\tau} \f$ */
            double pp;       /**< \f$ \left[ \frac{\partial^2 \gamma}{\partial \pi^2} \right]_{\tau} \f$ */
            double t;        /**< \f$ \left[ \frac{\partial \gamma}{\partial \tau} \right]_{\pi} \f$ */
            double tt;       /**< \f$ \left[ \frac{\partial^2 \gamma}{\partial \tau^2} \right]_{\pi} \f$ */
            double pt;       /**< \f$ \frac{\partial^2 \gamma}{\partial \pi \partial \tau} \f$ */
            friend std::ostream& operator<<(std::ostream& os, const GibbsEnergy_dimensionless& gamma)
            {
                return os << "v: "<<gamma.value<<", p: "<<gamma.p<<", pp: "<<gamma.pp<<", t: "<<gamma.t<<", tt: "<<gamma.tt<<", pt: "<<gamma.pt;
            }
            bool operator== (const GibbsEnergy_dimensionless &gamma) const
            {
                return (fabs(value-gamma.value) + fabs(p-gamma.p) + fabs(pp-gamma.pp) + fabs(t-gamma.t) + fabs(tt-gamma.tt) + fabs(pt-gamma.pt))/6.0 < 1E-7;
            }
            GibbsEnergy_dimensionless operator+ (const GibbsEnergy_dimensionless &gamma) const
            {
                GibbsEnergy_dimensionless tmp = {value + gamma.value, p + gamma.p, pp + gamma.pp, t+gamma.t, tt+gamma.tt, pt+gamma.pt};
                return tmp;
            }
            GibbsEnergy_dimensionless operator- (const GibbsEnergy_dimensionless &gamma) const
            {
                GibbsEnergy_dimensionless tmp = {value - gamma.value, p - gamma.p, pp-gamma.pp, t-gamma.t, tt-gamma.tt, pt-gamma.pt};
                return tmp;
            }
        };

        enum BasicThermodynamicProperties
        {
            Prop_v, Prop_u, Prop_s, Prop_h, Prop_cp, Prop_cv, Prop_w, Prop_rho
        };

        struct State_Region1
        {
            double T, P, pi, tau, RT;
            GibbsEnergy_dimensionless gamma;
        };
        struct State_Region2
        {
            double T, P, pi, tau, RT;
            GibbsEnergy_dimensionless gammao, gammar, gamma;
        };
        struct State_Region5
        {
            double T, P, pi, tau, RT;
            GibbsEnergy_dimensionless gammao, gammar, gamma;
        };

        class cIAPWS_IF97 : public cxThermal
        {
        private:
            /**
         * @brief Numerical values of the coefficients and exponents of the dimensionless Gibbs free energy for region 1, Eq. (7) of \cite IF97
         *
         */
            struct Coeff_Region1
            {
                const int num = 34;
                double I[34] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                                2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
                double J[34] = {-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0,
                                6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};
                double n[34] = {0.14632971213167, -0.84548187169114, -0.37563603672040e1,
                                0.33855169168385e1, -0.95791963387872, 0.15772038513228,
                                -0.16616417199501e-1, 0.81214629983568e-3, 0.28319080123804e-3,
                                -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,
                                -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,
                                -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5,
                                -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5,
                                -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6,
                                -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,
                                -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,
                                0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23,
                                -0.93537087292458e-25};
                double pStar = 16.53E6;     /**< \f$ p^*  \f$ */
                double TStar = 1386;        /**< \f$ T^*  \f$ */
            }coeff_region1;
            /**
             * @brief Numerical values of the coefficients and exponents of the ideal-gas part \f$ \gamma^o \f$ o of the dimensionless Gibbs free energy for region 2, Eq. (16)a, and the residual part \f$ \gamma^r \f$
             *
             */
            struct Coeff_Region2
            {
                const int numr = 43;
                const int numo = 9;
                // coeffs for ideal-gass part
                double Jo[9] = {0, 1, -5, -4, -3, -2, -1, 2, 3};
                double no[9] = {-0.96927686500217E+01, 0.10086655968018E+02, -0.56087911283020E-02,
                                0.71452738081455E-01, -0.40710498223928E+00, 0.14240819171444E+01,
                                -0.43839511319450E+01, -0.28408632460772E+00, 0.21268463753307E-01};
                // coeffs for residual part
                double Ir[43] = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7,
                                 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24,
                                 24, 24};
                double Jr[43] = {0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35,
                                 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39,
                                 26, 40, 58};
                double nr[43] = {-0.0017731742473212999, -0.017834862292357999, -0.045996013696365003,
                                 -0.057581259083432, -0.050325278727930002, -3.3032641670203e-05,
                                 -0.00018948987516315, -0.0039392777243355001, -0.043797295650572998,
                                 -2.6674547914087001e-05, 2.0481737692308999e-08,
                                 4.3870667284435001e-07, -3.2277677238570002e-05, -0.0015033924542148,
                                 -0.040668253562648998, -7.8847309559367001e-10,
                                 1.2790717852285001e-08, 4.8225372718507002e-07,
                                 2.2922076337661001e-06, -1.6714766451061001e-11,
                                 -0.0021171472321354998, -23.895741934103999, -5.9059564324270004e-18,
                                 -1.2621808899101e-06, -0.038946842435739003, 1.1256211360459e-11,
                                 -8.2311340897998004, 1.9809712802088e-08, 1.0406965210174e-19,
                                 -1.0234747095929e-13, -1.0018179379511e-09, -8.0882908646984998e-11,
                                 0.10693031879409, -0.33662250574170999, 8.9185845355420999e-25,
                                 3.0629316876231997e-13, -4.2002467698208001e-06,
                                 -5.9056029685639003e-26, 3.7826947613457002e-06,
                                 -1.2768608934681e-15, 7.3087610595061e-29, 5.5414715350778001e-17,
                                 -9.4369707241209998e-07};
                double TStar = 540; /**< \f$ T^*  \f$ */
                double pStar = 1E6; /**< \f$ p^*  \f$ */
            }coeff_region2;
            /**
             * @brief Numerical values of the coefficients and exponents of the residual part  r of the dimensionless Gibbs free energy for region 5. See Table 37, 38 in \cite IF97.
             *
             */
            struct Coeff_Region5
            {
                const int numr = 6, numo = 6;
                // ideal-gas part
                double Jo[6] = {0, 1, -3, -2, -1, 2};
                double no[6] = {-0.13179983674201e2, 0.68540841634434e1, -0.24805148933466e-1,
                                0.36901534980333, -0.31161318213925e1, -0.32961626538917};
                // residual part
                double Ir[6] = {1, 1, 1, 2, 2, 3};
                double Jr[6] = {1, 2, 3, 3, 9, 7};
                double nr[6] = {0.15736404855259e-2, 0.90153761673944e-3, -0.50270077677648e-2,
                                0.22440037409485e-5, -0.41163275453471e-5, 0.37919454822955e-7};
                double TStar = 1000;
                double pStar = 1E6;
            }coeff_region5;
            struct Coeff_Boundary3ab
            {
                double num=14;
                double I[14] = {0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36};
                double J[14] = {0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24};
                double n[14] = {0.600073641753024, -0.936203654849857e1, 0.246590798594147e2,
                                -0.107014222858224e3, -0.915821315805768e14, -0.862332011700662e4,
                                -0.235837344740032e2, 0.252304969384128e18, -0.389718771997719e19,
                                -0.333775713645296e23, 0.356499469636328e11, -0.148547544720641e27,
                                0.330611514838798e19, 0.813641294467829e38};
                double pStar = 22E6;    /**< \f$ p^*  \f$ */
                double hStar = 2600E3; /**< \f$ h^*  \f$ */
            }coeff_boundary3ab;

            /**
             * @brief Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2a, Eq. (22). See Table 20 of \cite IF97.
             *
             */
            struct Coeff_Backward_T_PH_Region2a
            {
                static const int num=34;
                double I[num] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                                 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7};
                double J[num] = {0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40,
                                 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28};
                double n[num] = {0.10898952318288e4, 0.84951654495535e3, -0.10781748091826e3,
                                 0.33153654801263e2, -0.74232016790248e1, 0.11765048724356e2,
                                 0.18445749355790e1, -0.41792700549624e1, 0.62478196935812e1,
                                 -0.17344563108114e2, -0.20058176862096e3, 0.27196065473796e3,
                                 -0.45511318285818e3, 0.30919688604755e4, 0.25226640357872e6,
                                 -0.61707422868339e-2, -0.31078046629583, 0.11670873077107e2,
                                 0.12812798404046e9, -0.98554909623276e9, 0.28224546973002e10,
                                 -0.35948971410703e10, 0.17227349913197e10, -0.13551334240775e5,
                                 0.12848734664650e8, 0.13865724283226e1, 0.23598832556514e6,
                                 -0.13105236545054e8, 0.73999835474766e4, -0.55196697030060e6,
                                 0.37154085996233e7, 0.19127729239660e5, -0.41535164835634e6,
                                 -0.62459855192507e2};
                double pStar = 1E6;     /**< \f$ p^*  \f$ */
                double hStar = 2000E3;  /**< \f$ h^{*} \f$ */
                // double TStar = 1;    /**< \f$ T^* \f$ */
            }coeff_backward_T_PH_region2a;

            /**
             * @brief Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2b, Eq. (23). See Table 21 of \cite IF97.
             *
             */
            struct Coeff_Backward_T_PH_Region2b
            {
                static const int num=38;
                double I[num] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3,
                                 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9};
                double J[num] = {0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18,
                                 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40};
                double n[num] = {0.14895041079516e4, 0.74307798314034e3, -0.97708318797837e2,
                                 0.24742464705674e1, -0.63281320016026, 0.11385952129658e1,
                                 -0.47811863648625, 0.85208123431544e-2, 0.93747147377932,
                                 0.33593118604916e1, 0.33809355601454e1, 0.16844539671904,
                                 0.73875745236695, -0.47128737436186, 0.15020273139707,
                                 -0.21764114219750e-2, -0.21810755324761e-1, -0.10829784403677,
                                 -0.46333324635812e-1, 0.71280351959551e-4, 0.11032831789999e-3,
                                 0.18955248387902e-3, 0.30891541160537e-2, 0.13555504554949e-2,
                                 0.28640237477456e-6, -0.10779857357512e-4, -0.76462712454814e-4,
                                 0.14052392818316e-4, -0.31083814331434e-4, -0.10302738212103e-5,
                                 0.28217281635040e-6, 0.12704902271945e-5, 0.73803353468292e-7,
                                 -0.11030139238909e-7, -0.81456365207833e-13, -0.25180545682962e-10,
                                 -0.17565233969407e-17, 0.86934156344163e-14};
                double pStar = 1E6;     /**< \f$ p^*  \f$ */
                double hStar = 2000E3;  /**< \f$ h^{*} \f$ */
                // double TStar = 1;    /**< \f$ T^* \f$ */
            }coeff_backward_T_PH_region2b;
            /**
             * @brief Numerical values of the coefficients and exponents of the backward equation T ( p, h ) for subregion 2c, Eq. (24). See Table 22 of \cite IF97.
             *
             */
            struct Coeff_Backward_T_PH_Region2c
            {
                static const int num=23;
                double I[num] = {-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6};
                double J[num] = {0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22};
                double n[num] = {-0.32368398555242e13, 0.73263350902181e13, 0.35825089945447e12,
                                 -0.58340131851590e12, -0.10783068217470e11, 0.20825544563171e11,
                                 0.61074783564516e6, 0.85977722535580e6, -0.25745723604170e5,
                                 0.31081088422714e5, 0.12082315865936e4, 0.48219755109255e3,
                                 0.37966001272486e1, -0.10842984880077e2, -0.45364172676660e-1,
                                 0.14559115658698e-12, 0.11261597407230e-11, -0.17804982240686e-10,
                                 0.12324579690832e-6, -0.11606921130984e-5, 0.27846367088554e-4,
                                 -0.59270038474176e-3, 0.12918582991878e-2};
                double pStar = 1E6;     /**< \f$ p^*  \f$ */
                double hStar = 2000E3;  /**< \f$ h^{*} \f$ */
                // double TStar = 1;    /**< \f$ T^* \f$ */
            }coeff_backward_T_PH_region2c;
            struct Coeff_Backward_T_PH_Region1
            {
                static const int num=20;
                double I[num] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6};
                double J[num] = {0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32};
                double n[num] = {-0.23872489924521e3, 0.40421188637945e3, 0.11349746881718e3,
                                 -0.58457616048039e1, -0.15285482413140e-3, -0.10866707695377e-5,
                                 -0.13391744872602e2, 0.43211039183559e2, -0.54010067170506e2,
                                 0.30535892203916e2, -0.65964749423638e1, 0.93965400878363e-2,
                                 0.11573647505340e-6, -0.25858641282073e-4, -0.40644363084799e-8,
                                 0.66456186191635e-7, 0.80670734103027e-10, -0.93477771213947e-12,
                                 0.58265442020601e-14, -0.15020185953503e-16};
                double pStar = 1E6;     /**< \f$ p^*  \f$ */
                double hStar = 2500E3;  /**< \f$ h^{*} \f$ */
                // double TStar = 1;    /**< \f$ T^* \f$ */
            }coeff_backward_T_PH_region1;
            /**
             * @brief Coefficients and exponents of the backward equation \f$ T_{3a}(p,h)\f$ for subregion 3a in its dimensionless form, Eq. (2) in \cite IF97-Region3
             *
             */
            struct Coeff_Backward_T_PH_Region3a
            {
                const int num = 31;
                double I[31] = {-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8,
                                -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12};
                double J[31] = {0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0,
                                2, 0, 1, 1, 0, 1, 0, 3, 4, 5};
                double n[31] = {-0.133645667811215e-6, 0.455912656802978e-5, -0.146294640700979e-4,
                                0.639341312970080e-2, 0.372783927268847e3, -0.718654377460447e4,
                                0.573494752103400e6, -0.267569329111439e7, -0.334066283302614e-4,
                                -0.245479214069597e-1, 0.478087847764996e2, 0.764664131818904e-5,
                                0.128350627676972e-2, 0.171219081377331e-1, -0.851007304583213e1,
                                -0.136513461629781e-1, -0.384460997596657e-5, 0.337423807911655e-2,
                                -0.551624873066791, 0.729202277107470, -0.992522757376041e-2,
                                -.119308831407288, .793929190615421, .454270731799386,
                                .20999859125991, -0.642109823904738e-2, -0.235155868604540e-1,
                                0.252233108341612e-2, -0.764885133368119e-2, 0.136176427574291e-1,
                                -0.133027883575669e-1};
                double pStar = 100E6;
                double TStar = 760;
                double hStar = 2300E3;
            }coeff_backward_T_PH_region3a;
            /**
             * @brief Coefficients and exponents of the backward equation \f$ T_{3a}(p,h)\f$ for subregion 3b in its dimensionless form, Eq. (3) in \cite IF97-Region3
             *
             */
            struct Coeff_Backward_T_PH_Region3b
            {
                const int num = 33;
                double I[33] = {-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4,
                                -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8};
                double J[33] = {0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4,
                                6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1};
                double n[33] = {0.323254573644920e-4, -0.127575556587181e-3, -0.475851877356068e-3,
                                0.156183014181602e-2, 0.105724860113781, -0.858514221132534e2,
                                0.724140095480911e3, 0.296475810273257e-2, -0.592721983365988e-2,
                                -0.126305422818666e-1, -0.115716196364853, 0.849000969739595e2,
                                -0.108602260086615e-1, 0.154304475328851e-1, 0.750455441524466e-1,
                                0.252520973612982e-1, -0.602507901232996e-1, -0.307622221350501e1,
                                -0.574011959864879e-1, 0.503471360939849e1, -0.925081888584834,
                                0.391733882917546e1, -0.773146007130190e2, 0.949308762098587e4,
                                -0.141043719679409e7, 0.849166230819026e7, 0.861095729446704,
                                0.323346442811720, 0.873281936020439, -0.436653048526683,
                                0.286596714529479, -0.131778331276228, 0.676682064330275e-2};
                double pStar = 100E6;
                double TStar = 860;
                double hStar = 2800E3;
            }coeff_backward_T_PH_region3b;

            /**
             * @brief Table 34 of \cite IF97 for saturation curve calculation.
             *
             */
            struct Table34
            {
                double n[10] = {0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02,
                                0.12020824702470E+05, -0.32325550322333E+07, 0.14915108613530E+02,
                                -0.48232657361591E+04, 0.40511340542057E+06, -0.23855557567849E+00,
                                0.65017534844798E+03};
            }table34;
            ThermodynamicProperties m_prop_l, m_prop_v;
            CONSTENTS_Thermo m_constants;
        public:
            cIAPWS_IF97(/* args */);
            ~cIAPWS_IF97();
        private:
            void initialize_data();

        public:
            double Boundary_region23_P2T(double P);
            double Boundary_region2b2c_P2H(double P);
            double Boundary_region2b2c_H2P(double H);
            double Boundary_region3ab_P2H(double P);
            double Boundary_region3ab_H2P(double H);
            void getState_Region1(double P, double T, State_Region1& state);
            double Prop_Region1(const State_Region1& state, BasicThermodynamicProperties which);
            void getState_Region2(double P, double T, State_Region2& state);
            double Prop_Region2(State_Region2 state, BasicThermodynamicProperties which);
            void getState_Region5(double P, double T, State_Region5& state);
            double Prop_Region5(double P, double T, State_Region5 state, BasicThermodynamicProperties which);
            int    GetRegion_PH(double P, double H);
            int    GetRegion_PT(double P, double T);
            double Backward_T_PH_region1(double P, double H);
            double Backward_T_PH_region2a(double P, double H);
            double Backward_T_PH_region2b(double P, double H);
            double Backward_T_PH_region2c(double P, double H);
            double Backward_T_PH_region3a(double P, double H);
            double Backward_T_PH_region3b(double P, double H);

            double T_sat_P(const double& P);
            void H_sat_P(const double& P, double& H_l, double& H_v);
            double P_sat_T(const double& T);
            double T_PH(const double& P, const double& H);

            // vector version
            void GetRegion_PT(const std::vector<double> P, const std::vector<double> T, std::vector<int>& res);
            void GetRegion_PH(const std::vector<double> P, const std::vector<double> H, std::vector<int>& res);
            void P_sat_T(const std::vector<double> T, std::vector<double>& res);
            void T_sat_P(const std::vector<double> P, std::vector<double>& res);
            void H_sat_P(const std::vector<double> P, std::vector<double>& H_l,std::vector<double>& H_v);
            void T_PH(const std::vector<double> P, const std::vector<double> H, std::vector<double>& res);
            void Boundary_region3ab_H2P(const std::vector<double> H, std::vector<double>& res);
            void Boundary_region3ab_P2H(const std::vector<double> P, std::vector<double>& res);
        public:
            std::string name(){return Name_Backend_IAPWS_IF97;};
            double Tmin(){return m_constants.Tmin; };                  /**< Get the minimum temperature in K */
            double Tmax(){return m_constants.Tmax;};                  /**< Get the maximum temperature in K */
            double pmin(){return m_constants.pmin; };                    /**< Get the minimum pressure in Pa */
            double pmax(){return m_constants.pmax;};                  /**< Get the maximum pressure in Pa */
            double Ttriple(){return m_constants.Ttriple;};               /**< Get the triple point temperature in K */
            double T_critical(){return m_constants.T_critical;};            /**< Return the critical temperature in K */
            double p_critical(){return m_constants.p_critical;};            /**< Return the critical pressure in Pa */
            double rhomass_critical(){return m_constants.rhomass_critical;};      /**< Return the critical mass density in \f$ kg/m^3 \f$ */
            double molar_mass(){return m_constants.molar_mass;};       /**< Return the molar mass in kg/mol */

            void UpdateState_TP(ThermodynamicProperties& props, const double& T, const double& p) {throw NotImplementedError(name() + " does not implement UpdateState_TP function");};
            double Boiling_p(const double& T) {throw NotImplementedError(name() + " does not implement Boiling_p(const double& T) function");};
        };
    
    };
};

#endif