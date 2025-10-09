/**
 * @file IAPWS95.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Implementation of IAPWS95 EOS.
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "IAPWS95.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

namespace xThermal
{
    namespace IAPWS95
    {
        void cIAPWS95::initialize_data()
        {
            m_constants.Ttriple = 273.16;
            m_constants.Tmin = 273.165; // see section 5 Range of validity (pp.5) of IAPWS95 release., but I set it to 0 deg.C
            m_constants.Tmax = 2273.15;
            m_constants.pmin = 208.566; // see section 5 Range of validity (pp.5) of IAPWS95 release.
            m_constants.pmax = 5000E5;
            m_constants.T_critical = 647.096;
            m_constants.p_critical = 22.064E6;
            m_constants.rhomass_critical = 322.0;
            m_constants.molar_mass = 0.018015268; //kg/mol
        }
        cIAPWS95::cIAPWS95(/* args */)
        {
            initialize_data();
        }

        cIAPWS95::~cIAPWS95()
        {
        }

        /**
     * @brief The ideal-gas part \f$ \phi^o \f$ of the dimensionless Helmholtz free energy and its derivatives. See Table 4 in \cite IAPWS-95
     *
     * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
     * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
     * @param phio \f$ \phi^o \f$
     * @param update_derivatives Combination of the bitmask value: #Update_phi_d, #Update_phi_dd, #Update_phi_t, #Update_phi_tt, #Update_phi_dt, #Update_phi_all
     */
        void cIAPWS95::phi_o(const double& delta, const double& tau, HelmholtzEnergy_dimensionless& phio, unsigned int update_derivatives)
        {
            double one_minus_expGammaTau[5];
            phio.value = 0;
            // terms 1
            phio.value += log(delta) + m_coeff_phio.n0_term1[0] + m_coeff_phio.n0_term1[1]*tau + m_coeff_phio.n0_term1[2]*log(tau);
            // term 2
            for (int i = 0; i < m_coeff_phio.n2; i++)
            {
                one_minus_expGammaTau[i] = 1 - exp(-m_coeff_phio.gamma0_term2[i]*tau);
                phio.value += m_coeff_phio.n0_term2[i] * log(one_minus_expGammaTau[i]);
            }
            // calculate derivertives of phio
            // 1. ∂φ0/∂δ
            if((update_derivatives & Update_phi_d) == Update_phi_d)
            {
                phio.d = 1.0/delta;
            }
            // 2. ∂2φ0/∂δ2
            if((update_derivatives & Update_phi_dd) == Update_phi_dd)
            {
                phio.dd = -1.0/(delta * delta);
            }
            // 3. ∂φ0/∂τ
            if((update_derivatives & Update_phi_t) == Update_phi_t)
            {
                phio.t = m_coeff_phio.n0_term1[1] + m_coeff_phio.n0_term1[2]/tau;
                for (int i = 0; i < m_coeff_phio.n2; i++)
                {
                    phio.t += m_coeff_phio.n0_term2[i] * m_coeff_phio.gamma0_term2[i] * (1.0/one_minus_expGammaTau[i] - 1.0);
                }
            }
            // 4. ∂2φ0/∂τ2
            if((update_derivatives & Update_phi_tt) == Update_phi_tt)
            {
                phio.tt = -m_coeff_phio.n0_term1[2]/(tau*tau);
                for (int i = 0; i < m_coeff_phio.n2; i++)
                {
                    phio.tt -= m_coeff_phio.n0_term2[i] * m_coeff_phio.gamma0_term2[i] * m_coeff_phio.gamma0_term2[i] * exp(-m_coeff_phio.gamma0_term2[i]*tau)/(one_minus_expGammaTau[i]*one_minus_expGammaTau[i]);
                }
            }
            // ∂2φ0/∂τ∂δ = 0
            phio.dt = 0;
        }

        /**
         * @brief Implementation of residual part \f$ \phi^r \f$ of the dimensionless Helmholtz free energy: Eq. 6 in \cite IAPWS-95 and its derivatives, see Table 5 in \cite IAPWS-95
         * \f{matrix}
             * \phi^r & = & \sum\limits_{i=1}^{7} n_i \delta^{d_i}\tau^{t_i} + \sum\limits_{i=8}^{51} n_i \delta^{d_i} e^{-\delta^{c_i}}  \tau^{t_i} + \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} + \sum\limits_{i=55}^{56} n_i \delta \Delta^{b_i}\psi
             * , \left(
             * \Delta  =  \theta^2 + B_i [\color{red}{(\delta -1)^2}]^{a_i} ,
             * \theta  =  (1-\tau) + A_i [\color{red}{(\delta -1)^2}]^{1/2\beta_i} ,
             * \psi  =  e^{-C_i (\delta -1)^2 - D_i (\tau - 1)^2}
             * \right)
             * \f}
             *
         * \image html IAPWS96_phi/phir_der_T300K.svg "The residual part of Helmholtz free energy and its derivatives." width=100%.
         * \note The red dashed line is calculated from a python package of IAPWS \cite pyIAPWS.
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r(const double& delta, const double& tau)
        {
            double phir = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir += m_coeff_phir.ni_term1[i] * pow(delta, m_coeff_phir.di_term1[i]) * pow(tau, m_coeff_phir.ti_term1[i]);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir += m_coeff_phir.ni_term2[i] * pow(delta, m_coeff_phir.di_term2[i]) * pow(tau, m_coeff_phir.ti_term2[i]) * exp(-pow(delta,m_coeff_phir.ci_term2[i]));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir += m_coeff_phir.ni_term3[i] * pow(delta, m_coeff_phir.d_term3) * pow(tau, m_coeff_phir.ti_term3[i]) * exp(-m_coeff_phir.alpha_term3 * pow(delta - m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i] * pow(tau - m_coeff_phir.gammai_term3[i], 2.0));
            }
            // terms 4: Nonanalytical term
            double Delta, theta, psi, delta_1_sqr;
            for (int i = 0; i < m_coeff_phir.n4; i++)
            {
                delta_1_sqr = (delta-1)*(delta-1);
                theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
                Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
                psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));

                phir += m_coeff_phir.ni_term4[i] * pow(Delta, m_coeff_phir.bi_term4[i]) * delta * psi;
            }
            return phir;
        }
        void cIAPWS95::phi_r(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r(delta[i], tau[i]);
        }

        /**
         * @brief Calculate partial derivative of residual part \f$ \left[ \frac{\partial \phi^r}{\partial \delta} \right]_{\tau} \f$
         *
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r_d(const double& delta, const double& tau)
        {
            double phir_d = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir_d += m_coeff_phir.ni_term1[i]*m_coeff_phir.di_term1[i]*pow(delta, m_coeff_phir.di_term1[i]-1)*pow(tau, m_coeff_phir.ti_term1[i]);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir_d += m_coeff_phir.ni_term2[i]*exp(-pow(delta, m_coeff_phir.ci_term2[i]))*(pow(delta, m_coeff_phir.di_term2[i]-1)*pow(tau, m_coeff_phir.ti_term2[i])*(m_coeff_phir.di_term2[i] - m_coeff_phir.ci_term2[i]*pow(delta, m_coeff_phir.ci_term2[i])));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir_d += m_coeff_phir.ni_term3[i]*pow(delta, m_coeff_phir.d_term3)*pow(tau, m_coeff_phir.ti_term3[i])*exp(-m_coeff_phir.alpha_term3*pow(delta - m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i]*pow(tau - m_coeff_phir.gammai_term3[i], 2.0))*(m_coeff_phir.d_term3/delta - 2.0*m_coeff_phir.alpha_term3*(delta - m_coeff_phir.epsilon_term3));
            }
            // terms 4: Nonanalytical term
            double Delta, theta, psi, delta_1_sqr, psi_d, Delta_bi_d, Delta_d;
            for (int i = 0; i < m_coeff_phir.n4; i++)
            {
                delta_1_sqr = (delta-1)*(delta-1);
                theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
                Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
                psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));
                psi_d = -2.0*m_coeff_phir.Ci_term4[i]*(delta-1)*psi;
                Delta_d = (delta - 1)*(m_coeff_phir.A_term4 * theta * 2.0/m_coeff_phir.beta_term4*pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4 - 1.0) + 2.0*m_coeff_phir.B_term4*m_coeff_phir.a_term4*pow(delta_1_sqr, m_coeff_phir.a_term4-1.0));
                Delta_bi_d = m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i]-1.0)*Delta_d;

                phir_d += m_coeff_phir.ni_term4[i]*(pow(Delta, m_coeff_phir.bi_term4[i])*(psi + delta*psi_d) + Delta_bi_d * delta * psi);
            }
            return phir_d;
        }
        void cIAPWS95::phi_r_d(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r_d(delta[i], tau[i]);
        }

        /**
         * @brief Calculate partial derivative of residual part \f$ \left[ \frac{\partial^2 \phi^r}{\partial \delta^2} \right]_{\tau} \f$
         *
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r_dd(const double& delta, const double& tau)
        {
            double phir_dd = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir_dd += m_coeff_phir.ni_term1[i]*m_coeff_phir.di_term1[i]*(m_coeff_phir.di_term1[i] - 1)*pow(delta, m_coeff_phir.di_term1[i]-2)*pow(tau, m_coeff_phir.ti_term1[i]);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir_dd += m_coeff_phir.ni_term2[i]*exp(-pow(delta, m_coeff_phir.ci_term2[i]))*(pow(delta, m_coeff_phir.di_term2[i]-2)*pow(tau, m_coeff_phir.ti_term2[i])*((m_coeff_phir.di_term2[i] - m_coeff_phir.ci_term2[i]*pow(delta, m_coeff_phir.ci_term2[i]))*(m_coeff_phir.di_term2[i] - 1.0 - m_coeff_phir.ci_term2[i]*pow(delta, m_coeff_phir.ci_term2[i])) - m_coeff_phir.ci_term2[i]*m_coeff_phir.ci_term2[i]*pow(delta, m_coeff_phir.ci_term2[i])));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir_dd += m_coeff_phir.ni_term3[i]*pow(tau, m_coeff_phir.ti_term3[i])*exp(-m_coeff_phir.alpha_term3*pow(delta - m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i]*pow(tau - m_coeff_phir.gammai_term3[i], 2.0)) * (-2.0*m_coeff_phir.alpha_term3*pow(delta, m_coeff_phir.d_term3) + 4.0*pow(m_coeff_phir.alpha_term3, 2.0)*pow(delta, m_coeff_phir.d_term3)*pow(delta - m_coeff_phir.epsilon_term3, 2.0) - 4.0*m_coeff_phir.d_term3*m_coeff_phir.alpha_term3*pow(delta, m_coeff_phir.d_term3-1)*(delta-m_coeff_phir.epsilon_term3) + m_coeff_phir.d_term3*(m_coeff_phir.d_term3-1)*pow(delta, m_coeff_phir.d_term3-2));
            }
            // terms 4: Nonanalytical term
            double Delta, theta, psi, delta_1_sqr, psi_d, psi_dd, Delta_bi_d, Delta_bi_dd, Delta_d, Delta_dd;
            for (int i = 0; i < m_coeff_phir.n4; i++)
            {
                delta_1_sqr = (delta-1)*(delta-1);
                theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
                Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
                psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));
                psi_d = -2.0*m_coeff_phir.Ci_term4[i]*(delta-1)*psi;
                Delta_d = (delta - 1)*(m_coeff_phir.A_term4 * theta * 2.0/m_coeff_phir.beta_term4*pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4 - 1.0) + 2.0*m_coeff_phir.B_term4*m_coeff_phir.a_term4*pow(delta_1_sqr, m_coeff_phir.a_term4-1.0));
                Delta_bi_d = m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i]-1.0)*Delta_d;
                psi_dd = (2.0*m_coeff_phir.Ci_term4[i]*pow(delta - 1.0, 2.0) - 1.0)*2.0*m_coeff_phir.Ci_term4[i]*psi;
                Delta_dd = 1.0/(delta-1.0)*Delta_d + pow(delta-1, 2.0)*(4.0*m_coeff_phir.B_term4*m_coeff_phir.a_term4*(m_coeff_phir.a_term4-1)*pow(delta_1_sqr, m_coeff_phir.a_term4-2) + 2.0*pow(m_coeff_phir.A_term4, 2.0)*pow(m_coeff_phir.beta_term4, -2)*pow(pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4-1), 2.0) + m_coeff_phir.A_term4*theta*4/m_coeff_phir.beta_term4*(0.5/m_coeff_phir.beta_term4 - 1.0)*pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4 - 2));
                Delta_bi_dd = m_coeff_phir.bi_term4[i] * (pow(Delta, m_coeff_phir.bi_term4[i]-1)*Delta_dd + (m_coeff_phir.bi_term4[i]-1)*pow(Delta, m_coeff_phir.bi_term4[i]-2)*pow(Delta_d, 2.0));

                phir_dd += m_coeff_phir.ni_term4[i]*(pow(Delta, m_coeff_phir.bi_term4[i])*(2*psi_d + delta*psi_dd) + 2.0*Delta_bi_d*(psi + delta*psi_d) + Delta_bi_dd * delta * psi);
            }
            return phir_dd;
        }
        void cIAPWS95::phi_r_dd(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r_dd(delta[i], tau[i]);
        }

        /**
         * @brief Calculate partial derivative of residual part \f$ \left[ \frac{\partial \phi^r}{\partial \tau} \right]_{\delta} \f$
         *
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r_t(const double& delta, const double& tau)
        {
            double phir_t = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir_t += m_coeff_phir.ni_term1[i]*m_coeff_phir.ti_term1[i]*pow(delta, m_coeff_phir.di_term1[i])*pow(tau, m_coeff_phir.ti_term1[i] - 1.0);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir_t += m_coeff_phir.ni_term2[i]*m_coeff_phir.ti_term2[i]*pow(delta, m_coeff_phir.di_term2[i])*pow(tau, m_coeff_phir.ti_term2[i] - 1.0)*exp(-pow(delta, m_coeff_phir.ci_term2[i]));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir_t += m_coeff_phir.ni_term3[i]*pow(delta, m_coeff_phir.d_term3)*pow(tau, m_coeff_phir.ti_term3[i])*exp(-m_coeff_phir.alpha_term3*pow(delta - m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i]*pow(tau-m_coeff_phir.gammai_term3[i], 2.0)) * (m_coeff_phir.ti_term3[i]/tau - 2.0*m_coeff_phir.betai_term3[i]*(tau - m_coeff_phir.gammai_term3[i]));
            }
            // terms 4: Nonanalytical term
            double Delta, theta, psi, delta_1_sqr, psi_t, Delta_bi_t;
            for (int i = 0; i < m_coeff_phir.n4; i++)
            {
                delta_1_sqr = (delta-1)*(delta-1);
                theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
                Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
                psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));
                psi_t = -2.0*m_coeff_phir.Di_term4[i]*(tau - 1.0)*psi;
                Delta_bi_t = -2.0*theta*m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i] - 1.0);

                phir_t += m_coeff_phir.ni_term4[i]*delta*(Delta_bi_t*psi + pow(Delta, m_coeff_phir.bi_term4[i])*psi_t);
            }
            return phir_t;
        }
        void cIAPWS95::phi_r_t(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r_t(delta[i], tau[i]);
        }

        /**
         * @brief Calculate partial derivative of residual part \f$ \left[ \frac{\partial^2 \phi^r}{\partial \tau^2} \right]_{\delta} \f$
         *
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r_tt(const double& delta, const double& tau)
        {
            double phir_tt = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir_tt += m_coeff_phir.ni_term1[i]*m_coeff_phir.ti_term1[i]*(m_coeff_phir.ti_term1[i]-1.0)*pow(delta, m_coeff_phir.di_term1[i])*pow(tau, m_coeff_phir.ti_term1[i]-2.0);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir_tt += m_coeff_phir.ni_term2[i]*m_coeff_phir.ti_term2[i]*(m_coeff_phir.ti_term2[i]-1.0)*pow(delta, m_coeff_phir.di_term2[i])*pow(tau, m_coeff_phir.ti_term2[i]-2.0)*exp(-pow(delta, m_coeff_phir.ci_term2[i]));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir_tt += m_coeff_phir.ni_term3[i]*pow(delta, m_coeff_phir.d_term3)*pow(tau, m_coeff_phir.ti_term3[i])*exp(-m_coeff_phir.alpha_term3*pow(delta-m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i]*pow(tau-m_coeff_phir.gammai_term3[i], 2.0))*(pow(m_coeff_phir.ti_term3[i]/tau - 2.0*m_coeff_phir.betai_term3[i]*(tau - m_coeff_phir.gammai_term3[i]), 2.0) - m_coeff_phir.ti_term3[i]/pow(tau, 2.0) - 2.0*m_coeff_phir.betai_term3[i]);
            }
            // // terms 4: Nonanalytical term
            // double Delta, theta, psi, delta_1_sqr, psi_t, psi_tt, Delta_bi_t, Delta_bi_tt;
            // for (int i = 0; i < m_coeff_phir.n4; i++)
            // {
            //     delta_1_sqr = (delta-1)*(delta-1);
            //     theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
            //     Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
            //     psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));
            //     psi_t = -2.0*m_coeff_phir.Di_term4[i]*(tau - 1.0)*psi;
            //     psi_tt = (2.0*m_coeff_phir.Di_term4[i]*pow(tau-1.0, 2.0) - 1.0)*2.0*m_coeff_phir.Di_term4[i]*psi;
            //     Delta_bi_t = -2.0*theta*m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i] - 1.0);
            //     Delta_bi_tt = 2.0*m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i]-1.0) + 2.0*pow(theta, 2.0)*m_coeff_phir.bi_term4[i]*(m_coeff_phir.bi_term4[i]-1.0)*pow(Delta, m_coeff_phir.bi_term4[i]-2.0);

            //     phir_tt += m_coeff_phir.ni_term4[i]*delta*(Delta_bi_tt*psi + 2.0*Delta_bi_t*psi_t + pow(Delta, m_coeff_phir.bi_term4[i])*psi_tt);
            // }
            return phir_tt;
        }
        void cIAPWS95::phi_r_tt(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r_tt(delta[i], tau[i]);
        }

        /**
         * @brief Calculate partial derivative of residual part \f$ \frac{\partial^2 \phi^r}{\partial \delta \partial \tau} \f$
         *
         * @param delta \f$ \delta = \rho / \rho_{critical} \f$ is the reduced density.
         * @param tau \f$ \tau = T_{critical}/T \f$ is the inverse reduced temperature. \note The unit of temperature is \f$ K \f$. The unit of density if \f$ kg/m^3\f$
         * @return double
         */
        double cIAPWS95::phi_r_dt(const double& delta, const double& tau)
        {
            double phir_dt = 0;
            // terms 1: Polinomial
            for (int i = 0; i < m_coeff_phir.n1; i++)
            {
                phir_dt += m_coeff_phir.ni_term1[i]*m_coeff_phir.di_term1[i]*m_coeff_phir.ti_term1[i]*pow(delta, m_coeff_phir.di_term1[i]-1.0)*pow(tau, m_coeff_phir.ti_term1[i]-1.0);
            }
            // terms 2: Exponential
            for (int i = 0; i < m_coeff_phir.n2; i++)
            {
                phir_dt += m_coeff_phir.ni_term2[i]*m_coeff_phir.ti_term2[i]*pow(delta, m_coeff_phir.di_term2[i]-1.0)*pow(tau, m_coeff_phir.ti_term2[i]-1.0)*(m_coeff_phir.di_term2[i] - m_coeff_phir.ci_term2[i]*pow(delta, m_coeff_phir.ci_term2[i]))*exp(-pow(delta, m_coeff_phir.ci_term2[i]));
            }
            // terms 3: Gaussian Bell-shaped term
            for (int i = 0; i < m_coeff_phir.n3; i++)
            {
                phir_dt += m_coeff_phir.ni_term3[i]*pow(delta, m_coeff_phir.d_term3)*pow(tau, m_coeff_phir.ti_term3[i])*exp(-m_coeff_phir.alpha_term3*pow(delta-m_coeff_phir.epsilon_term3, 2.0) - m_coeff_phir.betai_term3[i]*pow(tau-m_coeff_phir.gammai_term3[i], 2.0))*(m_coeff_phir.d_term3/delta - 2.0*m_coeff_phir.alpha_term3*(delta - m_coeff_phir.epsilon_term3))*(m_coeff_phir.ti_term3[i]/tau - 2.0*m_coeff_phir.betai_term3[i]*(tau - m_coeff_phir.gammai_term3[i]));
            }
            // terms 4: Nonanalytical term
            double Delta, theta, psi, delta_1_sqr, psi_t, psi_d, psi_dt, Delta_d, Delta_bi_t, Delta_bi_d, Delta_bi_dt;
            for (int i = 0; i < m_coeff_phir.n4; i++)
            {
                delta_1_sqr = (delta-1)*(delta-1);
                theta = (1-tau) + m_coeff_phir.A_term4 * pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4);
                Delta = theta*theta + m_coeff_phir.B_term4 * pow(delta_1_sqr, m_coeff_phir.a_term4);
                psi   = exp(-m_coeff_phir.Ci_term4[i]*delta_1_sqr - m_coeff_phir.Di_term4[i] * (tau-1)*(tau-1));
                psi_t = -2.0*m_coeff_phir.Di_term4[i]*(tau - 1.0)*psi;
                Delta_bi_t = -2.0*theta*m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i] - 1.0);
                psi_d = -2.0*m_coeff_phir.Ci_term4[i]*(delta-1)*psi;
                Delta_d = (delta - 1)*(m_coeff_phir.A_term4 * theta * 2.0/m_coeff_phir.beta_term4*pow(delta_1_sqr, 0.5/m_coeff_phir.beta_term4 - 1.0) + 2.0*m_coeff_phir.B_term4*m_coeff_phir.a_term4*pow(delta_1_sqr, m_coeff_phir.a_term4-1.0));
                Delta_bi_d = m_coeff_phir.bi_term4[i]*pow(Delta, m_coeff_phir.bi_term4[i]-1.0)*Delta_d;
                psi_dt = 4.0*m_coeff_phir.Ci_term4[i]*m_coeff_phir.Di_term4[i]*(delta - 1.0)*(tau - 1.0)*psi;
                Delta_bi_dt = -m_coeff_phir.A_term4*m_coeff_phir.bi_term4[i]*2.0/m_coeff_phir.beta_term4*pow(Delta, m_coeff_phir.bi_term4[i]-1.0)*(delta-1.0)*pow(pow(delta-1, 2.0), 0.5/m_coeff_phir.beta_term4-1.0) - 2.0*theta*m_coeff_phir.bi_term4[i]*(m_coeff_phir.bi_term4[i]-1.0)*pow(Delta, m_coeff_phir.bi_term4[i]-2.0)*Delta_d;

                phir_dt += m_coeff_phir.ni_term4[i]*(pow(Delta, m_coeff_phir.bi_term4[i]) * (psi_t + delta*psi_dt) + delta*Delta_bi_d*psi_t + Delta_bi_t*(psi + delta*psi_d) + Delta_bi_dt*delta*psi);
            }
            return phir_dt;
        }
        void cIAPWS95::phi_r_dt(const std::vector<double> delta, const std::vector<double> tau, std::vector<double>& res)
        {
            res.clear();
            res.resize(tau.size());
            for (size_t i = 0; i < tau.size(); i++)res[i] = phi_r_dt(delta[i], tau[i]);
        }

        void cIAPWS95::phi_r(const double& delta, const double& tau, HelmholtzEnergy_dimensionless& phir)
        {
            phir.value  = phi_r(delta, tau);
            phir.d      = phi_r_d(delta, tau);
            phir.dd     = phi_r_dd(delta, tau);
            phir.t      = phi_r_t(delta, tau);
            phir.tt     = phi_r_tt(delta, tau);
            phir.dt     = phi_r_dt(delta, tau);
        }

        /**
     * @brief Vapor-pressure equation. See Eq. 2.5 (pp. 399) in \cite IAPWS-95-art
     * \note This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. \cite IAPWS-95-art.
     * \warning The unit of this function depends on unit of H2O::P_c, I set critical pressure H2O::P_c to Pa, so this function return pressure in Pa.
     * @param T_K [K]
     * @return double [Pa]
     */
        double cIAPWS95::P_Sat_estimate(const double& T_K)
        {
            double tau = m_constants.T_critical/T_K;
            double theta = (1.0 - T_K/m_constants.T_critical);
            double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};

            return exp(tau * (a[0] * theta + a[1] * pow(theta,1.5) + a[2] * pow(theta, 3.0) + a[3] * pow(theta, 3.5) + a[4] * pow(theta,4.0) + a[5] * pow(theta ,7.5))) * m_constants.p_critical;
        }

        /**
         * @brief Vapor-pressure and its derivative equation. See Eq. 2.5, 2.5a (pp. 399) in \cite IAPWS-95-art
         *
         * @param T_K
         * @param p
         * @param dpdT
         * @return double
         */
        void cIAPWS95::P_Sat_estimate(const double& T_K, double& p, double& dpdT)
        {
            double tau = m_constants.T_critical/T_K;
            double theta = (1.0 - T_K/m_constants.T_critical);
            double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
            double log_p_pc = tau * (a[0] * theta + a[1] * pow(theta,1.5) + a[2] * pow(theta, 3.0) + a[3] * pow(theta, 3.5) + a[4] * pow(theta,4.0) + a[5] * pow(theta ,7.5));
            p = exp(log_p_pc)*m_constants.p_critical;
            // dpdT
            dpdT = -p/T_K * (log_p_pc + a[0] + 1.5*a[1]*pow(theta, 0.5) + 3.0*a[2]*pow(theta, 2.0) + 3.5*a[3]*pow(theta, 2.5) + 4.0*a[4]*pow(theta, 3.0) + 7.5*a[5]*pow(theta, 6.5));
        }

        double func_T_Sat_estimate(double T_K, void *params)
        {
            Params_T_Sat_estimate* param = (Params_T_Sat_estimate*)params;
            cIAPWS95* iapws = param->iapws;
            // printf("x: %f\n",T_K);
            return iapws->P_Sat_estimate(T_K) - param->p;
        }

        /**
         * @brief Solve nonlinear equation of #P_Sat_estimate
         *
         * @param P
         * @return double
         */
        double cIAPWS95::T_Sat_estimate(const double& P)
        {
            // Could use IF97 equation instead
            return m_IF97.T_sat_P(P);
            // double T_K = 460.0; //use center value between 273.15 and T_c as initial guess
            // int status;
            // int iter = 0;
            // const gsl_root_fsolver_type *T;
            // gsl_root_fsolver *s;
            // double x_lo = T_MIN, x_hi = T_c;
            // gsl_function F;
            // Params_T_Sat_estimate params={this, P};

            // F.function = &func_T_Sat_estimate;
            // F.params = &params;

            // T = gsl_root_fsolver_brent;
            // s = gsl_root_fsolver_alloc (T);
            // gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
            // do
            // {
            //     iter++;
            //     status = gsl_root_fsolver_iterate (s);
            //     T_K = gsl_root_fsolver_root (s);
            //     x_lo = gsl_root_fsolver_x_lower (s);
            //     x_hi = gsl_root_fsolver_x_upper (s);
            //     status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure);
            //     // if (status == GSL_SUCCESS) printf ("Converged:\n");
            //     // printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,T_K);
            // }
            // while (status == GSL_CONTINUE && iter < ITERATION_MAX);

            // if(!((status == GSL_SUCCESS)))
            // {
            //     printf ("status = %s\n\n", gsl_strerror (status));
            //     printf("P = %.3E Pa\n", P);
            //     ERROR("Fatal error in cIAPWS95::T_Sat_estimate(double P)");
            // }
            // gsl_root_fsolver_free (s);

            // return T_K;
        }

        /**
     * @brief Saturated liquid density. See Eq. 2.6 (pp. 399) in \cite IAPWS-95-art
     * \note This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. \cite IAPWS-95-art.
     * @param T_K [K]
     * @return double [kg/m3]
     */
        double cIAPWS95::Rho_Liquid_Sat_estimate(const double& T_K)
        {
            double tau = m_constants.T_critical/T_K;
            double theta = (1.0 - T_K/m_constants.T_critical);
            double b[6] = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45};

            return (1
                    + b[0] * pow(theta, 1/3.0)
                    + b[1] * pow(theta, 2/3.0)
                    + b[2] * pow(theta, 5/3.0)
                    + b[3] * pow(theta, 16/3.0)
                    + b[4] * pow(theta, 43/3.0)
                    + b[5] * pow(theta, 110/3.0)) * m_constants.rhomass_critical;
        }

        /**
         * @brief Saturated vapor density. See Eq. 2.7 (pp. 399) in \cite IAPWS-95-art
         * \note This equation is only auxiliary equation for calculating properties along the vapor–liquid phase boundary. Although the differences be- tween the results from these equations and the corresponding results from the IAPWS-95 formulation, Eq. 6.4, are extremely small, these equations are not thermodynamically consistent with IAPWS-95. Nevertheless, the application of these auxiliary equations might be useful for simple calculations of these saturation properties or for the determination of starting values when iteratively calculating the saturation properties from IAPWS-95 according to the phase-equilibrium condition. \cite IAPWS-95-art.
         * @param T_K [K]
         * @return double [kg/m3]
         */
        double cIAPWS95::Rho_Vapor_Sat_estimate(const double& T_K)
        {
            double tau = m_constants.T_critical/T_K;
            double theta = (1.0 - T_K/m_constants.T_critical);
            double c[6] = {-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581, -63.9201063};
            double rho_V_sat =0;
            return exp(c[0] * pow(theta, 1/3.0)  +
                       c[1] * pow(theta, 2/3.0)  +
                       c[2] * pow(theta, 4/3.0)  +
                       c[3] * pow(theta, 3.0)    +
                       c[4] * pow(theta, 37/6.0) +
                       c[5] * pow(theta, 71/6.0)) * m_constants.rhomass_critical;
        }

        void print_state_TP2Rho(size_t iter, gsl_multiroot_fsolver * s)
        {
            printf ("iter = %3lu Rho = % 15.8f,  f(x) = % .3e \n",iter,
                    gsl_vector_get (s->x, 0), gsl_vector_get (s->f, 0));
        }

        /**
         * @brief Construct nonlinear functions according to phase-equilibrium condition. See Table 3 (pp.12) in \cite IAPWS-95
         *
         * @param x
         * @param params
         * @param f
         * @return int
         */
        int func_PhaseEquilibrium(const gsl_vector * x, void *params, gsl_vector * f)
        {
            Params_SolvePhaseEquilibrium* param = (Params_SolvePhaseEquilibrium*)params;
            cIAPWS95* iapws = param->iapws;
            double y0, y1, y2;
            switch (param->Solve_PorT)
            {
                case SOLVE_SATURATED_P:
                {
                    double tau = param->TorP.tau;
                    double p        = gsl_vector_get (x, 0);
                    const double rho_l    = gsl_vector_get (x, 1);
                    const double rho_v    = gsl_vector_get (x, 2);
                    const double delta_l  = rho_l/iapws->m_constants.rhomass_critical;
                    const double delta_v  = rho_v/iapws->m_constants.rhomass_critical;
                    y0     = p/(param->RT*rho_l) - 1.0 - delta_l * iapws->phi_r_d(delta_l, tau);
                    y1     = p/(param->RT*rho_v) - 1.0 - delta_v * iapws->phi_r_d(delta_v, tau);
                    y2     = p/param->RT*(1.0/rho_v - 1.0/rho_l) - log(rho_l/rho_v) - iapws->phi_r(delta_l, tau) + iapws->phi_r(delta_v, tau);
                }
                    break;
                case SOLVE_SATURATED_T:
                {
                    double p = param->TorP.P;
                    double T_K        = gsl_vector_get (x, 0);
                    double tau          = iapws->m_constants.T_critical/T_K;
                    double RT           = iapws->m_constants.R*T_K;
                    const double rho_l    = gsl_vector_get (x, 1);
                    const double rho_v    = gsl_vector_get (x, 2);
                    const double delta_l  = rho_l/iapws->m_constants.rhomass_critical;
                    const double delta_v  = rho_v/iapws->m_constants.rhomass_critical;
                    y0     = p/(RT*rho_l) - 1.0 - delta_l * iapws->phi_r_d(delta_l, tau);
                    y1     = p/(RT*rho_v) - 1.0 - delta_v * iapws->phi_r_d(delta_v, tau);
                    y2     = p/RT*(1.0/rho_v - 1.0/rho_l) - log(rho_l/rho_v) - iapws->phi_r(delta_l, tau) + iapws->phi_r(delta_v, tau);
                }
                    break;

                default:
                ERROR("func_PhaseEquilibrium: param->Solve_PorT is not one of SOLVE_SATURATED_P, SOLVE_SATURATED_T: "+std::to_string(param->Solve_PorT));
                    break;
            }

            gsl_vector_set (f, 0, y0);
            gsl_vector_set (f, 1, y1);
            gsl_vector_set (f, 2, y2);

            return GSL_SUCCESS;
        }
        void print_state_PhaseEquilibrium(size_t iter, gsl_multiroot_fsolver * s)
        {
            printf ("iter = %3lu x = % 15.8f, rho_l=% 15.8f, rho_v=%15.8f,  f(x) = % .3e % .3e %.3e\n",iter,
                    gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->x, 2),
                    gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1), gsl_vector_get (s->f, 2));
        }

        /**
         * @brief Solve saturated pressure, liquid density and vapor density on phase boundary of boiling curve by given temperature.
         * Using #P_Sat_estimate, #Rho_Liquid_Sat_estimate and #Rho_Vapor_Sat_estimate calculate initial values of \f$ p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime} \f$, respectively; then use GSL solve three nonlinear quations with three unknowns (\f$ p_{\sigma}, \rho^{\prime}, \rho^{\prime\prime} \f$) according to the phase-equilibrium condition (See Table 3, pp. 11 in \cite IAPWS-95).
         * \note The pressure unit in the phase-equilibrium condition(Maxwell criterion) is Pa only if xThermal::R unit is J/kg/K !!! see Table 3 (pp.11) \cite IAPWS-95. I set xThermal::R in unit of J/kg/K.
         * @param T_K
         * @param P
         * @param rho_l
         * @param rho_v
         */
        void cIAPWS95::Boiling_p(const double& T_K, double& P, double& rho_l, double& rho_v)
        {
            if(T_K == m_constants.T_critical)
            {
                P = m_constants.p_critical;
                rho_l = m_constants.rhomass_critical;
                rho_v = m_constants.rhomass_critical;
                return;
            }
            double P0 = P_Sat_estimate(T_K);
            double Rho_l0 = Rho_Liquid_Sat_estimate(T_K);
            double Rho_v0 = Rho_Vapor_Sat_estimate(T_K);

            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            int status;
            size_t i, iter = 0;
            const size_t n = 3;
            Params_SolvePhaseEquilibrium params={this, T_K*m_constants.R, m_constants.T_critical/T_K, SOLVE_SATURATED_P};
            gsl_multiroot_function f = {&func_PhaseEquilibrium, n, &params};
            double x_init[n] = {P0, Rho_l0, Rho_v0}; //initial guess
            gsl_vector *x = gsl_vector_alloc (n);
            gsl_vector_set (x, 0, x_init[0]);
            gsl_vector_set (x, 1, x_init[1]);
            gsl_vector_set (x, 2, x_init[2]);
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc (T, n);
            gsl_multiroot_fsolver_set (s, &f, x);
            // print_state_PhaseEquilibrium (iter, s);
            // iteration: solve
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);
                // print_state_PhaseEquilibrium (iter, s);
                if (status) /* check if solver is stuck */
                    break;
                status = gsl_multiroot_test_residual (s->f, TOL_PTRro);
            }while (status == GSL_CONTINUE && iter < ITERATION_MAX);
            // printf ("status = %s, iter: %d\n", gsl_strerror (status),iter);
            if(!((status == GSL_SUCCESS))) // || (status == GSL_ENOPROG)
            {
                if(status == GSL_ENOPROG)
                {
                    print_state_TP2Rho (iter, s);
                    printf("T = %.3f K\n", T_K);
                    WARNING("Boiling_p: " + std::string(gsl_strerror (status)));
                }else
                {
                    print_state_TP2Rho (iter, s);
                    printf ("status = %s\n\n", gsl_strerror (status));
                    printf("T = %.3f K\n", T_K);
                    ERROR("Fatal error in cIAPWS95::Boiling_P(const double T_K, double& P, double& rho_l, double& rho_v)");
                }
            }
            P = gsl_vector_get (s->x, 0);
            rho_l = gsl_vector_get (s->x, 1);
            rho_v = gsl_vector_get (s->x, 2);
            // release pointer
            gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
        }

        void cIAPWS95::Boiling_T(const double& P, double& T_K, double& rho_l, double& rho_v)
        {
            //critical point
            if(P == m_constants.p_critical)
            {
               T_K = m_constants.T_critical;
               rho_l = m_constants.rhomass_critical;
               rho_v = m_constants.rhomass_critical;
               return;
            }
            // STATUS("Start calculate");
            // T_K = T_Sat_estimate(P); //have problem when P=P_MIN
            T_K = m_IF97.T_sat_P(P); //USE IAPWS-IF97 equation calculate a initial guess of T
            // printf("estimated T = %f\n", T_K);
            // T_K = T_Sat_estimate_df(P); //unreliable
            double Rho_l0 = Rho_Liquid_Sat_estimate(T_K);
            double Rho_v0 = Rho_Vapor_Sat_estimate(T_K);

            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            int status;
            size_t i, iter = 0;
            const size_t n = 3;
            Params_SolvePhaseEquilibrium params={this, 0, P, SOLVE_SATURATED_T};
            gsl_multiroot_function f = {&func_PhaseEquilibrium, n, &params};
            double x_init[n] = {T_K, Rho_l0, Rho_v0}; //initial guess
            gsl_vector *x = gsl_vector_alloc (n);
            gsl_vector_set (x, 0, x_init[0]);
            gsl_vector_set (x, 1, x_init[1]);
            gsl_vector_set (x, 2, x_init[2]);
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc (T, n);
            gsl_multiroot_fsolver_set (s, &f, x);
            // print_state_PhaseEquilibrium (iter, s);
            // iteration: solve
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);
                // print_state_PhaseEquilibrium (iter, s);
                if (status) /* check if solver is stuck */
                    break;
                status = gsl_multiroot_test_residual (s->f, TOL_PTRro);
            }while (status == GSL_CONTINUE && iter < ITERATION_MAX);
            if(!((status == GSL_SUCCESS))) // || (status == GSL_ENOPROG)
            {
                if(status == GSL_ENOPROG)
                {
                    print_state_TP2Rho (iter, s);
                    printf("P = %.3E Pa\n", P);
                    WARNING("Boiling_T: " + std::string(gsl_strerror (status)));
                }else
                {
                    print_state_TP2Rho (iter, s);
                    printf ("status = %s\n\n", gsl_strerror (status));
                    printf("P = %.3E Pa\n", P);
                    ERROR("Fatal error in Boiling_T(const double P, double& T_K, double& rho_l, double& rho_v)");
                }
            }
            T_K = gsl_vector_get (s->x, 0);
            rho_l = gsl_vector_get (s->x, 1);
            rho_v = gsl_vector_get (s->x, 2);
            // release pointer
            gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
        }

        double cIAPWS95::P_delta_tau(const double& delta, const double& tau)
        {
            return delta*m_constants.rhomass_critical*m_constants.R*m_constants.T_critical/tau * (1 + delta * phi_r_d(delta, tau));
        }

        /**
         * @brief Construct nonlinear function of \f$ f(\rho, p, T) \f$, which is used to solve density by given Temperautre and pressure.
         *
         * @param x
         * @param params
         * @param f
         * @return int
         */
        int func_TP2Rho(const gsl_vector * x, void *params, gsl_vector * f)
        {
            Params_TP2Rho* param = (Params_TP2Rho*)params;
            cIAPWS95* iapws = param->iapws;
            double T_K = param->T_K;
            double tau = param->tau;
            double p   = param->p;
            double RhocRT = param->RhocRT;
            const double rho    = gsl_vector_get (x, 0);
            const double delta  = rho/iapws->m_constants.rhomass_critical;
            const double y0     = RhocRT*delta*(1.0 + delta * iapws->phi_r_d(delta, tau)) - p;
            gsl_vector_set (f, 0, y0);
            return GSL_SUCCESS;
        }
        double func_TP2Rho(double rho, void *params)
        {
            Params_TP2Rho* param = (Params_TP2Rho*)params;
            cIAPWS95* iapws = param->iapws;
            const double delta  = rho/iapws->m_constants.rhomass_critical;
            return param->RhocRT*delta*(1.0 + delta * iapws->phi_r_d(delta, param->tau)) - param->p;
        }

        /**
         * @brief Solve nonlinear equation of \f$ p(\delta, \tau(T_0)) = \rho R T_0 (1 + \delta \phi^r_{\delta}) = P_0 \f$ by given \f$(P_0, T_0) \f$ to get density in single phase region. See Eq. 1 in Table 3 of \cite IAPWS-95.
         * \warning Has some issue close to phase boundary. But #Rho_bisection is much more reliable.
         * @param T_K
         * @param P
         * @return double
         */
        double cIAPWS95::Rho_Newton(const double T_K, const double P)
        {
            double Rho0 = m_constants.rhomass_critical;
            if (T_K<m_constants.T_critical) //if T<T_critical, use saturation density as initial guess (For T_min<T<T_c, the rhol_sat, rhov_sat alway exist); otherwise, use the critical density as initial guess
            {
                double p_sat, rhol_sat, rhov_sat;
                Boiling_p(T_K, p_sat, rhol_sat, rhov_sat);
                // printf("T: %f C, sat p: %f bar, rhol_sat: %f kg/m3, rhov_sat: %f kg/m3\n",T_K-273.15, P/1E5, rhol_sat, rhov_sat);
                if(P>p_sat)
                {
                    Rho0 = rhol_sat;
                }else
                {
                    Rho0 = rhov_sat;
                }
            }
            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            int status;
            size_t i, iter = 0;
            const size_t n = 1;
            Params_TP2Rho params={this, T_K, m_constants.T_critical/T_K, P, m_constants.rhomass_critical*m_constants.R*T_K};
            gsl_multiroot_function f = {&func_TP2Rho, n, &params};
            double x_init[n] = {Rho0}; //initial guess
            gsl_vector *x = gsl_vector_alloc (n);
            gsl_vector_set (x, 0, x_init[0]);
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc (T, n);
            gsl_multiroot_fsolver_set (s, &f, x);
            // print_state_TP2Rho (iter, s);
            // iteration: solve
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);
                // print_state_TP2Rho (iter, s);
                if (status) /* check if solver is stuck */
                    break;
                status = gsl_multiroot_test_residual (s->f, TOL_Pressure);
            }while (status == GSL_CONTINUE && iter < ITERATION_MAX);
            if(!((status == GSL_SUCCESS))) // || (status == GSL_ENOPROG)
            {
                if(status == GSL_ENOPROG)
                {
                    print_state_TP2Rho (iter, s);
                    printf("T_K = %.3f K, P = %.3E Pa\n", T_K, P);
                    WARNING("cIAPWS95::Rho : " + std::string(gsl_strerror (status)));
                }else
                {
                    print_state_TP2Rho (iter, s);
                    printf ("status = %s\n\n", gsl_strerror (status));
                    printf("T_K = %.3f K, P = %.3E Pa\n", T_K, P);
                    ERROR("Fatal error in cIAPWS95::Rho(const double T_K, const double P)");
                }
            }
            Rho0 = gsl_vector_get (s->x, 0);
            // release pointer
            gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
            return Rho0;
        }

        /**
         * @brief Use bisection method to solve one-dimension nonlinear equation for Rho(p,T).
         *
         * @param T_K
         * @param P
         * @return double
         */
        double cIAPWS95::Rho_bisection(const double T_K, const double P, double Rho_guess, double Rho_min, double Rho_max)
        {
            double rho = m_constants.rhomass_critical;
            double x_lo = 1E-4, x_hi = 1400;
            if (T_K<m_constants.T_critical) //if T<T_critical, use saturation density as initial guess (For T_min<T<T_c, the rhol_sat, rhov_sat alway exist); otherwise, use the critical density as initial guess
            {
                double p_sat, rhol_sat, rhov_sat;
                Boiling_p(T_K, p_sat, rhol_sat, rhov_sat);
                // printf("T: %f C, sat p: %f bar, rhol_sat: %f kg/m3, rhov_sat: %f kg/m3\n",T_K-273.15, P/1E5, rhol_sat, rhov_sat);
                if(P>p_sat)
                {
                    rho = rhol_sat;
                    x_lo = rhol_sat;
                }else
                {
                    rho = rhov_sat;
                    x_hi = rhov_sat;
                }
            }
//            else
//            {
//                Rho_lookup(rho, x_lo, x_hi, T_K, P);
//            }
            int status;
            int iter = 0;
            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;

            gsl_function F;
            Params_TP2Rho params={this, T_K, m_constants.T_critical/T_K, P, m_constants.rhomass_critical*m_constants.R*T_K};

            F.function = &func_TP2Rho;
            F.params = &params;

            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc (T);
            gsl_root_fsolver_set (s, &F, x_lo, x_hi);
            // printf ("using %s method\n", gsl_root_fsolver_name (s));
            // printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
//            printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,rho);
            do
            {
                iter++;
                status = gsl_root_fsolver_iterate (s);
                rho = gsl_root_fsolver_root (s);
                x_lo = gsl_root_fsolver_x_lower (s);
                x_hi = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (x_lo, x_hi, 0, TOL_Pressure/100);
                // if (status == GSL_SUCCESS) printf ("Converged:\n");
//                 printf ("%5d [%.7f, %.7f] %.7f\n",iter, x_lo, x_hi,rho);
            }
            while (status == GSL_CONTINUE && iter < ITERATION_MAX);
//            printf ("%5d [%.7f, %.7f] %.7f  [%.7f, %.7f]\n",iter, T_K, P, rho,x_lo,x_hi);
            if(!((status == GSL_SUCCESS)))
            {
                printf ("status = %s\n\n", gsl_strerror (status));
                printf("T_K = %.3f K, P = %.3E Pa\n", T_K, P);
                ERROR("Fatal error in cIAPWS95::Rho(const double T_K, const double P)");
            }
            gsl_root_fsolver_free (s);

            return rho;
        }
        /**
         * @brief Calculate density \f$ \rho \f$ by given temperature and pressure.
         *
         * \image html comparOtherCodes/TP/Rho.svg "Compare result with other codes." width=100%.
         * \note The python package **pyIAPWS95** and **pyIAPWS97** can be found at https://github.com/jjgomera/iapws, which implements both IAPWS95 and IAPWS-IF97 formula. **freesteam**(http://freesteam.sourceforge.net) and **PROST** (http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas) implement IAPWS-IF97 and IAPWS84, respectively.
         * @param T_K [K]
         * @param P [Pa]
         * @param method
         * @return double
         */
        double cIAPWS95::Rho(const double& T_K, const double& P, std::string method)
        {
            // printf("method: %s\n", method.c_str());
            if(method == "newton")
            {
                return Rho_Newton(T_K, P);
            }else if(method=="bisection")
            {
                return Rho_bisection(T_K, P);
            }else
            {
                WARNING("The method for Rho calculation only support [newton, bisection]: "+method+" is not supported, use default bisection.");
                return Rho_bisection(T_K, P);
            }
        }

        void cIAPWS95::_enthalpy(const double& T_K, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir, double& h)
        {
            h = (1.0 + tau*(phio.t + phir.t) + delta*phir.d)*m_constants.R*T_K;
        }

        /**
         * @brief Calculate specific enthalpy by given temperature and pressure.
         *
         * \image html comparOtherCodes/TP/H.svg "Compare result with other codes." width=100%.
         * \note The python package **pyIAPWS95** and **pyIAPWS97** can be found at https://github.com/jjgomera/iapws, which implements both IAPWS95 and IAPWS-IF97 formula. **freesteam**(http://freesteam.sourceforge.net) and **PROST** (http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas) implement IAPWS-IF97 and IAPWS84, respectively.
         * @param T_K
         * @param P
         * @param method
         * @return double
         */
        double cIAPWS95::_enthalpy(const double& T_K, const double& P, std::string method)
        {
            double rho = Rho(T_K, P, method);
            HelmholtzEnergy_dimensionless phio;
            HelmholtzEnergy_dimensionless phir;
            double delta = rho/m_constants.rhomass_critical;
            double tau   = m_constants.T_critical/T_K;
            phir.t = phi_r_t(delta, tau);
            phir.d = phi_r_d(delta, tau);
            phi_o(delta, tau, phio, Update_phi_t);
            // return (1.0 + tau*(phio.t + phir_t) + delta*phir_d)*m_constants.R*T_K;
            double h;
            _enthalpy(T_K, delta, tau, phio, phir, h);
            return h;
        }

        /**
         * @brief Adopt the simplified form (\f$ c_p \f$) in Table 3 of \cite IAPWS-95.
         * \note This equation is equivalent to \f$ \frac{\partial h}{\partial T}_p \f$ in Table 5 of \cite PartialDer
         * @param prop
         * @param rho
         * @param T
         * @param delta
         * @param tau
         * @param phio
         * @param phir
         */
        void cIAPWS95::_dhdT_P(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir)
        {
            // double dhdT_rho, dhdrho_T, dpdT_rho, dpdrho_T;
            // dhdT_rho = R*(-tau*tau * (phio.tt + phir.tt) + (1 + delta*phir.d - tau*delta*phir.dt));
            // dhdrho_T = T*R/rho * (tau*delta*phir.dt + delta*phir.d + delta*delta*phir.dd);
            // dpdT_rho = rho*R*(1.0 + delta*phir.d - tau*delta*phir.dt);
            // dpdrho_T = T*R*(1.0 + 2.0*delta*phir.d + delta*delta*phir.dd);
            // prop = dhdT_rho - dhdrho_T*dpdT_rho/dpdrho_T;
            prop = R*(-tau*tau * (phio.tt + phir.tt) + pow(1+delta*phir.d - delta*tau*phir.dt, 2.0)/(1+2.0*delta*phir.d + delta*delta*phir.dd));
        }
        /**
         * @brief Implementation of \f$ p(\delta, \tau), h(\delta, \tau) \f$ used in SinglePhase_HP.
         *
         * @param x
         * @param params
         * @param f
         * @return int
         */
        int func_HP2RhoT(const gsl_vector * x, void *params, gsl_vector * f)
        {
            Params_HP2RhoT* param = (Params_HP2RhoT*)params;
            cIAPWS95* iapws = param->iapws;
            double y0, y1;
            double p = param->p;
            double h = param->h;
            const double rho    = gsl_vector_get (x, 0);
            const double delta  = rho/iapws->m_constants.rhomass_critical;
            double T_K        = gsl_vector_get (x, 1);
            double tau          = iapws->m_constants.T_critical/T_K;
            double RT           = R*T_K;
            HelmholtzEnergy_dimensionless phio;
            iapws->phi_o(delta, tau, phio, Update_phi_all);
            y0     = h/RT -1.0 - tau*(phio.t + iapws->phi_r_t(delta, tau)) - delta*iapws->phi_r_d(delta, tau);      /**< Equation \f$ h(\delta, \tau) \f$ */
            y1     = p/(RT*rho) - 1.0 - delta * iapws->phi_r_d(delta, tau);                                         /**< Equation \f$ p(\delta, \tau) \f$ */
            gsl_vector_set (f, 0, y0);
            gsl_vector_set (f, 1, y1);

            return GSL_SUCCESS;
        }
        void print_state_HP2RhoT(size_t iter, gsl_multiroot_fsolver * s)
        {
            printf ("iter = %3lu Rho = % 15.8f, T = % 15.8f, err_h(x) = % .3e, err_p(x) = % .3e \n",iter,
                    gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
                    gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
        }

        /**
         * @brief Solving temperature and density in single phase region by combining pressure equation \f$ p(\delta, \tau) \f$ and enthalpy equation \f$ h(\delta, \tau) \f$ in Table 3 of \cite IAPWS-95
         *
         * @param H Specific enthalpy [J/kg]
         * @param P Pressure [Pa]
         * @param rho Return density [kg/m3]
         * @param T_K Return temperature [K]. Have to give a initial value for T_K.
         * @param method
         */
        void cIAPWS95::SinglePhase_HP(const double& H, const double& P, double& rho, double& T_K, std::string method)
        {
            rho = Rho(T_K, P, method); //Initial value of rho

            const gsl_multiroot_fsolver_type *T;
            gsl_multiroot_fsolver *s;
            int status;
            size_t i, iter = 0;
            const size_t n = 2;
            Params_HP2RhoT params={this, H, P};
            gsl_multiroot_function f = {&func_HP2RhoT, n, &params};
            double x_init[n] = {rho, T_K}; //initial guess
            gsl_vector *x = gsl_vector_alloc (n);
            gsl_vector_set (x, 0, x_init[0]);
            gsl_vector_set (x, 1, x_init[1]);
            T = gsl_multiroot_fsolver_hybrids;
            s = gsl_multiroot_fsolver_alloc (T, n);
            gsl_multiroot_fsolver_set (s, &f, x);
            // print_state_HP2RhoT (iter, s);
            // iteration: solve
            do
            {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s);
                // print_state_HP2RhoT (iter, s);
                if (status) /* check if solver is stuck */
                    break;
                status = gsl_multiroot_test_residual (s->f, TOL_PTRro);
            }while (status == GSL_CONTINUE && iter < ITERATION_MAX);
            // printf ("status = %s\n", gsl_strerror (status));
            if(!((status == GSL_SUCCESS))) // || (status == GSL_ENOPROG)
            {
                if(status == GSL_ENOPROG)
                {
                    print_state_HP2RhoT (iter, s);
                    printf("H = %.3f J/kg, P = %.3f Pa\n", H, P);
                    WARNING("SinglePhase_HP: " + std::string(gsl_strerror (status)));
                }else
                {
                    print_state_HP2RhoT (iter, s);
                    printf ("status = %s\n\n", gsl_strerror (status));
                    printf("H = %.3f J/kg, P = %.3f Pa\n", H, P);
                    ERROR("Fatal error in cIAPWS95::SinglePhase_HP(const double& H, const double& P, double& rho, double& T_K, std::string method)");
                }
            }
            rho = gsl_vector_get (s->x, 0);
            T_K = gsl_vector_get (s->x, 1);
            // release pointer
            gsl_multiroot_fsolver_free (s); gsl_vector_free (x);
        }

        void cIAPWS95::UpdateState_TP(ThermodynamicProperties& props, State& state, const double& T, const double& P)
        {
            props.fluidName = name();

            if(T<m_constants.Tmin || T>m_constants.Tmax)
            {
                printf("T = %f K, P = %f Pa\n", T, P);
                throw OutOfRangeError("Fatal error in void cIAPWS95::UpdateState_TP(const double& T, const double& P, State& state)\nT out of bound: T["
                +std::to_string(m_constants.Tmin)+", "+std::to_string(m_constants.Tmax)+"], P["+std::to_string(m_constants.pmin)+", "+std::to_string(m_constants.pmax)+"]");
            }
            props.T = T;
            props.p = P;
            props.phase = findPhaseRegion_TPX(T,P);
            state.phase = props.phase;
            // calculate Rho
            //props.Rho=Rho(T, P);
            props.Rho = Rho_bisection(T, P);
            // calculate others based on Rho
            state.delta.singlePhase = props.Rho/m_constants.rhomass_critical;
            state.tau   = m_constants.T_critical/props.T;
            phi_o(state.delta.singlePhase, state.tau, state.phi.singlePhase.phio);
            phi_r(state.delta.singlePhase, state.tau, state.phi.singlePhase.phir);
            // 1. H
            _h(props.H,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
            // 2. Cp
            _dhdT_P(props.Cp,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
            // 3. Mu
            Viscosity_H2O_IAPWS2008(props.T, props.Rho, props.Mu);
            // 4. Derivatives
            // 4.1 dPdRho_T
            double tmp_dPdRho_T, tmp_dPdT_Rho;
            _dPdRho_T(tmp_dPdRho_T,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
            // 4.2 dPdT_Rho
            _dPdT_Rho(tmp_dPdT_Rho,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
            // ============== derivatives ==============
            props.IsothermalCompressibility = 1.0/(props.Rho * tmp_dPdRho_T);
            props.IsobaricExpansivity = tmp_dPdT_Rho/tmp_dPdRho_T/props.Rho; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
            props.dRhodP = 1.0/tmp_dPdRho_T;
            props.dRhodT = -tmp_dPdT_Rho * props.dRhodP;
            //==========================================
            //==========================================
            switch (props.phase) {
                case SinglePhase_V:
                    props.Rho_v = props.Rho;
                    props.H_v = props.H;
                    props.Mu_v = props.Mu;
                    props.Cp_v = props.Cp;
                    break;
                default:
                    props.Rho_l = props.Rho;
                    props.H_l = props.H;
                    props.Mu_l = props.Mu;
                    props.Cp_l = props.Cp;
                    break;
            }
        }

        /**
         * @brief Calculate state of H2O from IAPWS-95 formula by given specific enthalpy and pressure.
         *
         * @param H Specific enthalpy [J/kg]
         * @param P Pressure [Pa]
         * @param state Thermodynamic state, see #State for data struct definition. It contains the basic information for describing a state and calculating other properties and their derivatives.
         * @param method Method for nonlinear equation solving, default is "bisection".
         */
        void cIAPWS95::UpdateState_HP(ThermodynamicProperties& props, State& state, const double& H, const double& P, std::string method)
        {
            props.fluidName = name();

            props.H = H;
            props.p = P;
            // calculate H_min first, according to m_constants.Tmin and P
            double hmin = _enthalpy(m_constants.Tmin, P);
            double hmax = _enthalpy(m_constants.Tmax, P);
            if(H<hmin)props.H=hmin; //Force H = hmin if it is less than the minimum bound
            else if(H>hmax)props.H=hmax;

            // 1. check phase region and get initial guess of T
            if(P>=m_constants.p_critical && P<=m_constants.pmax) //MUST BE in the supercritical region
            {
                // STATUS("Pressure greater than critical pressure");
                props.phase = Supercritical;
                if(P<=CONST_IF97_Pmax_Region1) //in valid IF97 region
                {
                    int region = m_IF97.GetRegion_PH(P, props.H);
                    switch (region)
                    {
                        case IF97_REGION1:
                            props.T = m_IF97.Backward_T_PH_region1(P,props.H);
                            break;
                        case IF97_REGION3a:
                            props.T = m_IF97.Backward_T_PH_region3a(P, props.H);
                            break;
                        case IF97_REGION3b:
                            props.T = m_IF97.Backward_T_PH_region3b(P, props.H);
                            break;
                        case IF97_REGION2c:
                            props.T = m_IF97.Backward_T_PH_region2c(P, props.H);
                            break;
                        case IF97_REGION2b:
                            props.T = m_IF97.Backward_T_PH_region2b(P, props.H);
                            break;
                        case IF97_REGION5:
                            props.T = 1273.13; //For IF97 region 5, there is not an explicit backward equation to calculate T, so here use a random number in region 5 as a initial guess.
                            break;
                        case IF97_RegionUndefined:
                            props.T = 1273.13; //For the upper part (P>50MPa) of region 5, it is not defined in IF97, so here use a random number in region 5 as a initial guess.
                            break;
                        default:
                        ERROR("Unknown IF97 region index in supper critical region: "+std::to_string(region));
                            break;
                    }
                }else
                {
                    // For high pressure, H-T relation is close to a linear trend, so here use a polyfit equation of T(H)(P=100MPa) to calculate initial guess of T
                    double coeffs[6]={2.52879000e+02, 2.32953512e-04, 2.66741583e-11 -1.67239088e-17, 2.48755621e-25, 5.73883079e-31};
                    props.T = 0;
                    for (int i = 0; i < 6; i++)
                    {
                        props.T += coeffs[i]*pow(props.H, i);
                    }
                }
            }else if(P<m_constants.p_critical && P>=m_constants.pmin)
            {
                double T_sat;
                // STATUS("Boiling_T");
                // in liquid-phase region, two phase region, or vapor-phase region
                Boiling_T(P, T_sat, props.Rho_l, props.Rho_v);
                // STATUS("Boiling_T solved")
                // update phio, phir
                state.tau = m_constants.T_critical/T_sat;
                state.delta.twoPhase.delta_l = props.Rho_l/m_constants.rhomass_critical;
                state.delta.twoPhase.delta_v = props.Rho_v/m_constants.rhomass_critical;
                phi_o(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, Update_phi_all);
                phi_o(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, Update_phi_all);
                phi_r(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phir_l);
                phi_r(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phir_v);
                _enthalpy(T_sat, state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l, props.H_l);
                _enthalpy(T_sat, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v, props.H_v);

                // STATUS("Saturated temperature and rhol, rhov are solved");
                // printf("T_sat = %f K\n", T_sat);
                // check phase region
                // Note that the initial guess of T depends on IF97 region, but the region is not consistent with IAPWS95, there are slight differences between two equations. e.g., for a H<state.liquid.H.value , the region index given by IF97 would be two phase region(region 4), for this case we use a small shift of T_sat as the initial guess
                if(props.H<props.H_l)               //liquid phase region
                {
                    props.phase = SinglePhase_L;
                    int region = m_IF97.GetRegion_PH(P, props.H);
                    // STATUS("IF97 phase region: "+std::to_string(region));
                    switch (region)
                    {
                        case IF97_REGION1:
                            // STATUS("initial T: "+std::to_string(Backward_T_PH_region1(P, state.H.value)));
                            props.T = m_IF97.Backward_T_PH_region1(P, props.H);
                            break;
                        case IF97_REGION3a:
                            props.T = m_IF97.Backward_T_PH_region3a(P, props.H);
                            break;
                        case IF97_REGION4:
                            props.T = T_sat - 0.1; //Relax: Because the small difference on phase boundary between IF97 and IAPWS95 equations, here use a little bit smaller value than T_sat as the initial guess.
                            break;
                        default:
                        ERROR("Unknown IF97 phase region in pure liquid region in cIAPWS95::UpdateState_HP: "+std::to_string(region));
                            break;
                    }
                    if(props.T > T_sat)props.T = T_sat - 0.01; //Relax: Make sure it is indeed in the liquid phase region
                }else if(props.H>props.H_v)         //vapor phase region
                {
                    props.phase = SinglePhase_V;
                    int region = m_IF97.GetRegion_PH(P, props.H);
                    switch (region)
                    {
                        case IF97_REGION2a:
                            props.T = m_IF97.Backward_T_PH_region2a(P, props.H);
                            break;
                        case IF97_REGION2b:
                            props.T = m_IF97.Backward_T_PH_region2b(P, props.H);
                            break;
                        case IF97_REGION2c:
                            props.T = m_IF97.Backward_T_PH_region2c(P, props.H);
                            break;
                        case IF97_REGION3b:
                            props.T = m_IF97.Backward_T_PH_region3b(P, props.H);
                            break;
                        case IF97_REGION5:    //For IF97 region 5, there is not an explicit backward equation to calculate T, so here use a random number in region 5 as a initial guess.
                            props.T = 1273.15;
                            break;
                        case IF97_REGION4:
                            props.T = T_sat + 0.01; //Relax
                            break;
                        default:
                        ERROR("Unknown IF97 phase region in pure vapor region in cIAPWS95::UpdateState_HP: "+std::to_string(region));
                            break;
                    }
                    if(props.T < T_sat)props.T = T_sat + 0.01; //Relax: Make sure it is indeed in the vapor phase region
                    if(props.T>T_critical())props.phase = Supercritical_vapor;
                }else                               //Two phase region
                {
                    // STATUS("Enthalpy between h_l, h_v");
                    props.phase = TwoPhase_VL_Water;
                    props.T = T_sat;
                    state.x = (props.H - props.H_l)/(props.H_v - props.H_l); //vapor mass fraction
                    props.Rho = 1.0/(state.x/props.Rho_v + (1.0 - state.x)/props.Rho_l);
                }
            }
            else
            {
                printf("H = %f J/kg, P = %f Pa\n", props.H, P);
                ERROR("Fatal error in cIAPWS95::UpdateState_HP(const double H, const double P, State& state, std::string method)\nP out of bound");
            }

            // calculate single phase T and Helmholtz free energy and properties based on Rho
            if(props.phase!=TwoPhase_VL_Water)
            {
                SinglePhase_HP(props.H, P, props.Rho, props.T, method);
                state.delta.singlePhase = props.Rho/m_constants.rhomass_critical;
                state.tau   = m_constants.T_critical/props.T;
                phi_o(state.delta.singlePhase, state.tau, state.phi.singlePhase.phio);
                phi_r(state.delta.singlePhase, state.tau, state.phi.singlePhase.phir);
                //Mu
                Viscosity_H2O_IAPWS2008(props.T, props.Rho, props.Mu);
                // ============== derivatives ==============
                //Cp
                _dhdT_P(props.Cp,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                // dPdRho_P
                double tmp_dPdRho_T, tmp_dPdT_Rho;
                _dPdRho_T(tmp_dPdRho_T,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                // dPdT_Rho
                _dPdT_Rho(tmp_dPdT_Rho,props.Rho,props.T,state.delta.singlePhase,state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                props.IsothermalCompressibility = 1.0/(props.Rho * tmp_dPdRho_T);
                props.IsobaricExpansivity = tmp_dPdT_Rho/tmp_dPdRho_T/props.Rho; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
                props.dRhodP = 1.0/tmp_dPdRho_T;
                props.dRhodT = -tmp_dPdT_Rho * props.dRhodP;
                //==========================================
                switch (props.phase) {
                    case SinglePhase_V:
                        props.Rho_v = props.Rho;
                        props.H_v = props.H;
                        props.Mu_v = props.Mu;
                        props.Cp_v = props.Cp;
                        break;
                    default:
                        props.Rho_l = props.Rho;
                        props.H_l = props.H;
                        props.Mu_l = props.Mu;
                        props.Cp_l = props.Cp;
                        break;
                }
            } else
            {
                HelmholtzEnergy_dimensionless phio, phir;
                double delta = props.Rho/m_constants.rhomass_critical;
                phi_o(delta, state.tau, phio);
                phi_r(delta, state.tau, phir);
                // ============== derivatives ==============
                // Some derivatives to T is not defined in the two phase region, use their initial value 0

                //Cp: dHdT_P
                // _dhdT_P(props.Cp_l,props.Rho_l,props.T,state.delta.twoPhase.delta_l,state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
                // _dhdT_P(props.Cp_v,props.Rho_v,props.T,state.delta.twoPhase.delta_v,state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);
                // _dhdT_P(props.Cp, props.Rho, props.T, delta, state.tau, phio, phir);
                // // dPdRho_P
                // double tmp_dPdRho_T, tmp_dPdRho_T_l, tmp_dPdRho_T_v, tmp_dPdT_Rho, tmp_dPdT_Rho_l, tmp_dPdT_Rho_v;
                // _dPdRho_T(tmp_dPdRho_T_l,props.Rho_l,props.T,state.delta.twoPhase.delta_l,state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
                // _dPdRho_T(tmp_dPdRho_T_v,props.Rho_v,props.T,state.delta.twoPhase.delta_v,state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);
                // _dPdRho_T(tmp_dPdRho_T, props.Rho, props.T, delta, state.tau, phio, phir);
                // // dPdT_Rho
                // _dPdT_Rho(tmp_dPdT_Rho_l,props.Rho_l,props.T,state.delta.twoPhase.delta_l,state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
                // _dPdT_Rho(tmp_dPdT_Rho_v,props.Rho_v,props.T,state.delta.twoPhase.delta_v,state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);
                // _dPdT_Rho(tmp_dPdT_Rho, props.Rho, props.T, delta, state.tau, phio, phir);
                // // commonly used derivatives
                // props.IsothermalCompressibility_l = 1.0/(props.Rho_l * tmp_dPdRho_T_l);
                // props.IsothermalCompressibility_v = 1.0/(props.Rho_v * tmp_dPdRho_T_v);
                // props.IsothermalCompressibility = 1.0/(props.Rho * tmp_dPdRho_T);
                // props.IsobaricExpansivity_l = tmp_dPdT_Rho_l/tmp_dPdRho_T_l/props.Rho_l;
                // props.IsobaricExpansivity_v = tmp_dPdT_Rho_v/tmp_dPdRho_T_v/props.Rho_v;
                // props.IsobaricExpansivity = tmp_dPdT_Rho/tmp_dPdRho_T/props.Rho; // See Table 5 of Thorade et al.(2013, doi:10.1007/s12665-013-2394-z)
                // props.dRhodP_l = 1.0/tmp_dPdRho_T_l;
                // props.dRhodT_l = -tmp_dPdT_Rho_l * props.dRhodP_l;
                // props.dRhodP_v = 1.0/tmp_dPdRho_T_v;
                // props.dRhodT_v = -tmp_dPdT_Rho_v * props.dRhodP_v;
                // props.dRhodP = 1.0/tmp_dPdRho_T;
                // props.dRhodT = -tmp_dPdT_Rho * props.dRhodP;
                //==========================================
                // Mu
                Viscosity_H2O_IAPWS2008(props.T, props.Rho_l, props.Mu_l);
                Viscosity_H2O_IAPWS2008(props.T, props.Rho_v, props.Mu_v);
                // saturation
                props.S_v = state.x;
                props.S_l = 1.0 - state.x;
                //Mu
                Viscosity_H2O_IAPWS2008(props.T, props.Rho, props.Mu);
            }

            // STATUS("phase index: "+std::to_string(state.phase));
        }

        void cIAPWS95::_h(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir)
        {
            prop = (1.0 + tau*(phio.t + phir.t) + delta*phir.d)*R*T;
        }

        void cIAPWS95::h(ThermodynamicProperties& props,const State& state)
        {
            switch (props.phase)
            {
                case SinglePhase_L:
                {
                    _h(props.H_l, props.Rho_l, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.H_v = 0;
                    props.H = props.H_l;
                }
                    break;
                case SinglePhase_V:
                {
                    _h(props.H_v, props.Rho_v, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.H_l = 0;
                    props.H = props.H_v;
                }
                    break;
                case Supercritical:    //Regard super critical fluid as liquid.
                {
                    _h(props.H_l, props.Rho_l, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.H_v = 0;
                    props.H = props.H_l;
                }
                    break;
                case TwoPhase_VL_Water:
                {
                    _h(props.H_l, props.Rho_l, props.T, state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
                    _h(props.H_v, props.Rho_v, props.T, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);
                    props.H = (1.0 - state.x)*props.H_l + props.H_v*state.x;
                }
                    break;
                default:
                ERROR("Unknown Phase index of IAPWS95 in void cIAPWS95::h(double& h, double& prop_l, double& prop_v): "+std::to_string(props.phase));
                    break;
            }
        }

        void cIAPWS95::dhdT_P(ThermodynamicProperties& props,const State& state)
        {
            switch (props.phase)
            {
                case SinglePhase_L:
                {
                    _dhdT_P(props.Cp_l, props.Rho_l, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.Cp_v = 0;
                    props.Cp = props.Cp_l;
                }
                    break;
                case SinglePhase_V:
                {
                    _dhdT_P(props.Cp_v, props.Rho_v, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.Cp_l = 0;
                    props.Cp = props.Cp_v;
                }
                    break;
                case Supercritical:    //Regard super critical fluid as liquid.
                {
                    _dhdT_P(props.Cp_l, props.Rho_l, props.T, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
                    props.Cp_v = 0;
                    props.Cp = props.Cp_l;
                }
                    break;
                case TwoPhase_VL_Water:
                {
                    _dhdT_P(props.Cp_l, props.Rho_l, props.T, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
                    _dhdT_P(props.Cp_v, props.Rho_v, props.T, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);
                    props.Cp = (1.0 - state.x)*props.Cp_l + props.Cp_v*state.x;
                }
                    break;
                default:
                ERROR("Unknown Phase index of IAPWS95 in void cIAPWS95::_dhdT_P(double& h, double& prop_l, double& prop_v): "+std::to_string(props.phase));
                    break;
            }
        }

        double cIAPWS95::Boiling_p(const double &T) {
            double P, rhol, rhov;
            Boiling_p(T, P, rhol, rhov);
            return P;
        }

        double cIAPWS95::Boiling_p(const double &T, double &rho_l, double &rho_v) {
            double P;
            Boiling_p(T, P, rho_l, rho_v);
            return P;
        }

        void
        cIAPWS95::UpdateState_TPX(ThermodynamicProperties &props, const double &T, const double &P, const double &X) {
            State state;
            UpdateState_TP(props, state, T, P);
        }

        PhaseRegion cIAPWS95::findPhaseRegion_TPX(const double &T, const double &P, const double &X) {
            if((P>=m_constants.p_critical && P<=m_constants.pmax) || (T>=m_constants.T_critical && T<=m_constants.Tmax))
            {
                if(T<=m_constants.T_critical)
                {
                    return Supercritical_liquid;
                }else if(P<=m_constants.p_critical)
                {
                    return Supercritical_vapor;
                }else
                {
                    return Supercritical;
                }
            }else if(P<m_constants.p_critical && P>=m_constants.pmin)
            {
                double T_sat = m_IF97.T_sat_P(P);
                if(T< (T_sat - 5))
                {
                    return SinglePhase_L;
                }else if(T> (T_sat+5))
                {
                    return SinglePhase_V;
                }else
                {
                    double rhol, rhov;
                    Boiling_T(P, T_sat, rhol, rhov);
                    if(T<=T_sat)
                    {
                        return SinglePhase_L;
                    }else
                    {
                        return SinglePhase_V;
                    }
                }
            }else
            {
                printf("T = %f K, P = %f Pa\n", T, P);
                throw OutOfRangeError ("Fatal error in void cIAPWS95::findPhaseRegion_TPX(const double &T, const double &P, const double &X)\nP out of bound");
            }
        }

        double cIAPWS95::Boiling_p(const double &T, ThermodynamicProperties &props) {
            props.fluidName = name();
            props.T = T;
            Boiling_p(T, props.p, props.Rho_l, props.Rho_v);
            // calculate H_l, H_v
            State state;
            state.tau = m_constants.T_critical/T;
            state.delta.twoPhase.delta_l = props.Rho_l/m_constants.rhomass_critical;
            state.delta.twoPhase.delta_v = props.Rho_v/m_constants.rhomass_critical;
            phi_o(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, Update_phi_all);
            phi_o(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, Update_phi_all);
            phi_r(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phir_l);
            phi_r(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phir_v);
            _enthalpy(T, state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l, props.H_l);
            _enthalpy(T, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v, props.H_v);
            _dhdT_P(props.Cp_l, props.Rho_l,props.T,state.delta.twoPhase.delta_l,state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l);
            _dhdT_P(props.Cp_v, props.Rho_v,props.T,state.delta.twoPhase.delta_v,state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v);

            return props.p;
        }

        double cIAPWS95::Boiling_T(const double &p) {
            double T_sat, rhol, rhov;
            Boiling_T(p, T_sat, rhol, rhov);
            return T_sat;
        }

        double cIAPWS95::Boiling_T(const double &p, double &rho_l, double &rho_v) {
            double T_sat;
            Boiling_T(p, T_sat, rho_l, rho_v);
            return T_sat;
        }

        double cIAPWS95::Boiling_T(const double &p, ThermodynamicProperties &props) {
            props.p = p;
            Boiling_T(p, props.T, props.Rho_l, props.Rho_v);
            // calculate H_l, H_v
            State state;
            state.tau = m_constants.T_critical/props.T;
            state.delta.twoPhase.delta_l = props.Rho_l/m_constants.rhomass_critical;
            state.delta.twoPhase.delta_v = props.Rho_v/m_constants.rhomass_critical;
            phi_o(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, Update_phi_all);
            phi_o(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, Update_phi_all);
            phi_r(state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phir_l);
            phi_r(state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phir_v);
            _enthalpy(props.T, state.delta.twoPhase.delta_l, state.tau, state.phi.twoPhase.phio_l, state.phi.twoPhase.phir_l, props.H_l);
            _enthalpy(props.T, state.delta.twoPhase.delta_v, state.tau, state.phi.twoPhase.phio_v, state.phi.twoPhase.phir_v, props.H_v);

            return props.p;
        }

        void
        cIAPWS95::UpdateState_HPX(ThermodynamicProperties &props, const double &H, const double &p, const double &X) {
            State state;
            UpdateState_HP(props, state, H, p);
        }

        PhaseRegion cIAPWS95::findPhaseRegion_HPX(const double &H, const double &p, const double &X) {
            ThermodynamicProperties props;
            props.H = H;
            props.p = p;
            // calculate H_min first, according to m_constants.Tmin and P
            double hmin = _enthalpy(m_constants.Tmin, p);
            double hmax = _enthalpy(m_constants.Tmax, p);
            if(H<hmin)props.H=hmin; //Force H = hmin if it is less than the minimum bound
            else if(H>hmax)props.H=hmax;

            //
            if(p>=m_constants.p_critical && p<=m_constants.pmax) //MUST BE in the supercritical region
            {
                double h_crit = _enthalpy(m_constants.T_critical, p);
                if(H>h_crit)
                {
                    return  Supercritical;
                } else
                {
                    return Supercritical_liquid;
                }
            }else if(p<m_constants.p_critical && p>=m_constants.pmin)
            {
                Boiling_T(p, props);
                if(H<props.H_l)               //liquid phase region
                {
                    return SinglePhase_L;

                }else if(props.H>props.H_v)         //vapor phase region
                {
                    double h_crit = _enthalpy(m_constants.T_critical, p);
                    if(H>h_crit)
                    {
                        return Supercritical_vapor;
                    } else
                    {
                        return SinglePhase_V;
                    }
                }else                               //Two phase region
                {
                    return TwoPhase_VL_Water;
                }
            }
            else
            {
                printf("H = %f J/kg, P = %f Pa\n", props.H, p);
                ERROR("Fatal error in cIAPWS95::UpdateState_HP(const double H, const double P, State& state, std::string method)\nP out of bound");
            }
        }

        /**
        * @brief Calculate dynamic viscosity of H2O based on IAPWS release 2008, see \cite mu2008.
        * @param T [K]
        * @param Rho [kg/m^3]
        * @param dRhodP [SI]
        * @param Mu [Pa s]
        */
        void cIAPWS95::Viscosity_H2O_IAPWS2008(const double &T, const double &Rho, double &Mu) {
            double Tbar = T / m_constants_mu2008.T_star; //eq.(5)
            double Rhobar = Rho / m_constants_mu2008.rho_star; //eq.(6)
            // calculate mu0, eq(11)
            double OnebyTbar_Mu0[m_constants_mu2008.N_Hi];
            OnebyTbar_Mu0[0] = 1;
            OnebyTbar_Mu0[1] = 1.0/Tbar;
            OnebyTbar_Mu0[2] = OnebyTbar_Mu0[1] * OnebyTbar_Mu0[1];
            OnebyTbar_Mu0[3] = OnebyTbar_Mu0[2] * OnebyTbar_Mu0[1];
            double H_by_Tbar = m_constants_mu2008.H_i[0] * OnebyTbar_Mu0[0]
                               + m_constants_mu2008.H_i[1] * OnebyTbar_Mu0[1]
                               + m_constants_mu2008.H_i[2] * OnebyTbar_Mu0[2]
                               + m_constants_mu2008.H_i[3] * OnebyTbar_Mu0[3];
            double Mu0bar = 100.0 * sqrt(Tbar) / H_by_Tbar;
            // calculate mu1, eq(12)
            double Tbar_Mu1[m_constants_mu2008.row_Hij];
            Tbar_Mu1[0] = 1.0;
            Tbar_Mu1[1] = (OnebyTbar_Mu0[1] - 1.0);
            Tbar_Mu1[2] = Tbar_Mu1[1] * Tbar_Mu1[1];
            Tbar_Mu1[3] = Tbar_Mu1[2] * Tbar_Mu1[1];
            Tbar_Mu1[4] = Tbar_Mu1[3] * Tbar_Mu1[1];
            Tbar_Mu1[5] = Tbar_Mu1[4] * Tbar_Mu1[1];
            double Rhobar_Mu1[m_constants_mu2008.col_Hij];
            Rhobar_Mu1[0] = 1.0;
            Rhobar_Mu1[1] = Rhobar - 1.0;
            Rhobar_Mu1[2] = Rhobar_Mu1[1] * Rhobar_Mu1[1];
            Rhobar_Mu1[3] = Rhobar_Mu1[2] * Rhobar_Mu1[1];
            Rhobar_Mu1[4] = Rhobar_Mu1[3] * Rhobar_Mu1[1];
            Rhobar_Mu1[5] = Rhobar_Mu1[4] * Rhobar_Mu1[1];
            Rhobar_Mu1[6] = Rhobar_Mu1[5] * Rhobar_Mu1[1];
            double sum_Mu1 = 0;
            for (int i = 0; i < m_constants_mu2008.row_Hij; ++i) {
                for (int j = 0; j < m_constants_mu2008.col_Hij; ++j) {
                    sum_Mu1 += Tbar_Mu1[i] * m_constants_mu2008.H_ij[i][j] * Rhobar_Mu1[j];
                }
            }
            double Mu1bar = exp(Rhobar * sum_Mu1);

            // calculate Mu2, equation (14, 15, )
            // -- calculate dRhodP_T
            double dRhodP = 0, dRhodP_TR = 0;
            // 1. calculate dRhodP_T for T = T, rho = Rho
            State state;
            state.delta.singlePhase = Rho/m_constants.rhomass_critical;
            state.tau   = m_constants.T_critical/T;
            phi_o(state.delta.singlePhase, state.tau, state.phi.singlePhase.phio);
            phi_r(state.delta.singlePhase, state.tau, state.phi.singlePhase.phir);
            _dPdRho_T(dRhodP, Rho, m_constants_mu2008.T_R, state.delta.singlePhase, state.tau, state.phi.singlePhase.phio, state.phi.singlePhase.phir);
            dRhodP = 1.0/dRhodP;
            // 2. calculate dRhodP_T for T = T_R, rho = Rho
            State state_TR;
            state_TR.delta.singlePhase = Rho/m_constants.rhomass_critical;
            state_TR.tau   = m_constants.T_critical/m_constants_mu2008.T_R;
            phi_o(state_TR.delta.singlePhase, state_TR.tau, state_TR.phi.singlePhase.phio);
            phi_r(state_TR.delta.singlePhase, state_TR.tau, state_TR.phi.singlePhase.phir);
            _dPdRho_T(dRhodP_TR, Rho, m_constants_mu2008.T_R, state_TR.delta.singlePhase, state_TR.tau, state_TR.phi.singlePhase.phio, state_TR.phi.singlePhase.phir);
            dRhodP_TR = 1.0/dRhodP_TR;
            // -- calculate Y
            double Delta_Xbar = std::max(0.0, Rhobar * m_constants_mu2008.p_starByRho_star * (dRhodP - dRhodP_TR * m_constants_mu2008.T_R_bar/Tbar)); // see comment below equation(21a)
            double xi = m_constants_mu2008.xi0 * pow(Delta_Xbar/m_constants_mu2008.Gamma0, m_constants_mu2008.nuBygamma);
            double Y = 0;
            double qC_times_xi = m_constants_mu2008.qC * xi;
            double qC_times_xi_sqr = qC_times_xi * qC_times_xi;
            double qD_times_xi = m_constants_mu2008.qD * xi;
            if (xi <= 0.3817016416)
            {
               Y = 0.2 * qC_times_xi * pow(qD_times_xi, 5.0) * (1.0 - qC_times_xi + qD_times_xi*qD_times_xi - 765.0/504.0*qD_times_xi*qD_times_xi);
            } else
            {
                double psiD = acos(1.0/sqrt(1 + qD_times_xi*qD_times_xi));
                double w = sqrt(fabs((qC_times_xi - 1.0)/(qC_times_xi + 1.0))) * tan(psiD/2.0);
                double Lw = 0;
                if(qC_times_xi>1)
                {
                    Lw = log((1.0 + w)/(1.0 - w));
                } else
                {
                    Lw = 2.0* atan(fabs(w));
                }
                Y = 1.0/12.0 * sin(3.0 * psiD) - 0.25/qC_times_xi*sin(2.0 * psiD) + 1.0/qC_times_xi_sqr*(1.0 - 1.25*qC_times_xi_sqr)*sin(psiD) - 1.0/pow(qC_times_xi,3.0)*((1.0 - 1.5*qC_times_xi_sqr)*psiD - pow(fabs(qC_times_xi_sqr - 1.0), 1.5) * Lw);
            }
            double Mu2bar = exp(m_constants_mu2008.x_mu * Y);
            Mu = Mu0bar * Mu1bar * Mu2bar * m_constants_mu2008.mu_star;
        }

        void cIAPWS95::_dPdRho_T(double &prop, const double& rho, const double& T, const double &delta, const double &tau, const HelmholtzEnergy_dimensionless &phio, const HelmholtzEnergy_dimensionless &phir) const
        {
            prop = T * m_constants.R * (1.0 + 2.0 * delta * phir.d + delta*delta*phir.dd);
        }

        void cIAPWS95::_dPdT_Rho(double& prop, const double& rho, const double& T, const double& delta, const double& tau, const HelmholtzEnergy_dimensionless& phio, const HelmholtzEnergy_dimensionless& phir) const
        {
            prop = rho * m_constants.R * (1.0 + delta * phir.d - tau * delta * phir.dt);
        }

        /**
        * @brief Calculate dynamic viscosity for give T,p
        * @param T [K]
        * @param P [Pa]
        * @return [Pa s]
        */
        double cIAPWS95::Mu(const double &T, const double &P) {
            double Rho = Rho_bisection(T, P);
            double mu = 0;
            Viscosity_H2O_IAPWS2008(T, Rho, mu);
            return mu;
        }

        // void cIAPWS95::UpdateState_TPX(ThermodynamicPropertiesArray &stateArray, const size_t &N, const double *T, const double *p, const double *X)
        // {
        //     ThermodynamicProperties props;
        //     stateArray.N = N;
        //     for (size_t i = 0; i < N; ++i)
        //     {
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

        /**
        * @brief Computer-program verification given by \cite mu2008.
        * @return
        */
        void cIAPWS95::Verification_Mu() {
            std::vector<double> T = {298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15};
            std::vector<double> rho = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};
            std::vector<double> mu0 = {889.7351, 1437.649467, 307.883622, 14.538324, 217.685358, 32.619287, 35.802262, 77.430175, 44.217245, 47.640433, 64.154608};
            double mu;
            printf("Verification of xThermal implementation of dynamic viscosity of H2O of IAPWS 2008 release\n");
            for (int i = 0; i < T.size(); ++i) {
                Viscosity_H2O_IAPWS2008(T[i], rho[i], mu);
                printf("T = %7.2f, rho = %5.0f, mu_verification = %12.6f uPa s, mu_xThermal: %12.6f uPa s, err = %10.6f uPa s\n", T[i], rho[i], mu0[i], mu*1E6, mu0[i]-mu*1E6);
            }
        }

        ThermodynamicProperties cIAPWS95::UpdateState_TPX(const double &T, const double &p, const double &X) {
            ThermodynamicProperties props;
            UpdateState_TPX(props, T, p);
            return props;
        }
    };

};