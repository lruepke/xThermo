/**
 * @file DataStructures.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef DATASTRUCTURES_xThermal_H
#define DATASTRUCTURES_xThermal_H

#include "stdfunc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>

/**
 * @defgroup GENERAL_CONSTS Some general constants definition.
 * 
 * @{
 */
/** Char array length of a name, e.g. density */
#define LENGTH_CHAR_NAME      100
/** @} */

/**
 * @defgroup BitMaskRequestProp Bit mask of requesting partial derivatives of thermodynamic properties. See \cite PartialDer for definition.
 * 
 * @{
 */
/** do not request any properties */
#define Request_None     0
/** \f$2^{0}\f$: \f$ \left( \frac{\partial p}{\partial \rho} \right)_T \f$ */
#define Request_dpdrho_T 1
/** \f$2^{1}\f$:  \f$ \left( \frac{\partial p}{\partial T} \right)_{\rho} \f$ */
#define Request_dpdT_rho 2
/** \f$2^{2}\f$:  \f$ \left( \frac{\partial s}{\partial \rho} \right)_T \f$ */
#define Request_dsdrho_T 4
/** \f$2^{3}\f$:  \f$ \left( \frac{\partial s}{\partial T} \right)_{\rho} \f$ */
#define Request_dsdT_rho 8
/** \f$2^{4}\f$:  \f$ \left( \frac{\partial u}{\partial \rho} \right)_T \f$ */
#define Request_dudrho_T 16
/** \f$2^{5}\f$:  \f$ \left( \frac{\partial u}{\partial T} \right)_{\rho} \f$ */
#define Request_dudT_rho 32
/** \f$2^{6}\f$:  \f$ \left( \frac{\partial h}{\partial \rho} \right)_T \f$ */
#define Request_dhdrho_T 64
/** \f$2^{7}\f$:  \f$ \left( \frac{\partial h}{\partial T} \right)_{\rho} \f$ */
#define Request_dhdT_rho 128
/** \f$2^{8}\f$:  \f$ \left( \frac{\partial g}{\partial \rho} \right)_T \f$ */
#define Request_dgdrho_T 256
/** \f$2^{9}\f$:  \f$ \left( \frac{\partial g}{\partial T} \right)_{\rho} \f$ */
#define Request_dgdT_rho 512
/** \f$2^{10}\f$:  \f$ \left( \frac{\partial^2 p}{\partial \rho^2} \right)_{T} \f$ */
#define Request_ddpddrho_T 1024
/** \f$2^{11}\f$:  \f$ \left( \frac{\partial^2 p}{\partial T^2} \right)_{\rho} \f$ */
#define Request_ddpddT_rho 2048
/** \f$2^{12}\f$:  \f$ \left( \frac{\partial^2 p}{\partial \rho \partial T} \right) \f$ */
#define Request_ddpdrhodT 4096
/** \f$2^{13}\f$:  \f$ \left( \frac{\partial^2 s}{\partial \rho^2} \right)_{T} \f$ */
#define Request_ddsddrho_T 8192
/** \f$2^{14}\f$:  \f$ \left( \frac{\partial^2 s}{\partial T^2} \right)_{\rho} \f$ */
#define Request_ddsddT_rho 16384
/** \f$2^{15}\f$:  \f$ \left( \frac{\partial^2 s}{\partial \rho \partial T} \right) \f$ */
#define Request_ddsdrhodT 32768
/** \f$2^{16}\f$:  \f$ \left( \frac{\partial^2 u}{\partial \rho^2} \right)_{T} \f$ */
#define Request_dduddrho_T 65536
/** \f$2^{17}\f$:  \f$ \left( \frac{\partial^2 u}{\partial T^2} \right)_{\rho} \f$ */
#define Request_dduddT_rho 131072
/** \f$2^{18}\f$:  \f$ \left( \frac{\partial^2 u}{\partial \rho \partial T} \right) \f$ */
#define Request_ddudrhodT 262144
/** \f$2^{19}\f$:  \f$ \left( \frac{\partial^2 h}{\partial \rho^2} \right)_{T} \f$ */
#define Request_ddhddrho_T 524288
/** \f$2^{20}\f$:  \f$ \left( \frac{\partial^2 h}{\partial T^2} \right)_{\rho} \f$ */
#define Request_ddhddT_rho 1048576
/** \f$2^{21}\f$:  \f$ \left( \frac{\partial^2 h}{\partial \rho \partial T} \right) \f$ */
#define Request_ddhdrhodT 2097152
/** \f$2^{22}\f$:  \f$ \left( \frac{\partial^2 g}{\partial \rho^2} \right)_{T} \f$ */
#define Request_ddgddrho_T 4194304
/** \f$2^{23}\f$:  \f$ \left( \frac{\partial^2 g}{\partial T^2} \right)_{\rho} \f$ */
#define Request_ddgddT_rho 8388608
/** \f$2^{24}\f$:  \f$ \left( \frac{\partial^2 g}{\partial \rho \partial T} \right) \f$ */
#define Request_ddgdrhodT 16777216
//Transformation of partial derivatives. See Table 5 of \cite PartialDer
/** \f$2^{25}\f$:  \f$ \left( \frac{\partial \rho}{\partial p} \right)_{T} \f$*/
#define Request_drhodp_T    33554432
/** \f$2^{26}\f$:  \f$ \left( \frac{\partial h}{\partial p} \right)_{T} \f$*/
#define Request_dhdp_T      67108864
/** \f$2^{27}\f$:  \f$ \left( \frac{\partial \rho}{\partial T} \right)_{p} \f$*/
#define Request_drhodT_p    134217728
/** \f$2^{28}\f$:  \f$ \left( \frac{\partial h}{\partial T} \right)_{p} \f$*/
#define Request_dhdT_p      268435456
/** \f$2^{29}\f$:  \f$ \left( \frac{\partial p}{\partial \rho} \right)_{s} \f$*/
#define Request_dpdrho_s    536870912
/** \f$2^{30}\f$:  \f$ \left( \frac{\partial T}{\partial p} \right)_{h} \f$*/
#define Request_dTdp_h      1073741824
/** \f$2^{31}\f$:  \f$ \left( \frac{\partial \rho}{\partial p} \right)_{h} \f$*/
#define Request_drhodp_h    2147483648
/** \f$2^{32}\f$:  \f$ \left( \frac{\partial \rho}{\partial h} \right)_{p} \f$*/
#define Request_drhodh_p    4294967296
// Transformation of second-order partial derivatives. See Table 6 of \cite PartialDer.
/** \f$2^{33}\f$:  \f$ \left( \frac{\partial^2 T}{\partial p^2} \right)_{\rho} \f$*/
#define Request_ddTddp_rho      8589934592
/** \f$2^{34}\f$:  \f$ \left( \frac{\partial^2 T}{\partial \rho^2} \right)_{p} \f$*/
#define Request_ddTddrho_p      17179869184
/** \f$2^{35}\f$:  \f$ \left( \frac{\partial^2 T}{\partial p \partial \rho} \right) \f$*/
#define Request_ddTdpdrho       34359738368
/** \f$2^{36}\f$:  \f$ \left( \frac{\partial^2 \rho}{\partial p^2} \right)_{T} \f$*/
#define Request_ddrhoddp_T      68719476736
/** \f$2^{37}\f$:  \f$ \left( \frac{\partial^2 \rho}{\partial T^2} \right)_{p} \f$*/
#define Request_ddrhoddT_p      137438953472
/** \f$2^{38}\f$:  \f$ \left( \frac{\partial^2 \rho}{\partial T \partial p} \right) \f$*/
#define Request_ddrhodTdp       274877906944
/** \f$2^{39}\f$:  \f$ \left( \frac{\partial c_p}{\partial T} \right)_{\rho} \f$*/
#define Request_dcpdT_rho       549755813888
/** \f$2^{40}\f$:  \f$ \left( \frac{\partial c_p}{\partial \rho} \right)_{T} \f$*/
#define Request_dcpdrho_T       1099511627776
/** \f$2^{41}\f$:  \f$ \left( \frac{\partial c_p}{\partial p} \right)_{h} \f$*/
#define Request_dcpdp_h         2199023255552
// Bit mask of requesting thermodynamic properties. See Table 3 of \cite IAPWS-95 for properties definition.
// \note There is no density request, because density will be solved anyway and it is included in the IAPWS95::State.
/** \f$2^{42}\f$:  Internal energy \f$ u \f$ */
#define Request_u   4398046511104
/** \f$2^{43}\f$:  Specific Entropy \f$ s \f$ */
#define Request_s   8796093022208
/** \f$2^{44}\f$:  Specific Enthalpy \f$ h \f$ */
#define Request_h   17592186044416
/** \f$2^{45}\f$:  Isochoric heat capacity \f$ c_v \f$ */
#define Request_cv  35184372088832
/** \f$2^{46}\f$:  Isobaric heat capacity \f$ c_p \f$ */
#define Request_cp  70368744177664
/** \f$2^{47}\f$:  Speed of sound \f$ w \f$ */
#define Request_w   140737488355328
/** \f$2^{48}\f$:  Specific Gibbs-energy \f$g \f$. See Table of \cite PartialDer */
#define Request_g   281474976710656
/** \f$2^{49}\f$:  Density \f$ \rho \f$ */
#define Request_rho   562949953421312
/** \f$2^{50}\f$:  Temperature \f$ T \f$ */
#define Request_T   1125899906842624
/** @} */

namespace xThermal
{
    /**
     * @brief Define flag of solving saturated pressure or saturated temperature. 
     * 
     */
    enum SOLVE_SATURATED_PorT
    {
        SOLVE_SATURATED_P,      /**< Solve \f$ p_{sat} \f$ */
        SOLVE_SATURATED_T       /**< Solve \f$ T_{sat} \f$ */
    };
    /**
     * @brief Which space of the EOS is calculated or defined. For H2O, the available space are Space_TP and Space_HP. 
     * For H2ONaCl, the available space are Space_TPX and Space_HPX.
     * 
     */
    enum Space_EOS {Space_TP, Space_HP, Space_TPX, Space_HPX};
    /**
     * @brief index of dimension vector.
     * 
     */
    enum dimensionType
    {
        MASS, LENGTH, TIME, TEMPERATURE, MOLES, CURRENT, LUMINOUS_INTENSITY
    };
    const int nDimensions = 7;
    /**
     * @brief Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is \f$ kg/m^3 \f$, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0].
     * This idea comes from OpenFOAM.
     *
     */
    struct PhysicalDimension
    {
        short value[nDimensions];
        // friend std::ostream& operator<<(std::ostream& os, const PhysicalDimension& rVar)
        // {
        //     return os << "Dimension: ["
        //               <<rVar.value[MASS]<<", "
        //               <<rVar.value[LENGTH]<<", "
        //               <<rVar.value[TIME]<<", "
        //               <<rVar.value[TEMPERATURE]<<", "
        //               <<rVar.value[MOLES]<<", "
        //               <<rVar.value[CURRENT]<<", "
        //               <<rVar.value[LUMINOUS_INTENSITY]<<"]";

        // }
        // bool operator== (const PhysicalDimension &rVar) const
        // {
        //     return (value[MASS] == rVar.value[MASS]) && (value[LENGTH]==rVar.value[LENGTH]) && (value[TIME]==rVar.value[TIME]) && (value[TEMPERATURE]==rVar.value[TEMPERATURE]) && (value[MOLES]==rVar.value[MOLES]) && (value[CURRENT]==rVar.value[CURRENT]) && (value[LUMINOUS_INTENSITY]==rVar.value[LUMINOUS_INTENSITY]);
        // }
        // PhysicalDimension operator* (const PhysicalDimension &rVar) const
        // {
        //     PhysicalDimension tmp;
        //     for (int i = 0; i < nDimensions; i++){tmp.value[i] = value[i] + rVar.value[i];}
        //     return tmp;
        // }
        // PhysicalDimension operator/ (const PhysicalDimension &rVar) const
        // {
        //     PhysicalDimension tmp;
        //     for (int i = 0; i < nDimensions; i++){tmp.value[i] = value[i] - rVar.value[i];}
        //     return tmp;
        // }
    };
    /**
     * @defgroup DIMENSION_COMMONLYUSED Predefined commonly used physical dimension.
     * 
     * @{
     */
    /** dimensionless */
#define dimension_SI_None           {0,  0,  0,  0,  0,  0,  0}
    /** SI unit of pressure [Pa] */
#define dimension_SI_P              {1, -1, -2,  0,  0,  0,  0}
    /** SI unit of temperature [K] */
#define dimension_SI_T              {0,  0,  0,  1,  0,  0,  0}
    /** SI unit of specific enthalpy [J/kg] */
#define dimension_SI_H              {0,  2, -2,  0,  0,  0,  0}
    /** SI unit of density [kg/m^3] */
#define dimension_SI_Rho            {1, -3,  0,  0,  0,  0,  0}
    /** SI unit of mass [kg] */
#define dimension_SI_mass           {1,  0,  0,  0,  0,  0,  0}
    /** SI unit of length [m] */
#define dimension_SI_length         {0,  1,  0,  0,  0,  0,  0}
    /** SI unit of time [s] */
#define dimension_SI_time           {0,  0,  1,  0,  0,  0,  0}
    /** @} */

    /**
     * @defgroup Initialization_CommonlyUsedProps Predefined initialization of commonly used thermodynamic properties.
     * 
     * @{
     */
    /** initialization of property pressure */
#define thermodynamicProperty_P     {0, "Pressure",             "$p$",      "Pa",    dimension_SI_P}
    /** initialization of property temperature */
#define thermodynamicProperty_T     {0, "Temperature",          "$T$",      "K",     dimension_SI_T}
    /** initialization of property specific enthalpy */
#define thermodynamicProperty_H     {0, "Specific enthalpy",    "$H$",      "J/kg",  dimension_SI_H}
    /** initialization of property density */
#define thermodynamicProperty_Rho   {0, "Density",              "$\rho$",   "kg/m3", dimension_SI_Rho}
    /** initialization of property vapor mass fraction: \f$ x=\frac{m^{\prime\prime}}{m^{\prime} + m^{\prime}} = \frac{h - h^{\prime}}{h^{\prime\prime} - h^{\prime}} =\frac{v - v^{\prime}}{v^{\prime\prime} - v^{\prime}} = \frac{1/\rho - 1/\rho^{\prime}}{1/\rho^{\prime\prime} - 1/\rho^{\prime}} \f$, where the superscript \f$\prime, \prime\prime\f$ denote liquid and vapor, respectively. */
#define thermodynamicProperty_X     {0, "Vapor mass fraction",  "$x$",      "1",     dimension_SI_None}
    /** Undefined property */
#define thermodynamicProperty_Undefined {0, "Undefined",  "Undefined",      "Undefined",     dimension_SI_None}
    /** @} */

    /**
     * @brief Data struct of a thermodynamic property.
     * 
     */
    struct ThermodynamicProperty
    {
        double value;                       /**< value of a property */
        char name[LENGTH_CHAR_NAME];        /**< name of a property */
        char name_math[LENGTH_CHAR_NAME];   /**< Latex format of the name */
        char unit_name[LENGTH_CHAR_NAME];   /**< Latex format of the unit */
        PhysicalDimension dimension;                 /**< Dimension of a property in OpenFOAM's format: [kg, m, s, K, mol, A, cd] */
    };


    /**
     * @brief Phase region index of H2O, NaCl, H2ONaCl system.
     * 
     */
    enum PhaseRegion
    {
        MixPhaseRegion=-1,      /**< MixPhaseRegion is used in AMR lookup table. */
        SinglePhase_L,          /**< Single phase: liquid. Regard single-phase "liquid-like" and "vapor-like" fluid as "Single phase liquid" as well. */
        SinglePhase_V,          /**< Single phase: vapor*/
        SinglePhase_S,          /**< Single phase: solid, e.g. NaCl */
        Supercritical,          /**< Supercritical (p > pc, T > Tc) for pure water */
        Supercritical_vapor,    /**< Supercritical vapor (p < pc, T > Tc) for pure water */
        Supercritical_liquid,   /**< Supercritical liquid (p > pc, T < Tc) for pure water */
        Critical,               /**< At the critical point of pure water or on critical curve of H2ONaCl */
        TwoPhase_VL_Water,           /**< Two phase: liquid + vapor, salinity = 0*/
        TwoPhase_LH,           /**< Two phase: liquid + halite */
        TwoPhase_VH,           /**< Two phase: vapor + halite */
        // TwoPhase_V_L_L,         /**< Two phase: vapor + liquid coexistence, liquid branch */
        // TwoPhase_V_L_V,         /**< Two phase: vapor + liquid coexistence, vapor branch */
        TwoPhase_VL,            /**< Two phase: VL region */
        ThreePhase_VLH,       /**< Three phase: vapor + liquid + halite */
        Unknown,                /**< Unknown phase region for bug report*/
        NotImposed              /**< Phase is not imposed */
    };
    /**
     * @brief Map of phase region index (enum) to phase region name (string).
     *
     */
    typedef std::map<PhaseRegion,std::string> MAP_PhaseRegion2Name;
    static MAP_PhaseRegion2Name map_phase2name={
            {MixPhaseRegion, "Mix phase region in AMR-LUT"},
            {SinglePhase_L,"Liquid"},
            {SinglePhase_V, "Vapor"},
            {SinglePhase_S, "Solid"},
            {Supercritical, "Supercritical"},
            {Supercritical_vapor, "Sup.crit. vapor"},
            {Supercritical_liquid, "Sup.crit. liquid"},
            {Critical, "Critical point"},
            {TwoPhase_VL_Water, "V+L(water)"},
            {TwoPhase_LH, "L+H"},
            {TwoPhase_VH, "V+H"},
            {TwoPhase_VL, "V+L"},
            {ThreePhase_VLH, "V+L+H"},
            {Unknown, "Unknown"},
            {NotImposed, "Phase is not imposed"}
    };

    /**
         * only used as input argument for \link Rho_phase \endlink;
         */
    enum PhaseType{Vapor, Liquid};



    enum INPUT_PAIR{TPX, HPX};
    struct ThermodynamicProperties
    {
        double T=0, p=0, X=0, H=0;
        PhaseRegion phase=Unknown;
        double S_l=0, S_v=0, S_h=0; // volume fraction for H2O-NaCl, but mass fraction for H2O. For two-phase of pure water convection modeling, need to calculate it using, e.g., Eq. 37 of Vehling et al.(2018).
        double X_l=0, X_v=0;
        double Rho_l=0, Rho_v=0, Rho_h=0;
        double H_l=0, H_v=0, H_h=0;
        double Cp_l = 0, Cp_v = 0, Cp_h = 0;       // dH/dT
        double dHdP_l = 0, dHdP_v = 0, dHdP_h = 0; // dH/dP
        double Mu_l = 0, Mu_v = 0;
        double Rho=0, Cp = 0, Mu = 1;
        // derivative
        double dRhodP = 0, dRhodP_l = 0, dRhodP_v = 0;
        double dRhodT = 0, dRhodT_l = 0, dRhodT_v = 0;
        // commonly used derivatives for hydrothermal convection simulation
        double IsothermalCompressibility = 0, IsothermalCompressibility_l =0, IsothermalCompressibility_v =0;
        double IsobaricExpansivity=0, IsobaricExpansivity_l =0, IsobaricExpansivity_v =0;
        std::string fluidName;
        // std out
        friend std::ostream& operator<<(std::ostream& os, const ThermodynamicProperties& self) {
            std::stringstream sin;
            sin<<self.fluidName<<std::endl;
            if (self.fluidName=="H2O-NaCl")
            {
                if (self.phase==SinglePhase_L || self.phase==SinglePhase_S || self.phase==SinglePhase_V || self.phase==Supercritical || self.phase==Supercritical_liquid || self.phase==Supercritical_vapor)
                {
                    sin<<"T="<<self.T-273.15<<" deg.C = "<<self.T<<" K, P="<<self.p/1E5<<" bar, X="<<self.X*100<<" wt% NaCl, H="<<self.H<<" J/kg"<<std::endl;
                    sin<<"  phase: "<<map_phase2name[self.phase]<<std::endl;
                    sin<<"  Rho="<<self.Rho<<" [kg/m^3]"<<std::endl;
                    sin<<"  H="<<self.H<<" [J/kg]"<<std::endl;
                    sin<<"  Cp="<<self.Cp<<" [J/kg/K]"<<std::endl;
                    sin<<"  Mu="<<self.Mu<<" [Pa S]"<<std::endl;
                    sin<<"  dRhodT="<<self.dRhodT<<" [-]"<<std::endl;
                    sin<<"  dRhodP="<<self.dRhodP<<" [-]"<<std::endl;
                    sin<<"  Isothermal compressibility="<<self.IsothermalCompressibility<<" [1/Pa]"<<std::endl;
                    sin<<"  Isobaric expansivity="<<self.IsobaricExpansivity<<" [1/K]"<<std::endl;
                } else
                {
                    sin<<"T="<<self.T-273.15<<" deg.C = "<<self.T<<" K, P="<<self.p/1E5<<" bar, X="<<self.X*100<<" wt% NaCl, H="<<self.H<<" J/kg"<<std::endl;
                    sin<<"  phase: "<<map_phase2name[self.phase]<<std::endl;
                    sin<<"  S_l="<<self.S_l<<", S_v="<<self.S_v<<", S_h="<<self.S_h<<std::endl;
                    sin<<"  X_l="<<self.X_l<<", X_v: "<<self.X_v<<std::endl;
                    sin<<"  Rho="<<self.Rho<<", Rho_l="<<self.Rho_l<<", Rho_v="<<self.Rho_v<<", Rho_h="<<self.Rho_h<<" [kg/m^3]"<<std::endl;
                    sin<<"  H="<<self.H<<", H_l="<<self.H_l<<", H_v="<<self.H_v<<", H_h="<<self.H_h<<" [J/kg]"<<std::endl;
                    sin<<"  Cp_l="<<self.Cp_l<<", Cp_v="<<self.Cp_v<<", Cp_h="<<self.Cp_h<<" [J/kg/K]"<<std::endl;
                    sin<<"  Mu_l="<<self.Mu_l<<", Mu_v="<<self.Mu_v<<" [Pa S]"<<std::endl;
                    // sin<<"  dRhodT="<<self.dRhodT<<", dRhodT_l="<<self.dRhodT_l<<", dRhodT_v="<<self.dRhodT_v<<" [kg/m3/K]"<<std::endl;
                    // sin<<"  dRhodP="<<self.dRhodP<<", dRhodP_l="<<self.dRhodP_l<<", dRhodP_v="<<self.dRhodP_v<<" [kg/m3/Pa]"<<std::endl;
                    sin<<"  Isothermal compressibility="<<self.IsothermalCompressibility<<" [1/Pa]"<<std::endl;
                    // sin<<"  Isobaric expansivity="<<self.IsobaricExpansivity<<" [1/K]"<<std::endl;
                }

            } else if(self.fluidName=="IAPS84" || self.fluidName=="IAPWS95" || self.fluidName=="IAPWS95_CoolProp" || self.fluidName == "NaCl")
            {
                if (self.phase==SinglePhase_L || self.phase==SinglePhase_S || self.phase==SinglePhase_V || self.phase==Supercritical || self.phase==Supercritical_liquid || self.phase==Supercritical_vapor)
                {
                    sin<<"T="<<self.T-273.15<<" deg.C, P="<<self.p/1E5<<" bar, H="<<self.H<<" J/kg"<<std::endl;
                    sin<<"  phase: "<<map_phase2name[self.phase]<<std::endl;
                    sin<<"  Rho="<<self.Rho<<" [kg/m^3]"<<std::endl;
                    sin<<"  H="<<self.H<<" [J/kg]"<<std::endl;
                    sin<<"  Cp="<<self.Cp<<" [J/kg/K]"<<std::endl;
                    sin<<"  Mu="<<self.Mu<<" [Pa S]"<<std::endl;
                    sin<<"  dRhodT="<<self.dRhodT<<" [kg/m3/K]"<<std::endl;
                    sin<<"  dRhodP="<<self.dRhodP<<" [kg/m3/Pa]"<<std::endl;
                    sin<<"  Isothermal compressibility="<<self.IsothermalCompressibility<<" [1/Pa]"<<std::endl;
                    sin<<"  Isobaric expansivity="<<self.IsobaricExpansivity<<" [1/K]"<<std::endl;
                } else
                {
                    sin<<"T="<<self.T-273.15<<" deg.C, P="<<self.p/1E5<<" bar, H="<<self.H<<" J/kg"<<std::endl;
                    sin<<"  phase: "<<map_phase2name[self.phase]<<std::endl;
                    sin<<"  S_l="<<self.S_l<<", S_v="<<self.S_v<<std::endl;
                    sin<<"  Rho="<<self.Rho<<", Rho_l="<<self.Rho_l<<", Rho_v="<<self.Rho_v<<" [kg/m^3]"<<std::endl;
                    sin<<"  H="<<self.H<<", H_l="<<self.H_l<<", H_v="<<self.H_v<<" [J/kg]"<<std::endl;
                    sin<<"  Cp_l="<<self.Cp_l<<", Cp_v="<<self.Cp_v<<" [J/kg/K]"<<std::endl;
                    sin<<"  Mu_l="<<self.Mu_l<<", Mu_v="<<self.Mu_v<<" [Pa S]"<<std::endl;
                    // sin<<"  dRhodT="<<self.dRhodT<<", dRhodT_l="<<self.dRhodT_l<<", dRhodT_v="<<self.dRhodT_v<<" [kg/m3/K]"<<std::endl;
                    // sin<<"  dRhodP="<<self.dRhodP<<", dRhodP_l="<<self.dRhodP_l<<", dRhodP_v="<<self.dRhodP_v<<" [kg/m3/Pa]"<<std::endl;
                    sin<<"  Isothermal compressibility="<<self.IsothermalCompressibility<<" [1/Pa]"<<std::endl;
                    // sin<<"  Isobaric expansivity="<<self.IsobaricExpansivity<<" [1/K]"<<std::endl;
                }

            }else
            {
                ERROR("The fluidName is not set in ThermodynamicProperties or is not identified: "+self.fluidName);
            }

            return os << sin.str();
        }
        // info for Python API
        static std::string info(const ThermodynamicProperties& self)
        {
            std::stringstream sin;
            sin<<self;
            return sin.str();
        }
    };
    // for matlab API or OpenFOAM
    struct ThermodynamicPropertiesArray
    {
        size_t N=0; // length of array
        // INPUT_PAIR input_pair;
        double *T=nullptr, *p=nullptr, *X=nullptr, *H=nullptr;
        double *phase=nullptr;
        double *S_l=nullptr, *S_v=nullptr, *S_h=nullptr;
        double *X_l=nullptr, *X_v=nullptr;
        double *Rho_l=nullptr, *Rho_v=nullptr, *Rho_h=nullptr;
        double *H_l=nullptr, *H_v=nullptr, *H_h=nullptr;
        double *Mu_l = nullptr, *Mu_v = nullptr;
        double *Cp_l = nullptr, *Cp_v = nullptr, *Cp_h = nullptr;
        double *Rho=nullptr, *Mu = nullptr, *Cp = nullptr;
        std::string fluidName;
        void fill(const ThermodynamicProperties& props, const size_t& i) const
        {
            T[i] = props.T;
            p[i] = props.p;
            X[i] = props.X;
            H[i] = props.H;
            phase[i] = props.phase;
            S_l[i] = props.S_l;
            S_v[i] = props.S_v;
            S_h[i] = props.S_h;
            X_l[i] = props.X_l;
            X_v[i] = props.X_v;
            Rho_l[i] = props.Rho_l;
            Rho_v[i] = props.Rho_v;
            Rho_h[i] = props.Rho_h;
            H_l[i] = props.H_l;
            H_v[i] = props.H_v;
            H_h[i] = props.H_h;
            Cp_l[i] = props.Cp_l;
            Cp_v[i] = props.Cp_v;
            Cp_h[i] = props.Cp_h;
            Mu_l[i] = props.Mu_l;
            Mu_v[i] = props.Mu_v;
            Rho[i] = props.Rho;
            Cp[i] = props.Cp;
            Mu[i] = props.Mu;
        }
        void create(const size_t N0)
        {
            // clean first
            clean();
            N = N0;
            T = new double[N];
            p = new double[N];
            X = new double[N];
            H = new double[N];
            phase = new double[N];
            S_l = new double[N];
            S_v = new double[N];
            S_h = new double[N];
            X_l = new double[N];
            X_v = new double[N];
            Rho_l = new double[N];
            Rho_v = new double[N];
            Rho_h = new double[N];
            H_l = new double[N];
            H_v = new double[N];
            H_h = new double[N];
            Mu_l = new double[N];
            Mu_v = new double[N];
            Cp_l = new double[N];
            Cp_v = new double[N];
            Cp_h = new double[N];
            Rho = new double[N];
            Mu = new double[N];
            Cp = new double[N];
        }
        void clean()
        {
            if(T) delete[] T, T = nullptr;
            if(p) delete[] p, p = nullptr;
            if(X) delete[] X, X = nullptr;
            if(H) delete[] H, H = nullptr;
            if(phase) delete[] phase, phase = nullptr;
            if(S_l) delete[] S_l, S_l = nullptr;
            if(S_v) delete[] S_v, S_v = nullptr;
            if(S_h) delete[] S_h, S_h = nullptr;
            if(X_l) delete[] X_l, X_l = nullptr;
            if(X_v) delete[] X_v, X_v = nullptr;
            if(Rho_l) delete[] Rho_l, Rho_l = nullptr;
            if(Rho_v) delete[] Rho_v, Rho_v = nullptr;
            if(Rho_h) delete[] Rho_h, Rho_h = nullptr;
            if(H_l) delete[] H_l, H_l = nullptr;
            if(H_v) delete[] H_v, H_v = nullptr;
            if(H_h) delete[] H_h, H_h = nullptr;
            if(Mu_l) delete[] Mu_l, Mu_l = nullptr;
            if(Mu_v) delete[] Mu_v, Mu_v = nullptr;
            if(Cp_l) delete[] Cp_l, Cp_l = nullptr;
            if(Cp_v) delete[] Cp_v, Cp_v = nullptr;
            if(Cp_h) delete[] Cp_h, Cp_h = nullptr;
            if(Rho) delete[] Rho, Rho = nullptr;
            if(Mu) delete[] Mu, Mu = nullptr;
            if(Cp) delete[] Cp, Cp = nullptr;
        }
    };
    // for python API
    struct ThermodynamicPropertiesVector
    {
        std::vector<double> T, p, X, H;
        std::vector<int> phase;
        std::vector<double> S_l, S_v, S_h;
        std::vector<double> X_l, X_v;
        std::vector<double> Rho_l, Rho_v, Rho_h;
        std::vector<double> H_l, H_v, H_h;
        std::vector<double> Cp_l, Cp_v, Cp_h;
        std::vector<double> Mu_l, Mu_v;
        std::vector<double> Rho, Cp, Mu;
        // derivative
        std::vector<double> dRhodP, dRhodP_l, dRhodP_v;
        std::vector<double> dRhodT, dRhodT_l, dRhodT_v;
        // commonly used derivatives for hydrothermal convection simulation
        std::vector<double> IsothermalCompressibility, IsothermalCompressibility_l, IsothermalCompressibility_v;
        std::vector<double> IsobaricExpansivity, IsobaricExpansivity_l, IsobaricExpansivity_v;
        std::string fluidName;
        // std out
        friend std::ostream& operator<<(std::ostream& os, const ThermodynamicPropertiesVector& self) {
            std::stringstream sin;
            sin<<"ThermodynamicPropertiesVector, size: "<<self.T.size()<<std::endl;
            sin<<"T     ["<<self.T[0]<<", ..., "<<self.T[self.T.size()-1]<<"]  [K]"<<std::endl;
            sin<<"P     ["<<self.p[0]<<", ..., "<<self.p[self.p.size()-1]<<"]  [Pa]"<<std::endl;
            sin<<"X     ["<<self.X[0]<<", ..., "<<self.X[self.X.size()-1]<<"]  [kg/kg]"<<std::endl;
            sin<<"H     ["<<self.H[0]<<", ..., "<<self.H[self.H.size()-1]<<"]  [J/kg]"<<std::endl;
            sin<<"phase ["<<self.phase[0]<<", ..., "<<self.phase[self.phase.size()-1]<<"]  [-]"<<std::endl;
            sin<<"S_l   ["<<self.S_l[0]<<", ..., "<<self.S_l[self.S_l.size()-1]<<"]  [-]"<<std::endl;
            sin<<"IsothermalCompressibility   ["<<self.IsothermalCompressibility[0]<<", ..., "<<self.IsothermalCompressibility[self.IsothermalCompressibility.size()-1]<<"]  [1/Pa]"<<std::endl;
            sin<<"IsobaricExpansivity   ["<<self.IsobaricExpansivity[0]<<", ..., "<<self.IsobaricExpansivity[self.IsobaricExpansivity.size()-1]<<"]  [1/K]"<<std::endl;
            sin<<"... ..."<<std::endl;
            sin<<"S_v, S_h, X_l, X_v, Rho_l, Rho_v, Rho_h, H_l, H_v, H_h, Cp_l, Cp_v, Cp_h, Mu_l, Mu_v, Rho, Cp, Mu, dRhodP, dRhodT"<<std::endl;
            return os << sin.str();
        }
        void write(std::string filename)
        {
            std::string extname = extname_file(filename);
            if (extname=="csv")
            {
                write_csv(filename);
            } else
            {
                WARNING("Unidentified file format: "<<COLOR_RED<<extname<<COLOR_DEFAULT<<", use default csv");
                write_csv(filename+".csv");
            }
        }
        void write_vtk(const std::string& filename, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::string& xTitle, const std::string& yTitle, const std::string& zTitle, bool isNormalize=false)
        {
            using namespace std;
            STATUS("Writing ThermodynamicPropertiesVector to structured vtk grid file : "<<filename);
            std::ofstream fpout(filename);
            if(!fpout)
            {
                ERROR("Can not open file: "<<filename);
                exit(0);
            }
            string fname_py=filename+".py";
            //  write vtk head
            fpout<<"# vtk DataFile Version 2.0"<<endl;
            fpout<<"Properties of seawater"<<endl;
            fpout<<"ASCII"<<endl;
            fpout<<"DATASET RECTILINEAR_GRID"<<endl;
            fpout<<"DIMENSIONS "<<x.size()<<" "<<y.size()<<" "<<z.size()<<endl;
            double len_x=1, len_y=1, len_z=1;
            double xMAX=*std::max_element(x.begin(), x.end());
            double xMIN=*std::min_element(x.begin(), x.end());
            double yMAX=*std::max_element(y.begin(), y.end());
            double yMIN=*std::min_element(y.begin(), y.end());
            double zMAX=*std::max_element(z.begin(), z.end());
            double zMIN=*std::min_element(z.begin(), z.end());
            double scale_x=1, scale_y=1, scale_z=1;
            len_x=(xMAX==xMIN ? 1: xMAX-xMIN);
            len_y=(yMAX==yMIN ? 1: yMAX-yMIN);
            len_z=(zMAX==zMIN ? 1: zMAX-zMIN);
            if(!isNormalize) // if not normalize the vtk file, write a python script to better visualize result
            {
                scale_y=len_x/len_y;
                scale_z=len_x/len_z;
                std::ofstream fout_py(fname_py);
                if(!fout_py)
                {
                    cout<<"Warning: cannot generate pvPython script for Paraview. "<<fname_py<<endl;
                }else
                {
                    fout_py<<"from paraview.simple import *"<<endl;
                    fout_py<<"xHvtk = LegacyVTKReader(FileNames=[\'"<<filename<<"\'])"<<endl;
                    fout_py<<"renderView1 = GetActiveViewOrCreate(\'RenderView\')"<<endl;
                    fout_py<<"xHvtkDisplay = Show(xHvtk, renderView1)"<<endl;
                    fout_py<<"xHvtkDisplay.Representation = \'Surface\'"<<endl;
                    fout_py<<"renderView1.AxesGrid.Visibility = 1"<<endl;
                    fout_py<<"xHvtkDisplay.Scale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                    fout_py<<"renderView1.AxesGrid.DataScale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                    // fout_py<<"renderView1.AxesGrid.DataBoundsInflateFactor = 0"<<endl;
                    fout_py<<"renderView1.AxesGrid.XTitle = \'"<<xTitle<<"\'"<<endl;
                    fout_py<<"renderView1.AxesGrid.YTitle = \'"<<yTitle<<"\'"<<endl;
                    fout_py<<"renderView1.AxesGrid.ZTitle = \'"<<zTitle<<"\'"<<endl;
                    if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleFontSize = 16"<<endl;
                    if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleBold = 1"<<endl;
                    if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleFontSize = 16"<<endl;
                    if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleBold = 1"<<endl;
                    if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleFontSize = 16"<<endl;
                    if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleBold = 1"<<endl;
                    // set default data source as PhaseRegion
                    fout_py<<"#set default data source as PhaseRegion"<<endl;
                    fout_py<<"paraview.simple._DisableFirstRenderCameraReset()"<<endl;
                    fout_py<<"legacyVTKReader1 = GetActiveSource()"<<endl;
                    fout_py<<"renderView1 = GetActiveViewOrCreate('RenderView')"<<endl;
                    fout_py<<"legacyVTKReader1Display = GetDisplayProperties(legacyVTKReader1, view=renderView1)"<<endl;
                    fout_py<<"ColorBy(legacyVTKReader1Display, ('POINTS', 'PhaseRegion'))"<<endl;
                    fout_py<<"legacyVTKReader1Display.RescaleTransferFunctionToDataRange(True, False)"<<endl;
                    fout_py<<"legacyVTKReader1Display.SetScalarBarVisibility(renderView1, True)"<<endl;
                    fout_py<<"phaseRegionLUT = GetColorTransferFunction('PhaseRegion')\n"<<endl;
                    fout_py<<"renderView1.ResetCamera()"<<endl;
                    fout_py.close();
                    STATUS("Paraview-python script is generated as : "<<fname_py);
                }
            }
            fpout<<"X_COORDINATES "<<x.size()<<" float"<<endl;
            if(isNormalize)
            {
                for(size_t i=0;i<x.size();i++)fpout<<(x[i]-xMIN)/len_x<<" ";
                fpout<<endl;
            }else
            {
                for(size_t i=0;i<x.size();i++)fpout<<x[i]<<" ";
                fpout<<endl;
            }
            fpout<<"Y_COORDINATES "<<y.size()<<" float"<<endl;
            if(isNormalize)
            {
                for(size_t i=0;i<y.size();i++)fpout<<(y[i]-yMIN)/len_y<<" ";
                fpout<<endl;
            }else
            {
                for(size_t i=0;i<y.size();i++)fpout<<y[i]<<" ";
                fpout<<endl;
            }
            fpout<<"Z_COORDINATES "<<z.size()<<" float"<<endl;
            if(isNormalize)
            {
                for(size_t i=0;i<z.size();i++)fpout<<(z[i]-zMIN)/len_z<<" ";
                fpout<<endl;
            }else
            {
                for(size_t i=0;i<z.size();i++)fpout<<z[i]<<" ";
                fpout<<endl;
            }
            fpout<<"POINT_DATA "<<T.size()<<endl;
            // 1. phase region
            fpout<<"SCALARS PhaseRegion int"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(int var : phase)fpout<<var<<" ";
            fpout<<endl;
            // T
            fpout<<"SCALARS T double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : T)fpout<<var<<" ";
            fpout<<endl;
            // 2. H
            fpout<<"SCALARS H double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : H)fpout<<var<<" ";
            fpout<<endl;
            // 3. Rho
            fpout<<"SCALARS Rho double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : Rho)fpout<<var<<" ";
            fpout<<endl;
            // Cp
            fpout<<"SCALARS Cp double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : Cp)fpout<<var<<" ";
            // Mu
            fpout<<"SCALARS Mu double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : Mu)fpout<<var<<" ";
            // Isothermal compressibility
            fpout<<"SCALARS IsothermalCompressibility double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : IsothermalCompressibility)fpout<<var<<" ";
            // Isobaric expansivity
            fpout<<"SCALARS IsobaricExpansivity double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : IsobaricExpansivity)fpout<<var<<" ";
            // S_l
            fpout<<"SCALARS S_l double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for(double var : S_l)fpout<<var<<" ";
            fpout<<endl;

            fpout.close();
            STATUS("You can use command of "<<COLOR_PURPLE<<"paraview --script="<<fname_py<<COLOR_DEFAULT<<" to visualize result in paraview")
        }
        void write_csv(const std::string& filename)
        {
            std::ofstream fout(filename);
            if(!fout)
            {
                if(!fout) ERROR("Open file failed in ThermodynamicPropertiesVector.write: "+filename);
                exit(0);
            }
            if (fluidName=="H2O-NaCl")
            {
                fout<<"T[C], P[Pa], X[kg/kg], Xl, Xv, Phase[index], Phase[name], Rho[kg/m3], Rhol, Rhov, Rhoh, H[J/kg], Hl, Hv, Hh, Cp[J/kg/K], Cpl, Cpv, Cph, Mu[Pa S], Mul, Muv, Sl[kg/kg], Sv, Sh"<<std::endl;
                for (size_t i = 0; i < T.size(); ++i) {
                    fout<<T[i]<<", "<<p[i]<<", "<<X[i]<<", "<<X_l[i]<<", "<<X_v[i]<<", "<<phase[i]<<", "<<map_phase2name[(PhaseRegion)phase[i]]<<", "
                        <<Rho[i]<<", "<<Rho_l[i]<<", "<<Rho_v[i]<<", "<<Rho_h[i]<<", "
                        <<H[i]<<", "<<H_l[i]<<", "<<H_v[i]<<", "<<H_h[i]<<", "
                        <<Cp[i]<<", "<<Cp_l[i]<<", "<<Cp_v[i]<<", "<<Cp_h[i]<<", "
                        <<Mu[i]<<", "<<Mu_l[i]<<", "<<Mu_v[i]<<", "
                        <<S_l[i]<<", "<<S_v[i]<<", "<<S_h[i]
                        <<std::endl;
                }
            } else if(fluidName=="IAPS84" || fluidName=="IAPWS95" || fluidName=="IAPWS95_CoolProp" || fluidName == "NaCl")
            {
                fout<<"T[C], P[Pa], Phase[index], Phase[name], Rho[kg/m3], Rhol, Rhov, H[J/kg], Hl, Hv, Cp[J/kg/K], Cpl, Cpv, Mu[Pa S], Mul, Muv, Sl[kg/kg], Sv"<<std::endl;
                for (size_t i = 0; i < T.size(); ++i) {
                    fout<<T[i]<<", "<<p[i]<<", "<<phase[i]<<", "<<map_phase2name[(PhaseRegion)phase[i]]<<", "
                        <<Rho[i]<<", "<<Rho_l[i]<<", "<<Rho_v[i]<<", "
                        <<H[i]<<", "<<H_l[i]<<", "<<H_v[i]<<", "
                        <<Cp[i]<<", "<<Cp_l[i]<<", "<<Cp_v[i]<<", "
                        <<Mu[i]<<", "<<Mu_l[i]<<", "<<Mu_v[i]<<", "
                        <<S_l[i]<<", "<<S_v[i]<<", "
                        <<std::endl;
                }
            }else
            {
                ERROR("The fluidName is not set in ThermodynamicPropertiesVector or is not identified: "+fluidName);
            }
            fout.close();
            STATUS("The file has been saved: "+filename)
        }
        void fill(const ThermodynamicProperties& props, const size_t& i)
        {
            T[i] = props.T;
            p[i] = props.p;
            X[i] = props.X;
            H[i] = props.H;
            phase[i] = props.phase;
            S_l[i] = props.S_l;
            S_v[i] = props.S_v;
            S_h[i] = props.S_h;
            X_l[i] = props.X_l;
            X_v[i] = props.X_v;
            Rho_l[i] = props.Rho_l;
            Rho_v[i] = props.Rho_v;
            Rho_h[i] = props.Rho_h;
            H_l[i] = props.H_l;
            H_v[i] = props.H_v;
            H_h[i] = props.H_h;
            Cp_l[i] = props.Cp_l;
            Cp_v[i] = props.Cp_v;
            Cp_h[i] = props.Cp_h;
            Mu_l[i] = props.Mu_l;
            Mu_v[i] = props.Mu_v;
            Rho[i] = props.Rho;
            Cp[i] = props.Cp;
            Mu[i] = props.Mu;
            // derivatives
            dRhodP[i] = props.dRhodP;
            dRhodP_l[i] = props.dRhodP_l;
            dRhodP_v[i] = props.dRhodP_v;
            dRhodT[i] = props.dRhodT;
            dRhodT_l[i] = props.dRhodT_l;
            dRhodT_v[i] = props.dRhodT_v;
            IsothermalCompressibility[i] = props.IsothermalCompressibility;
            IsothermalCompressibility_l[i] = props.IsothermalCompressibility_l;
            IsothermalCompressibility_v[i] = props.IsothermalCompressibility_v;
            IsobaricExpansivity[i] = props.IsobaricExpansivity;
            IsobaricExpansivity_l[i] = props.IsobaricExpansivity_l;
            IsobaricExpansivity_v[i] = props.IsobaricExpansivity_v;
        }
        void resize(size_t N)
        {
            phase.resize(N,Unknown);
            T.resize(N, 0);
            p.resize(N, 0);
            X.resize(N, 0);
            H.resize(N, 0);
            S_l.resize(N, 0);
            S_v.resize(N, 0);
            S_h.resize(N, 0);
            X_l.resize(N, 0);
            X_v.resize(N, 0);
            Rho_l.resize(N, 0);
            Rho_v.resize(N, 0);
            Rho_h.resize(N, 0);
            H_l.resize(N, 0);
            H_v.resize(N, 0);
            H_h.resize(N, 0);
            Cp_l.resize(N, 0);
            Cp_v.resize(N, 0);
            Cp_h.resize(N, 0);
            Mu_l.resize(N, 0);
            Mu_v.resize(N, 0);
            Rho.resize(N, 0);
            Cp.resize(N, 0);
            Mu.resize(N, 0);
            // derivatives
            dRhodP.resize(N, 0);
            dRhodP_l.resize(N, 0);
            dRhodP_v.resize(N, 0);
            dRhodT.resize(N, 0);
            dRhodT_l.resize(N, 0);
            dRhodT_v.resize(N, 0);
            IsothermalCompressibility.resize(N, 0);
            IsothermalCompressibility_l.resize(N, 0);
            IsothermalCompressibility_v.resize(N, 0);
            IsobaricExpansivity.resize(N, 0);
            IsobaricExpansivity_l.resize(N, 0);
            IsobaricExpansivity_v.resize(N, 0);
        }
        // info for Python API
        std::string info(const ThermodynamicPropertiesVector& self)
        {
            std::stringstream sin;
            sin<<self;
            return sin.str();
        }
    };
    /**
     * 2D Deformable linear "structured" mesh. Related to matplotlib ax.plot_wireframe()
     */
    struct DeformLinearMesh
    {
        std::vector<std::vector<double> > T, p, X;
    };

    /**
     * 2D triangle mesh structure. Relate to matplotlib ax.tri_plot() in 2D or ax.plot_trisurf() in 3D
     */
    struct TriMesh
    {
        std::vector<double> x,y,z;
        std::vector<std::vector<int> > connection; //size is 3
        std::vector<std::string>names={"x","y","z"};
        void append(const TriMesh& trimesh)
        {
            int baseIndex = (int)x.size();
            // append connection
            for (const auto & i : trimesh.connection) {
                connection.push_back({i[0] + baseIndex, i[1] + baseIndex, i[2] + baseIndex});
            }
            // append points
            for (size_t i = 0; i < trimesh.x.size(); ++i) {
                x.push_back(trimesh.x[i]);
                y.push_back(trimesh.y[i]);
                z.push_back(trimesh.z[i]);
            }
        }
    };
    /**
    * Data structure for a phase boundary surface
    */
    struct Surface
    {
        std::string name, shortName;
        DeformLinearMesh mesh;
        std::string color;
    };
    /**
    * Data structure for a phase boundary line
    */
    struct Line
    {
        std::string name,shortName;
        std::vector<double > T,p,X;
        std::string color = "black";
    };
    /**
    * Data structure for a phase boundary point
    */
    struct Point
    {
        std::string name, shortName;
        double T,p,X;
        std::string color = "black";
    };
#define SCALE_X_linear 1
#define SCALE_X_log 2
#define SCALE_X_loglinear 3
    /**
     * Struct of phase boundaries. Boundary surface, lines and points.
     */
    struct PhaseBoundaries
    {
        std::vector<Surface> surfaces;
        std::vector<Line> lines;
        std::vector<Point> points;
        int scale_X; //"linear", "log", "loglinear"
        double Xcenter; //The log and linear part is seperated by Xcenter: kg/kg
        double ratio_log_to_linear;
        double Tmin, Tmax, pmin,pmax, Xmin, Xmax;
    };
    /**
    * @brief 为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。
    */
    struct CONSTENTS_Thermo
    {
        double R=461.51805, Tmin, Tmax, pmin, pmax, Ttriple, T_critical, p_critical, rhomass_critical, molar_mass, rhomolar_critical;
    };

#define STR_LENGTH_PROPINFO 30
    // some general data structure
    /**
     * @brief Information of a thermodynamic property
     *
     */
    struct propInfo
    {
        // string shortName, longName, unit;
        // int dimension[7]; //OpenFOAM dimension
        char shortName[STR_LENGTH_PROPINFO], longName[STR_LENGTH_PROPINFO], unit[STR_LENGTH_PROPINFO];
    };
    // ============= bitmask of properties =======================================
#define Update_prop_Rho     2   // 2^1 = 2
#define Update_prop_H       4   // 2^2 = 4
// #define Update_prop_drhodh  8   // 2^3 = 8
#define Update_prop_T       16   // 2^4 = 16
#define Update_prop_Cp       32   // 2^5 = 32
#define Update_prop_Mu       64   // 2^6 = 64
#define Update_prop_IsothermalCompressibility       128   // 2^7 = 128
#define Update_prop_IsobaricExpansivity       256   // 2^8 = 256
#define Update_prop_CommonlyUsedTPX (Update_prop_Rho | Update_prop_Cp | Update_prop_Mu | Update_prop_IsothermalCompressibility |Update_prop_IsobaricExpansivity )
#define Update_prop_CommonlyUsedHPX (Update_prop_H | Update_prop_T | Update_prop_Rho | Update_prop_Cp | Update_prop_Mu | Update_prop_IsothermalCompressibility |Update_prop_IsobaricExpansivity )
    // ===========================================================================
    /**
     * Head info of AMR_LUT which are used for Python API
     */
    struct Head_AMR_LUT
    {
        int dim=0;
        int space=-1;
        std::string spaceName="Unknown";
        int constWhich=-1;
        std::string constWhich_Name="Unknown";
        double constValue=0; //SI unit
        std::vector<std::string> var_names;
        std::vector<std::vector<double> > var_ranges;//SI unit
        std::vector<double> var_maxResolution; //SI unit
        int min_level, max_level;
        int num_leaves, num_quads, num_nodes, num_need_refine;
        int num_leaves_nextRefine;
        int num_props;
        std::vector<std::string> prop_names;
        std::string memory_total, memory_leaves, memory_nonLeaves, memory_quads, memory_properties;
    };
    // define some commonly used color
#define COLOR_gold {1.000000, 0.843137, 0.000000}
#define COLOR_crimson {0.862745, 0.078431, 0.235294}
#define COLOR_deeppink {1.000000, 0.078431, 0.576471}
#define COLOR_k {0.000000, 0.000000, 0.000000}
#define COLOR_yellow {1.000000, 1.000000, 0.000000}
#define COLOR_lightblue {0.678431, 0.847059, 0.901961}
#define COLOR_lime {0.000000, 1.000000, 0.000000}
#define COLOR_blue {0.000000, 0.000000, 1.000000}
#define COLOR_red {1.000000, 0.000000, 0.000000}
#define COLOR_lightgray {0.827451, 0.827451, 0.827451}
#define COLOR_gainsboro {0.862745, 0.862745, 0.862745}
#define COLOR_darkgray {0.662745, 0.662745, 0.662745}
#define COLOR_plum {0.866667, 0.627451, 0.866667}
#define COLOR_darkgray {0.662745, 0.662745, 0.662745}
#define COLOR_navajowhite {1.0, 0.8705882352941177, 0.6784313725490196}
#define COLOR_orange {1.0, 0.6470588235294118, 0.0}
#define COLOR_mfc_twoPhaseWater COLOR_gold
#define COLOR_mec_twoPhaseWater COLOR_crimson
#define COLOR_lc_twoPhaseWater COLOR_yellow
#define COLOR_mfc_CriticalPoint COLOR_deeppink
#define COLOR_lc_VH {1.0, 0.5137254901960784, 1.0}
#define COLOR_mec_CriticalPoint COLOR_k
#define COLOR_lc_VL_V {0.0, 0.6784313725490196, 0.6745098039215687}
#define COLOR_lc_VL_L {0.0, 0.6745098039215687, 1.0}
#define COLOR_lc_VLH {1.0, 0.7137254901960784, 0.37254901960784315}
#define COLOR_lc_LH COLOR_darkgray
#define COLOR_VL COLOR_plum
#define COLOR_VH COLOR_darkgray
#define COLOR_LH COLOR_lightblue
#define COLOR_VLH_lowT COLOR_navajowhite
#define COLOR_VLH_highT COLOR_yellow
#define COLOR_VLH {1.0, 0.7137254901960784, 0.37254901960784315}

    struct COLOR
    {
        double r,g,b;
    };
    struct Polygon_slice
    {
        std::vector<double> T,P,X,H;
        std::string name;
        std::vector<double> fc; //face color: [R,G,B]
        std::vector<double> ec; //edge color: [R,G,B]
        Polygon_slice(){}
        Polygon_slice(std::vector<double> T0, std::vector<double> P0, std::vector<double> X0, std::vector<double> H0, std::string name0, COLOR fc0, COLOR ec0)
                : T(std::move(T0)), P(std::move(P0)), X(std::move(X0)), H(std::move(H0)), name(std::move(name0)), fc{fc0.r, fc0.g, fc0.b}, ec{ec0.r, ec0.g, ec0.b}
        {

        }
        bool operator<(const Polygon_slice &other) const { return T < other.T; }
        bool operator==(const Polygon_slice &other) const { return name == other.name; }
        // std out
        friend std::ostream& operator<<(std::ostream& os, const Polygon_slice& self) {
            std::stringstream sin;
            sin<<"Polygon_slice {'T': ["<<self.T[0]<<",...,"<<self.T[self.T.size()-1]<<"],";
            sin<<"'P': ["<<self.P[0]<<",...,"<<self.P[self.P.size()-1]<<"],";
            sin<<"'X': ["<<self.X[0]<<",...,"<<self.X[self.X.size()-1]<<"],";
            sin<<"'H': ["<<self.H[0]<<",...,"<<self.H[self.H.size()-1]<<"],";
            sin<<"'name': '"<<self.name<<"',";
            sin<<"'fc': ("<<self.fc[0]<<","<<self.fc[1]<<","<<self.fc[2]<<"),";
            sin<<"'ec': ("<<self.ec[0]<<","<<self.ec[1]<<","<<self.ec[2]<<")";
            sin<<"}"<<std::endl;
            return os << sin.str();
        }
        // info for Python API
        std::string info(const Polygon_slice& self)
        {
            std::stringstream sin;
            sin<<self;
            return sin.str();
        }
    };
    struct Line_slice
    {
        std::vector<double>T,P,X,H;
        std::string name;
        std::vector<double> color;
        std::string ls;
        double lw;
        Line_slice(){}
        Line_slice(std::vector<double> T0, std::vector<double> P0, std::vector<double> X0, std::vector<double> H0, std::string name, COLOR c, std::string ls="solid", double lw=1)
                : T(std::move(T0)),P(std::move(P0)),X(std::move(X0)),H(std::move(H0)),name(std::move(name)),color{c.r, c.g, c.b}, ls(std::move(ls)),lw(lw)
        {
        }
        bool operator<(const Line_slice &other) const { return T < other.T; }
        bool operator==(const Line_slice &other) const { return name == other.name; }
        // std out
        friend std::ostream& operator<<(std::ostream& os, const Line_slice& self) {
            std::stringstream sin;
            sin<<"Line_slice {'T': ["<<self.T[0]<<",...,"<<self.T[self.T.size()-1]<<"],";
            sin<<"'P': ["<<self.P[0]<<",...,"<<self.P[self.P.size()-1]<<"],";
            sin<<"'X': ["<<self.X[0]<<",...,"<<self.X[self.X.size()-1]<<"],";
            sin<<"'H': ["<<self.H[0]<<",...,"<<self.H[self.H.size()-1]<<"],";
            sin<<"'name': '"<<self.name<<"',";
            sin<<"'ls': '"<<self.ls<<"',";
            sin<<"'color': ("<<self.color[0]<<","<<self.color[1]<<","<<self.color[2]<<"),";
            sin<<"'lw': "<<self.lw;
            sin<<"}"<<std::endl;
            return os << sin.str();
        }
        // info for Python API
        std::string info(const Line_slice& self)
        {
            std::stringstream sin;
            sin<<self;
            return sin.str();
        }
    };
    struct Point_slice
    {
        double T,P,X,H;
        std::string name;
        std::string marker;
        std::vector<double> mfc;
        std::vector<double> mec;
        // the construct function is needed for Python API
        Point_slice(){}
        Point_slice(double T, double P, double X, double H, std::string name, std::string marker, COLOR mfc0, COLOR mec0)
                : T(T),P(P), X(X), H(H), name(std::move(name)), marker(std::move(marker)), mfc{mfc0.r, mfc0.g, mfc0.b}, mec{mec0.r, mec0.g, mec0.b}
        {
        }
        bool operator<(const Point_slice &other) const { return T < other.T; }
        bool operator==(const Point_slice &other) const { return name == other.name; }
        // std out
        friend std::ostream& operator<<(std::ostream& os, const Point_slice& self) {
            std::stringstream sin;
            sin<<"Point_slice {'T': "<<self.T<<",";
            sin<<"'P': "<<self.P<<",";
            sin<<"'X': "<<self.X<<",";
            sin<<"'H': "<<self.H<<",";
            sin<<"'name': '"<<self.name<<"',";
            sin<<"'marker': '"<<self.marker<<"',";
            sin<<"'mfc': ("<<self.mfc[0]<<","<<self.mfc[1]<<","<<self.mfc[2]<<"),";
            sin<<"'mec': ("<<self.mec[0]<<","<<self.mec[1]<<","<<self.mec[2]<<")";
            sin<<"}"<<std::endl;
            return os << sin.str();
        }
        // info for Python API
        std::string info(const Point_slice& self)
        {
            std::stringstream sin;
            sin<<self;
            return sin.str();
        }
    };
    /**
     * Struct of phase region of a slice at constant T, P, H, or X. Contains phase region polygon, and some visualization properties.
     */
    struct PhaseRegion_Slice
    {
        std::map<std::string, std::vector<Polygon_slice> > regions;
        std::map<std::string, std::vector<Line_slice> > lines;
        std::map<std::string, std::vector<Point_slice> > points;
        // std out
        friend std::ostream& operator<<(std::ostream& os, const PhaseRegion_Slice& slice) {
            std::stringstream sin;
            sin<<"regions [";
            for(auto &item : slice.regions)sin<<"'"<<item.first<<"' ";
            sin<<"]"<<std::endl;
            sin<<"lines   [";
            for(auto &item : slice.lines)sin<<"'"<<item.first<<"' ";
            sin<<"]"<<std::endl;
            sin<<"points  [";
            for(auto &item : slice.points)sin<<"'"<<item.first<<"' ";
            sin<<"]"<<std::endl;
            return os << sin.str();
        }
        // info for Python API
        std::string info(const PhaseRegion_Slice& slice)
        {
            std::stringstream sin;
            sin<<slice;
            return sin.str();
        }
    };
};

#endif