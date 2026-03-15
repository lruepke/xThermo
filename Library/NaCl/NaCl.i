%module NaCl
%{
    #include "thermo.h"
    #include "NaCl.h"
    using namespace xThermal;
%}

%include std_string.i
%include std_vector.i
// ===== This is required for python API========
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<string>;
   %template(ConstCharVector) vector<const char*>;
   %template(UnsignedLongVector) vector<unsigned long int>;
}
// ========================================================================
%inline %{ 
    using namespace std;
%}
%include "typemaps.i"
%apply std::vector<double> *OUTPUT {std::vector<double>& res};

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
struct ThermodynamicProperties
{
    double T=0, p=0, X=0, H=0;
    PhaseRegion phase=Unknown;
    double S_l=0, S_v=0, S_h=0;
    double X_l=0, X_v=0;
    double Rho_l=0, Rho_v=0, Rho_h=0;
    double H_l=0, H_v=0, H_h=0;
    double Cp_l = 0, Cp_v = 0, Cp_h = 0;
    double Mu_l = 1, Mu_v = 1;
    double Rho=0, Cp = 0, Mu = 1;
    // derivative
    double dRhodP = 0, dRhodP_l = 0, dRhodP_v = 0;
    double dRhodT = 0, dRhodT_l = 0, dRhodT_v = 0;
    // commonly used derivatives for hydrothermal convection simulation
    double IsothermalCompressibility = 0, IsothermalCompressibility_l =0, IsothermalCompressibility_v =0;
    double IsobaricExpansivity=0, IsobaricExpansivity_l =0, IsobaricExpansivity_v =0;
    std::string fluidName;
    std::string info(const ThermodynamicProperties& self);
};
%extend ThermodynamicProperties {
        std::string __str__() {
            return self->info(*self);
        }
};
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
    std::vector<double> dRhodP, dRhodP_l, dRhodP_v;
    std::vector<double> dRhodT, dRhodT_l, dRhodT_v;
    std::vector<double> IsothermalCompressibility, IsothermalCompressibility_l, IsothermalCompressibility_v;
    std::vector<double> IsobaricExpansivity, IsobaricExpansivity_l, IsobaricExpansivity_v;
    std::string fluidName;
    std::string info(const ThermodynamicPropertiesVector& self);
};
%extend ThermodynamicPropertiesVector {
        std::string __str__() {
            return self->info(*self);
        }
};

//%apply double *OUTPUT {double& rho_l, double& rho_v};

//ignore some functions
//%ignore Boiling_p(const double& T, double&, double&);
//%ignore Boiling_T(const double& p, double&, double&);
%ignore UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);

//rename overloaded function
%rename(Boiling_p_rho) Boiling_p(const double& T, double&, double&);
%rename(Boiling_T_rho) Boiling_T(const double& p, double&, double&);

//ignore some warnings
%warnfilter(560,562) ;

#define xTHERMO_VAR
%include "thermo.h"
%include "NaCl.h"
