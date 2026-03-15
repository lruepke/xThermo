

# File IAPWS95\_CoolProp.h

[**File List**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**IAPWS95\_CoolProp.h**](IAPWS95__CoolProp_8h.md)

[Go to the documentation of this file](IAPWS95__CoolProp_8h.md)


```C++

#ifdef USE_COOLPROP

#ifndef IAPWS95_CoolProp_xThermal_H
#define IAPWS95_CoolProp_xThermal_H

#include "thermo.h"

namespace xThermal
{
    // IAPWS95_CoolProp
    #define Name_Backend_IAPWS95_CoolProp "IAPWS95_CoolProp"
    namespace COOLPROP
    {

        class cIAPWS95_CoolProp : public cxThermal
        {
        private:
            CONSTENTS_Thermo m_constants;
        public:
            cIAPWS95_CoolProp(/* args */);
            ~cIAPWS95_CoolProp();

        private:
            void initialize_data();
            PhaseRegion phase_CoolProp2xThermal(const int& phase_);
        public:
            std::string name(){return Name_Backend_IAPWS95_CoolProp;};
            // thermodynamic constants of the model
            double Tmin(){return m_constants.Tmin; };                  
            double Tmax(){return m_constants.Tmax;};                  
            double pmin(){return m_constants.pmin; };                    
            double pmax(){return m_constants.pmax;};                  
            double Ttriple(){return m_constants.Ttriple;};               
            double T_critical(){return m_constants.T_critical;};            
            double p_critical(){return m_constants.p_critical;};            
            double rhomass_critical(){return m_constants.rhomass_critical;};      
            double molar_mass(){return m_constants.molar_mass;};       
            // update thermodynamic state and properties for given [T,p], or [T,p,X], X is the salinity of H2ONaCl. For pure fluid, X=0 as a default.
            PhaseRegion findPhaseRegion_TPX(const double& T, const double& p, const double& X=0);
            void UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);
            void UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X=0);
            ThermodynamicProperties UpdateState_TPX(const double& T, const double& p, const double& X=0); //for Python API
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

#endif
```


