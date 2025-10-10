//https://docs.juliahub.com/CxxWrap/WGIJU/0.13.0/

#include "IAPWS95.h"
#include "IAPWS-IF97.h"

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

using namespace xThermal::IAPWS95;
using namespace xThermal::IAPWS_IF97;

JLCXX_MODULE define_Julia_IAPWS95(jlcxx::Module& mod)
{
    using namespace jlcxx;
    mod.add_type<xThermal::ThermodynamicProperties>("ThermodynamicProperties")
    .constructor()
    .method("T", [](const xThermal::ThermodynamicProperties& prop) { return prop.T; });

    mod.method("info", [](const xThermal::ThermodynamicProperties& prop) { std::stringstream sin;
            sin<<prop;
            return sin.str(); })
    ;

    mod.add_type<cIAPWS95>("IAPWS95")
        .constructor()
        .method("Rho", [](cIAPWS95& instance, double T, double p) {
                return instance.Rho(T, p, "bisection");
            })
        .method("UpdateState_TPX", [](cIAPWS95& instance, double T, double p) {
                // std::cout<<instance.UpdateState_TPX(T, p)<<"\n"<<endl;
                
                return instance.UpdateState_TPX(T, p);
            })
        ;
}
