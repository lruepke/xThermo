# import module
# IAPWS95 module
module IAPWS95Mod
using CxxWrap
@wrapmodule("/Users/zguo/MyData/Research/3-CodeProject/Hydrothermal-OpenFOAM/xThermal/install/API/julia/xThermal/H2O/libH2O_jl",:define_Julia_IAPWS95)
end
# # IF97 module
# module IF97Mod
# using CxxWrap
# @wrapmodule("H2O/libH2O_jl",:define_Julia_IF97)
# end

# Demo of using module
f = IAPWS95Mod.IAPWS95()
print(IAPWS95Mod.Rho(f, 373, 1E5))
props = IAPWS95Mod.UpdateState_TPX(f, 373, 1E5)
print(IAPWS95Mod.info(props))