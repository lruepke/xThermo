

Todo List



#### Module [**PhysicalConstants\_H2ONaCl**](group__PhysicalConstants__H2ONaCl.md)  

Need to discuss how to deal with the discrepancy between the critical point of H2O in the original paper (Driesner2007Part1) and that in IAPWS95\_CoolProp release. 
Estimate density difference range based on the density difference between IAPWS95\_CoolProp and IAPS84.
How to deal with properties at high-T low-P region, e.g. sw.Rho\_phase(1.700000e+02+273.15,7.914706e+00\*1E5,0,Rhol, H2ONaCl::Liquid); see mmc4 of Driesner(2007b).
Merge Rho\_phase of liquid and vapor together to calculate them simultaneously, because property of saturated vapor and liquid always needed simultaneous, it could save computing time in this way. 
Confirm that \(\Delta H_{fus} (T_{melting}, P) = H_l(T_{melting}, P) - H_s(T_{melting}, P)\) ? 
Need to confirm r4 in Tab.5 is correct 
need to confirm the reference point of halite enthalpy. 




    


