import numpy as np
from xThermal import H2O
from xThermal import NaCl
from xThermal import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")

sw = sw_84
T_list = [373, 500, 650]
p_list = [300e5, 300e5, 300e5]
X_list = [0.2, 0.2, 0.2]
for T, p, X in zip(T_list, p_list, X_list):
    props = sw.UpdateState_TPX(T, p, X)
    print(f"T={T} K  p={p/1e5:.0f} bar  X={X}  Rho={props.Rho:.2f} kg/m3  H={props.H/1e3:.2f} kJ/kg  phase={props.phase}")