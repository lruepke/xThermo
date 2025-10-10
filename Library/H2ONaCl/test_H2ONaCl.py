import numpy as np
from xThermal import H2O
from xThermal import NaCl
from xThermal import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")

sw = sw_84
props = sw.UpdateState_TPX_vector([373], [300E5], [0.2])
print(np.array(props))