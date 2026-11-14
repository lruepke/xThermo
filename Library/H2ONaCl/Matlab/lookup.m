%LOOKUP Get thermodynamic properties from the input AMR-LUT(Adaptive mesh refined lookup table) for given T|H, P, X.
%
%   % Example: lookup properties of H2O-NaCl in T-P-X space
%   LOOKUP(filename_AMR_LUT, 'H2O-NaCl', T, P, X, backend_name='IAPS84');
%   - Input: filename_AMR_LUT (string), fluid_name (string), T (double) [K], P (double) [Pa], X (double) [kg/kg], backend_name (string: optional, default is 'IAPS84'. [IAPS84, IAPWS95, IAPWS95_CoolProp])
%   - Output: Properties struct. Properties are stored in each field of the struct.
%
%   % Example: lookup properties of H2O-NaCl in H-P-X space
%   LOOKUP(filename_AMR_LUT, 'H2O-NaCl', H, P, X, backend_name='IAPS84');
%   - Input: filename_AMR_LUT (string), fluid_name (string), H (double) [J/kg], P (double) [Pa], X (double) [kg/kg], backend_name (string: optional, default is 'IAPS84'. [IAPS84, IAPWS95, IAPWS95_CoolProp])
%   - Output: Properties struct. Properties are stored in each field of the struct.
%
%   % Example: lookup properties of H2O in T-P space
%   LOOKUP(filename_AMR_LUT, 'H2O', T, P, backend_name='IAPS84');
%   - Input: filename_AMR_LUT (string), fluid_name (string), T (double) [K], P (double) [Pa], backend_name (string: optional, default is 'IAPS84'. [IAPS84, IAPWS95, IAPWS95_CoolProp])
%   - Output: Properties struct. Properties are stored in each field of the struct.
%
%   % Example: lookup properties of H2O in H-P space
%   LOOKUP(filename_AMR_LUT, 'H2O', H, P, backend_name='IAPS84');
%   - Input: filename_AMR_LUT (string), fluid_name (string), H (double) [J/kg], P (double) [Pa], backend_name (string: optional, default is 'IAPS84'. [IAPS84, IAPWS95, IAPWS95_CoolProp])
%   - Output: Properties struct. Properties are stored in each field of the struct.
%
%   See also prop_HPX, prop_TPX.
%
% See <a href="lookup.html">online documentation for more information</a>.