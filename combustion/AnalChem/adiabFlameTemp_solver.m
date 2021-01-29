clc; clear;
% define constants
T_STP = 298;  % standard temp <K>
T_air = 475;  % <K>
T_ABS = T_STP;  % <K>

hf_air = 0;  % molar enthalpy of formation <kJ/kmol>
hf_ABS = 62628;  % molar enthalpy of formation <kJ/kmol>
hf_CO2 = -393520;  % molar enthalpy of formation <kJ/kmol>
hf_H2O = -241820;  % molar enthalpy of formation <kJ/kmol>
hf_N2 = 0;  % molar enthalpy of formation <kJ/kmol>

Cp_air = 29.67;  % molar specific heat  <kJ/kmol*K>
Cp_ABS = 29.67;  % molar specific heat  <kJ/kmol*K>
Cp_CO2 = 37.22;  % molar specific heat  <kJ/kmol*K>
Cp_H2O = 36.57;  % molar specific heat  <kJ/kmol*K>
Cp_N2 = 29.02;  % molar specific heat  <kJ/kmol*K>

CpMat = [Cp_air Cp_ABS Cp_CO2 Cp_H2O Cp_N2];
hfMat = [hf_air hf_ABS hf_CO2 hf_H2O hf_N2];
T_r = [T_air T_ABS T_STP T_STP T_STP]; 

% find adiabatic flame temperature for ABS air combustion
mol_r = [24.107 1 0 0 0];  % molar values of reactants
mol_p = [0 0 3.85 2.425 19.26];  % molar values of products

T_aft = (sum(mol_r.*hfMat) - sum(mol_p.*hfMat) + sum(mol_r.*CpMat.*(T_r-T_STP)) +...
    sum(mol_p.*CpMat.*T_STP)) / sum(mol_p.*CpMat);

fprintf("T_AFT = %.4f K",T_aft);

