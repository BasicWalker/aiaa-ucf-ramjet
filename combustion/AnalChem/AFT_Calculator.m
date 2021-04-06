%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name:  Adiabatic Flame Temperature Calculations
% 
% Program Description: 
%
% 
% File Name: AFT_Calculator.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Karam Paul      01/17/21  Initial Creation 
% --------------------------------------------------------------------- %

clc; clear;

% Define Constants
T_STP = 298;    % standard temp <K>
T_air = 475;    % <K>
T_ABS = T_STP;  % <K>

hf_ABS = 62628;     % molar enthalpy of formation <kJ/kmol>
hf_air = 0;         % molar enthalpy of formation <kJ/kmol>
hf_CO2 = -393520;   % molar enthalpy of formation <kJ/kmol>
hf_H2O = -241820;   % molar enthalpy of formation <kJ/kmol>
hf_N2  = 0;          % molar enthalpy of formation <kJ/kmol>
hf_O2  = 0; 

Cp_ABS = 29.67;  % molar specific heat  <kJ/kmol*K>
Cp_air = 29.67;  % molar specific heat  <kJ/kmol*K>
Cp_CO2 = 37.22;  % molar specific heat  <kJ/kmol*K>
Cp_H2O = 36.57;  % molar specific heat  <kJ/kmol*K>
Cp_N2  = 29.02;   % molar specific heat  <kJ/kmol*K>
Cp_O2  = 29.376;  % molar specific heat  <kJ/kmol*K> 

CpMat = [Cp_ABS Cp_air  Cp_CO2 Cp_H2O Cp_N2 Cp_O2];
hfMat = [hf_ABS hf_air  hf_CO2 hf_H2O hf_N2 hf_O2];
T_r   = [T_ABS T_air T_STP T_STP T_STP T_STP]; 


%% Balance chemical equations based on equivalence ratios

phi(1,1) = 1.0; 
mol_r = [1 24.107 0 0 0 0];                                         % molar values of reactants
mol_p = [0 0 3.85 2.425 19.26 0];                                   % molar values of products
AFT(1,1) = aft_equation(mol_r, mol_p, hfMat,CpMat, T_r, T_STP);

phi(2,1) = 0.75;
mol_r = [1 36.1605 0 0 0 0];                                        % molar values of reactants
mol_p = [0 0 3.85 2.425 28.782 2.531];                              % molar values of products
AFT(2,1) = aft_equation(mol_r, mol_p, hfMat,CpMat, T_r, T_STP);

phi(3,1) = 0.50;
mol_r = [1 48.214 0 0 0 0];                                         % molar values of reactants
mol_p = [0 0 3.85 2.425 38.3 5.062];                                % molar values of products
AFT(3,1) = aft_equation(mol_r, mol_p, hfMat,CpMat, T_r, T_STP);

fprintf("T_AFT = %.4f K \n",AFT);

T = table(phi, AFT);



