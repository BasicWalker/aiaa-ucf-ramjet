% --------- AIAA Internal Ballistic Simulator code for UCF HPR --------- %
% File Name: Chemistry.m 
% 
% File Description: 
% Chemistry model, calculates thermochemical reactions and stochiometric
% combustion
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  002  Chemical balance, T_AFT, gamma calculations
% ---------------------------------------------------------------------- %

% ---------------------------------------------------------------------- %
% ABS Chemical Formula: C8H8 C4H6 C3H3N 
% Air Chemical Formula: 0.79N2 0.21O2
% Reaction:             H2O CO2 N2
% ---------------------------------------------------------------------- %
function [gamma_nzl, T_AFT] = Chemistry(f)

f_yield = f;            

% Stoich combustion chemical balance equation
f_st = 0.079;

% Equivalence ratio
phi = f_yield/f_st;                 % Equivalence Ratio



    if(phi == 1)                      % Equivalence Ratio = 1    
        
        T_AFT = 3189.37;            % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;         % gamma at nozzle throat
        
    elseif(phi < 1)                   % Equivalence Ratio < 1, fuel lean
        
        T_AFT = 3189.37;            % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;         % gamma at nozzle throat
    
    else                            % Equivalence Ratio > 1, fuel rich
        
        T_AFT = 3189.37;            % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;         % gamma at nozzle throat        
        
    end 

end
