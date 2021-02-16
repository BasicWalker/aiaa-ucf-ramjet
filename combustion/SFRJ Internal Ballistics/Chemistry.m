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

% Calculate heat of combustion

    reacttemp = 475 %K
    molf = 1; %Moles of fuel (stoichiometric)
    molair = 24.107;
    molh2o = 2.425;
    molco2 = 3.85;
    moln2 = 19.26;
    
    molr = molf + molair;
    molp = molh2o + molco2 + moln2;
    
    %Molecular weights needed (kg/kmol)
    %Reactants
    mwabs = 57.07;
    mwair = 28.96;
    
    %Products
    mwh2o = 18.02;
    mwco2 = 44.01;
    mwn2 = 28.013;
    
    %HEAT OF FORMATION FOR EACH MOLECULE%

    %hf(molecule) = heat of formation for each molecule
    %Values found online

    hfn2 = 0;
    hfo2 = 0;
    hfco2 = -393520; %kJ/kmol
    hfh2o = -241820; %kJ/kmol
    hfabs = 62630; %kJ/kmol

    tstp = 273.15; %K

    %Molar Cp of molecules needed (CP VALUES CHANGE WITH TEMPERATURE,WILL BE INCORPORATED AT A LATER TIME)

    %Mass Cp values (kJ/kg*K)
    %Other values taken at 300K
    cpabs = 1.423; %Average value given, won't change with temperature
    cpair = 1.006;
    cph2o = 1.864;
    cpco2 = .846;
    cpn2 = 1.04;

    %Molar Cp values
    mcpabs = cpabs*mwabs
    mcpair = cpair*mwair
    mcph2o = cph2o*mwh2o
    mcpco2 = cpco2*mwco2
    mcpn2 = cpn2*mwn2

    %Mole Fractions
    mfabs = molf/molr;
    mfair = molair/molr;
    mfh2o = molh2o/molp;
    mfco2 = molco2/molp;
    mfn2 = moln2/molp;

    %Average Molar Cp values
    mcpr = (mfabs * mcpabs) + (mfair * mcpair); %Reactants
    mcpp = (mfh2o * mcph2o) + (mfco2 * mcpco2) + (mfn2 * mcpn2)

    %HEAT OF COMBUSTION*

    %Long equation will be simplified into terms to eliminate chances of parentheses errors.

    %Followiing terms are number of moles multiplied by heat of formation
    %HEATS OF FORMATION FOR N2 AND O2 ARE ZERO

    %Products
    h2oterm = molh2o * hfh2o;
    co2term = molco2 * hfco2;


    %Reactants
    absterm = molf * hfabs;

    %Heat of reaction determined, then multiplied by -1 to get heat of combustion.

    molhrxn = h2oterm + co2term - absterm + (molp * mcpp * (reacttemp - tstp)) - (molr * mcpr * (reacttemp - tstp))

    hrxn = molhrxn/mwabs;

    hoc = hrxn * -1;

    fprintf("Heat of Combustion: %f kJ/kg\n", hoc)
        

end
