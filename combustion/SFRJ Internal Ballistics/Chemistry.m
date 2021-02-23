% --------- AIAA Internal Ballistic Simulator code for UCF HPR ---------- %
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
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% ABS Chemical Formula: C3.85 H4.85 N0.43 
% Air Chemical Formula: 0.79N2 0.21O2
% Reaction:             H2O CO2 N2
% ----------------------------------------------------------------------- %
function [phi, gamma_nzl, T_AFT] = Chemistry(f)
% ----------------- Initialize Values ----------------- % 
f_yield = f;                                                        % Actual fuel to air ratio     
MW_ABS  = 12.011*3.85 + 1.008*4.85 + 14.007*0.43;                   % Molecular Weight of ABS
MW_Air  = 28.9647;                                                  % Molecular Weight of Air
MW_N2   = 14.007*2;                                                 % Molecular Weight of N2
MW_O2   = 15.999*2;                                                 % Molecular Weight of O2
MW_CO2  = 12.011 + 2*15.999;                                        % Molecular Weight of CO2

% ---------------------- Stoichiometric Combustion ---------------------- % 
% C3.85 H4.85 N0.43 + (A)(0.71N2 + 0.2102) -> (B)H2O + (D)CO2 + (E)N2
% 3.85 = D
% 4.85 = 2B
% 0.43 + A*2*0.71 = 2E
% A*0.21*2 = B + 2D
% ----------------------------------------------------------------------- % 
D = 3.85;                                                           % Systems of equations
B = 4.85/2;                                                         % Systems of equations
A = (B + 2*D)/(0.21*2);                                             % Systems of equations
E = (0.43 + A*2*0.71)/2;                                            % Systems of equations
mol_st = [A B D E];                                                 % Molar Matrix     
mf_st = MW_ABS;                                                     % Mass of ABS (kg)
ma_st = mol_st(1)*MW_Air;                                           % Mass of Air (kg)
f_st = mf_st/ma_st;                                                 % Stoichiometric fuel to air ratio
%fprintf('\nCHN + %.2f(0.71N2 + 0.21O2) -> %.2fH2O + %.2fCO2 + %.2fN2\n', mol_st(1),mol_st(2),mol_st(3),mol_st(4));

phi = f_yield/f_st;

if(phi == 1)                            %  Fuel lean Combustion 
        
        T_AFT = 3189.37;                % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;             % gamma at nozzle throat
        
    elseif(phi < 1)                    
% ------------------------ Fuel lean Combustion ------------------------- % 
% C3.85 H4.85 N0.43 + (1/phi)(A)(0.71N2 + 0.2102) -> (B)H2O + (D)CO2 + (E)N2 + (F)O2
% C: 3.85 = D
% H: 4.85 = 2B
% N: 0.43 + A*2*0.71 = 2E
% O: A*0.21*2 = B + 2D
% ----------------------------------------------------------------------- % 
        
        T_AFT = 3189.37;            % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;         % gamma at nozzle throat
    
    else                            
% ------------------------ Fuel Rich Combustion ------------------------- % 
% C3.85 H4.85 N0.43 + (1/phi)(A)(0.71N2 + 0.2102) -> (B)H2O + (D)CO2 + (E)N2 + (F)C3.85 H4.85 N0.43
% C: 3.85 = D + 3.85F
% H: 4.85 = 2B + 4.85F
% N: 0.43 + 2*0.79*A*(1/phi) = 2E + 0.43F
% O: 2*0.21*A*(1/phi) = B + 2D
% ----------------------------------------------------------------------- % 
molair = mol_st(1);
% syms B D E F
% eqn1 = D + 3.85*F == 3.85;
% eqn2 = 3*B + 4.85*F == 4.85;
% eqn3 = 2*E + 0.43*F == 0.43 + (2*0.79*molair)/phi;
% eqn4 = B + 2*D == (2*0.21*molair)/phi;
% [LeftMatrix,RightMatrix] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4],[B, D, E, F]);
% mol_fr = vpa(linsolve(LeftMatrix,RightMatrix));

        
        T_AFT = 3189.37;            % Adiabatic Flame Temperature (K)
        gamma_nzl = 1.3845;         % gamma at nozzle throat        
        
    end 

% Calculate heat of combustion

    reacttemp = 475; %K
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
    mcpabs = cpabs*mwabs;
    mcpair = cpair*mwair;
    mcph2o = cph2o*mwh2o;
    mcpco2 = cpco2*mwco2;
    mcpn2 = cpn2*mwn2;

    %Mole Fractions
    mfabs = molf/molr;
    mfair = molair/molr;
    mfh2o = molh2o/molp;
    mfco2 = molco2/molp;
    mfn2 = moln2/molp;

    %Average Molar Cp values
    mcpr = (mfabs * mcpabs) + (mfair * mcpair); %Reactants
    mcpp = (mfh2o * mcph2o) + (mfco2 * mcpco2) + (mfn2 * mcpn2);

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

    molhrxn = h2oterm + co2term - absterm + (molp * mcpp * (reacttemp - tstp)) - (molr * mcpr * (reacttemp - tstp));

    hrxn = molhrxn/mwabs;

    hoc = hrxn * -1;

    %fprintf("Heat of Combustion: %f kJ/kg\n", hoc);
        

end
