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
classdef Chemistry
    properties
        fuel_type
        w
        x
        y
        z
        Cp_Fuel
        hf_Fuel
        MW_Fuel
        MW_Air
        MW_N2
        MW_O2
        MW_CO2
        T_STP
        T_air
        T_ABS
        CpMat
        hfMat
        T_mat
        mol_st
        mf_st
        ma_st
        f_st
        
    end
    methods
        function init = Chemistry(fuel_type)
            if ~exist('fuel_type','var')
                % fuel type defaulted to ABS
                init.fuel_type = 'ABS';
                init.w = 3.85 ;  % fuel Carbon coefficient
                init.x = 4.85;  % fuel Hydrogen coefficient
                init.y = 0;  % fuel Oxygen coefficient
                init.z = 0.43;  % fuel Nitrogen coefficient
                init.Cp_Fuel = 29.67;  % molar specific heat  <kJ/kmol*K>
                init.hf_Fuel = 62628;     % molar enthalpy of formation <kJ/kmol>
            end
            init.MW_Fuel = 12.011*init.w + 1.008*init.x + 15.9995*init.y + 14.007*init.z;                % Molecular Weight of Fuel
            init.MW_Air  = 28.9647;                                                  % Molecular Weight of Air
            init.MW_N2   = 14.007*2;                                                 % Molecular Weight of N2
            init.MW_O2   = 15.999*2;                                                 % Molecular Weight of O2
            init.MW_CO2  = 12.011 + 2*15.999;                                        % Molecular Weight of CO2
            init.T_STP = 298;    % standard temp <K>
            init.T_air = 475;    % <K>
            init.T_ABS = init.T_STP;  % <K>
            
            hf_air = 0;         % molar enthalpy of formation <kJ/kmol>
            hf_CO2 = -393520;   % molar enthalpy of formation <kJ/kmol>
            hf_H2O = -241820;   % molar enthalpy of formation <kJ/kmol>
            hf_N2  = 0;         % molar enthalpy of formation <kJ/kmol>
            hf_O2  = 0; 
            
            Cp_air = 29.67;     % molar specific heat  <kJ/kmol*K>
            Cp_CO2 = 37.22;     % molar specific heat  <kJ/kmol*K>
            Cp_H2O = 36.57;     % molar specific heat  <kJ/kmol*K>
            Cp_N2  = 29.02;     % molar specific heat  <kJ/kmol*K>
            Cp_O2  = 29.376;    % molar specific heat  <kJ/kmol*K> 
            
            init.CpMat = [init.Cp_Fuel Cp_air Cp_H2O Cp_CO2 Cp_N2 Cp_O2 init.Cp_Fuel]';
            init.hfMat = [init.hf_Fuel hf_air hf_H2O hf_CO2 hf_N2 hf_O2 init.hf_Fuel]';
            
            % runs stoichiometric coefficients
% ---------------------- Stoichiometric Combustion ---------------------- % 
% Cw Hx Ny + (A)(0.71N2 + 0.2102) -> (B)H2O + (D)CO2 + (E)N2
% 3.85 = D
% 4.85 = 2B
% 0.43 + A*2*0.71 = 2E
% A*0.21*2 = B + 2D
% ----------------------------------------------------------------------- % 
            %    a     b  c  d
            A = [0     0 -1  0;...  % carbon
                 0    -2  0  0;...  % hydrogen
                 0.42 -1 -2  0;...  % oxygen
                 1.58  0  0 -2];    % nitrogen
            B = [-init.w;-init.x;-init.y;-init.z];

            init.mol_st = linsolve(A,B); 
            init.mf_st = init.MW_Fuel;  % Mass of fuel (kg)
            init.ma_st = init.mol_st(1)*init.MW_Air;  % Mass of Air (kg)
            init.f_st = init.mf_st/init.ma_st;  % Stoichiometric fuel to air ratio            
            
        end
        
        function T_AFT = AFT(self, mol_mat, T_mat)
            T_AFT = (sum(mol_mat.*self.hfMat) + sum(mol_mat.*self.CpMat.*(T_mat-self.T_STP)))...
                / abs(sum(mol_mat(3:end).*self.CpMat(3:end)));
        end
        
        function [phi, T_AFT] = phiSolver(self, f, T_air)
            % ----------------- Initialize Values ----------------- % 
            f_yield = f;  % Actual fuel to air ratio     
            phi = f_yield/self.f_st;

            if(phi == 1)  
% ------------------------ Stoichiometric Combustion -------------------- % 
% C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
% ----------------------------------------------------------------------- %     
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e  f g
                A = [ 0 -1  0  0 0;...  % carbon
                     -2  0  0  0 0;...  % hydrogen
                     -1 -2  0  0 0;...  % oxygen
                      0  0 -2  0 0];    % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';


                T_AFT = AFT(self, mol_mat, self.T_mat); 
            %     T_AFT = 3189.37;                % Adiabatic Flame Temperature (K)
            %     gamma_nzl = 1.3845;             % gamma at nozzle throat

            elseif(phi < 1)                    
% ------------------------ Fuel lean Combustion ------------------------- % 
% C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
% ----------------------------------------------------------------------- % 
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e  f g
                A = [ 0 -1  0  0 0;...  % carbon
                     -2  0  0  0 0;...  % hydrogen
                     -1 -2  0 -2 0;...  % oxygen
                      0  0 -2  0 0];    % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';


                T_AFT = AFT(self, mol_mat, self.T_mat);    
            else

% ------------------------ Fuel Rich Combustion ------------------------- %
% C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
% ----------------------------------------------------------------------- % 
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e f  g
                A = [ 0 -1  0 0 -3.85;...  % carbon
                     -2  0  0 0 -4.85;...  % hydrogen
                     -1 -2  0 0  0   ;...  % oxygen
                      0  0 -2 0 -0.43];    % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';

                T_AFT = AFT(self, mol_mat, self.T_mat);     
            end 
        end
    end
end

