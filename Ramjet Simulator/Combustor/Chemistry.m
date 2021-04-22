% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Chemistry.m 
% 
% File Description: 
% Chemistry model, calculates thermochemical reactions and stochiometric
% combustion
% ABS Chemical Formula: C3.85 H4.85 N0.43 
% Air Chemical Formula: 0.79N2 0.21O2
% Products:             H2O CO2 N2 (ABS) (02)
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  002  Chemical balance, T_AFT, gamma calculations
% Samer Armaly    03/12/21  002  Optimization Update
% Ethan Sherlock  04/02/21  ---  Added gamma/R solver
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
        MW_H2O
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
% ---------------------------------------------------------------------- %
% Function name: init
%
% Function Description: 
% Initialize Chemistry Model
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Samer Armaly    03/12/21  002  Optimization update
% ---------------------------------------------------------------------- %
        
        function init = Chemistry(fuel_type)
            if ~exist('fuel_type','var')
                % fuel type defaulted to ABS
                init.fuel_type = 'ABS';
                init.w = 3.85 ;                 % fuel Carbon coefficient
                init.x = 4.85;                  % fuel Hydrogen coefficient
                init.y = 0;                     % fuel Oxygen coefficient
                init.z = 0.43;                  % fuel Nitrogen coefficient
                init.Cp_Fuel = 81.33;           % molar specific heat  <kJ/kmol*K>
                init.hf_Fuel = 62628;           % molar enthalpy of formation <kJ/kmol>
            end
            init.MW_Fuel = 12.011*init.w + 1.008*init.x + 15.9995*init.y + 14.007*init.z;   % Molecular Weight of Fuel
            init.MW_Air  = 28.9647;             % Molecular Weight of Air
            init.MW_N2   = 14.007*2;            % Molecular Weight of N2
            init.MW_O2   = 15.999*2;            % Molecular Weight of O2
            init.MW_CO2  = 12.011 + 2*15.999;   % Molecular Weight of CO2
            init.MW_H2O  = 18.02;               % Molecular Weight of H2O
            init.T_STP = 298;                   % standard temp <K>
            init.T_air = 475;                   % <K>
            init.T_ABS = init.T_STP;            % <K>
            
            hf_air = 0;                         % molar enthalpy of formation <kJ/kmol>
            hf_CO2 = -393520;                   % molar enthalpy of formation <kJ/kmol>
            hf_H2O = -241820;                   % molar enthalpy of formation <kJ/kmol>
            hf_N2  = 0;                         % molar enthalpy of formation <kJ/kmol>
            hf_O2  = 0; 
            
            Cp_air = 29.67;                     % molar specific heat  <kJ/kmol*K>
            Cp_CO2 = 37.22;                     % molar specific heat  <kJ/kmol*K>
            Cp_H2O = 36.57;                     % molar specific heat  <kJ/kmol*K>
            Cp_N2  = 29.02;                     % molar specific heat  <kJ/kmol*K>
            Cp_O2  = 29.376;                    % molar specific heat  <kJ/kmol*K> 
            
            init.CpMat = [init.Cp_Fuel Cp_air Cp_H2O Cp_CO2 Cp_N2 Cp_O2 init.Cp_Fuel]';
            init.hfMat = [init.hf_Fuel hf_air hf_H2O hf_CO2 hf_N2 hf_O2 init.hf_Fuel]';
            
            % Stoichiometric Combustion 
            % Cw Hx Ny + (A)(0.71N2 + 0.2102) -> (B)H2O + (D)CO2 + (E)N2
            % 3.85 = D
            % 4.85 = 2B
            % 0.43 + A*2*0.71 = 2E
            % A*0.21*2 = B + 2D            
            %    a     b  c  d
            A = [0     0 -1  0;...              % carbon
                 0    -2  0  0;...              % hydrogen
                 0.42 -1 -2  0;...              % oxygen
                 1.58  0  0 -2];                % nitrogen
            B = [-init.w;-init.x;-init.y;-init.z];

            init.mol_st = linsolve(A,B); 
            init.mf_st = init.MW_Fuel;                  % Mass of fuel (kg)
            init.ma_st = init.mol_st(1)*init.MW_Air;    % Mass of Air (kg)
            init.f_st = init.mf_st/init.ma_st;          % Stoichiometric fuel to air ratio            
            
        end
% ---------------------------------------------------------------------- %
% Function name: T_AFT
%
% Function Description: 
% Adiabaitc flame temp calculator
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Samer Armaly    03/12/21  002  Optimization update
% ---------------------------------------------------------------------- %        
        function T_AFT = AFT(self, mol_mat, T_mat)
            T_AFT = (sum(mol_mat.*self.hfMat) + sum(mol_mat.*self.CpMat.*(T_mat-self.T_STP)))...
                / abs(sum(mol_mat(3:end).*self.CpMat(3:end)));
        end
% ---------------------------------------------------------------------- %
% Function name: gammaRSolver
%
% Function Description: 
% Gamma/R solver
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Nathan Starks   04/01/21  ---  Initial Creation
% Ethan Sherlock  04/02/21  ---  Update for fuel lean/rich conditions
% ---------------------------------------------------------------------- %  
        function [gamma, R] = gammaRSolver(self, mol_mat, T_AFT)
            Ru = 8.314;                                                     % kJ/kmol*k
               
            % Molar Cp as a function of temperature, see Cp calculation chart 
            theta = T_AFT/100;
            Cp_H2O = 143.05 - (183.54 * (theta^0.25)) + (82.751 * (theta^ 0.5)) - (3.6989 * theta);
            Cp_CO2 = -3.7357 + (30.529 * (theta^0.5)) - (4.1034 * theta) + (.024198 * (theta^2));
            Cp_N2 = 39.06 - (512.79 * (theta^(-1.5))) + (1072.7 * (theta^(-2))) - (820.4 * (theta^(-3)));
            Cp_O2 = 37.432 + (0.02102 *(theta^1.5)) - (178.57*(theta^-1.5)) + (236.88 * (theta^-2));
            
            % Molar basis to mass basis
            Cp_H2O_kg = Cp_H2O/self.MW_Air;                                 % kJ/Kgk
            Cp_CO2_kg = Cp_CO2/self.MW_CO2;                                 % kJ/Kgk
            Cp_N2_kg = Cp_N2/self.MW_N2;                                    % kJ/Kgk
            Cp_O2_kg = Cp_O2/self.MW_O2;                                    % kJ/Kgk
            Cp_ABS_kg = 1.423;                                              % kJ/Kgk
            
            % Convert products from moles to Kg
            mass_H2O = -1*mol_mat(3) * self.MW_H2O;                         % Mass of H2O
            mass_CO2 = -1*mol_mat(4) * self.MW_CO2;                         % Mass of CO2
            mass_N2 = -1*mol_mat(5) * self.MW_N2;                           % Mass of N2
            mass_O2 = -1*0.0 * self.MW_O2;                                  % Mass of O2
            mass_ABS = -1*0.0 * self.MW_Fuel;                               % Mass of ABS
            
            mass_p = mass_H2O + mass_CO2 + mass_N2 + mass_O2 + mass_ABS;    % Mass of Products
            
            % Mass Fraction of Products
            maf_H2O = mass_H2O/mass_p;                                      % Mass Fraction of H2O
            maf_CO2 = mass_CO2/mass_p;                                      % Mass Fraction of CO2
            maf_N2 = mass_N2/mass_p;                                        % Mass Fraction of N2
            maf_O2 = mass_O2/mass_p;                                        % Mass Fraction of 02
            maf_ABS = mass_ABS/mass_p;                                      % Mass Fraction of ABS
            
            % Sum mass fraction Cp Value of products
            m_Cpp = (maf_H2O*Cp_H2O_kg)+(maf_CO2*Cp_CO2_kg)+(maf_N2*Cp_N2_kg) + (maf_O2*Cp_O2_kg) + (maf_ABS*mass_ABS);
            
            % Gas Constant (R)
            R_H2O = Ru/self.MW_H2O;                                         % kJ/Kg*K
            R_CO2 = Ru/self.MW_CO2;                                         % kJ/Kg*K
            R_N2  = Ru/self.MW_N2;                                          % kJ/Kg*K
            R_O2  = Ru/self.MW_O2;                                          % kJ/Kg*K
            R_ABS = Ru/self.MW_Fuel;                                        % kJ/Kg*K
            
            % Cv values from Cp and R values
            Cv_H2O = Cp_H2O_kg - R_H2O;
            Cv_CO2 = Cp_CO2_kg - R_CO2;
            Cv_N2  = Cp_N2_kg - R_N2;
            Cv_O2  = Cp_O2_kg - R_O2;
            Cv_ABS = Cp_ABS_kg - R_ABS;
            
            % Sum mass fraction Cv value of products
            m_Cvp = (maf_H2O*Cv_H2O) + (maf_CO2*Cv_CO2) + (maf_N2*Cv_N2) + (maf_O2*Cv_O2) + (maf_ABS*Cv_ABS);
            
            % Gas constant (R) of products
            R = ((maf_H2O*R_H2O) + (maf_CO2*R_CO2) + (maf_N2*R_N2) + (maf_O2*R_O2) + (maf_ABS*R_ABS))*1e3;  % j/(kg*K)
            
            % Specific Heat Ratio
            gamma = m_Cpp/m_Cvp;
            
        end
% ---------------------------------------------------------------------- %
% Function name: phiSolver
%
% Function Description: 
% Equivalence ratio solver
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Samer Armaly    --------  ---  Initial Creation 
% ---------------------------------------------------------------------- %   
        function [phi, T_AFT, gamma, R] = phiSolver(self, f, T_air)
            % ----------------- Initialize Values ----------------- % 
            f_yield = f;                        % Actual fuel to air ratio     
            phi = f_yield/self.f_st;

            if(phi == 1)  
                
                % Stoichiometric Combustion
                % C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
                
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e  f g
                A = [ 0 -1  0  0 0;...          % carbon
                     -2  0  0  0 0;...          % hydrogen
                     -1 -2  0  0 0;...          % oxygen
                      0  0 -2  0 0];            % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';


                T_AFT = AFT(self, mol_mat, self.T_mat); 

                [gamma, R] = gammaRSolver(self, mol_mat, T_AFT);
            elseif(phi < 1)                    
                
                % Fuel lean Combustion
                % C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
                
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e  f g
                A = [ 0 -1  0  0 0;...          % carbon
                     -2  0  0  0 0;...          % hydrogen
                     -1 -2  0 -2 0;...          % oxygen
                      0  0 -2  0 0];            % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';


                T_AFT = AFT(self, mol_mat, self.T_mat);    
                
                [gamma, R] = gammaRSolver(self, mol_mat, T_AFT);
            else

                % Fuel Rich Combustion
                % C3.85 H4.85 N0.43 + (a)(0.71N2 + 0.2102) -> (b)H2O + (d)CO2 + (e)N2 + (f)O2 + (g)C3.85 H4.85 N0.43
                
                air_coeff = self.mol_st(1)/phi;
                %     b  d  e f  g
                A = [ 0 -1  0 0 -3.85;...       % carbon
                     -2  0  0 0 -4.85;...       % hydrogen
                     -1 -2  0 0  0   ;...       % oxygen
                      0  0 -2 0 -0.43];         % nitrogen
                B = [-self.w;-self.x;(-self.y - 0.42*air_coeff);(-self.z - 1.58*air_coeff)];

                mol_mat = linsolve(A,B);
                mol_mat = [1; air_coeff; mol_mat*-1];
                self.T_mat = [self.T_STP T_air 0 0 0 0 0]';

                T_AFT = AFT(self, mol_mat, self.T_mat);     
                
                [gamma, R] = gammaRSolver(self, mol_mat, T_AFT);
            end 
        end        
    end
end

