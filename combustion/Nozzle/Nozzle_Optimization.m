% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Nozzle_Optimization.m 
% 
% File Description: 
% Numerical solver to design the optimal nozzle geometry
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  03/28/21  ---  Initial Creation 
% ---------------------------------------------------------------------- %
% [mach,T,P,rho,area] = flowisentropic(gamma_nzlT(1),2,'sub') % Isentropic Relations

% ----------------------------- User Inputs ---------------------------- %
P_c = 400.0;                        % Pressure combustion chamber
T_c = 3000.0;                       % Temperature combustion chamber
mach_exit = 2.0;                    % Exit Mach 
m_dot = 1.3;                        % Total mass flow rate
g = 1.3;                            % Gamma at the throat

% ------------------------ Initialized Variables ----------------------- %
GM1 = g - 1;                        % Gamma + 1     
GP1 = g + 1;                        % Gamma - 1 
MW = 45.3;                          % Molecular weight of products
R_u = 8314.4621;                    % Universal gas constant (J/kmol-k)
R = R_u/MW;                         % Gas constant 

% ---------------------------- Calculations ---------------------------- %
P_t = P_c * (1 + GM1/2)^(-g/GM1);   % Pressure at Nozzle Throat
T_t = T_c / (1 + GM1/2);            % Temperature at Nozzle Throat

A_t = (M_dot/P_t) * sqrt((R*T_t)/g);

fprintf('-------------- Nozzle Results --------------\n')
fprintf('Area of Throat: %.2f  \n', A_t)
fprintf('--------------------------------------------\n')