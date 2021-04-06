% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Nozzle.m 
% 
% File Description: 
% Numerical solver to design the optimal nozzle geometry
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %
% [mach,T,P,rho,area] = flowisentropic(gamma_nzlT(1),2,'sub') % Isentropic Relations

% Pres_throat = min()
P_c = 400.0;                        % Pressure combustion chamber
T_c = 3000.0;                       % Temperature combustion chamber
g = 1.3;                            % Gamma at the throat
GM1 = g - 1;                        % Gamma + 1     
GP1 = g + 1;                        % Gamma - 1 
M_dot = 1.3;                        % Mass flow rate total 
R_u = 8314.4621;                    % Universal gas constant (J/kmol-k)
%MW =            % Molecular weight of products


P_t = P_c * (1 + GM1/2)^(-g/GM1);   % Pressure at Nozzle Throat
T_t = T_c / (1 + GM1/2);            % Temperature at Nozzle Throat

% A_t = (M_dot/P_t) * sqrt()