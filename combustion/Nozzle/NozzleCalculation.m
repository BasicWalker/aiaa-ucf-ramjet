%% yo what am i doing

clc; clear;

%[mach, T, P, rho, area] = flowisentropic(1.4, 2, 'mach')
%T_Tfinal = T/T0

gamma = 1.4;
mach_exit = 2;

%b = (gamma - 1)/(gamma + 1);

%T_Tfinal = (1 - b)/(1 + b(mach_exit^2 - 1));
T_Tfinal = (1 + ((gamma - 1)/2 * mach_exit^2))^-1;
P_Pfinal = (T_Tfinal)^(gamma/(gamma-1));
A_Astar_expo = (gamma + 1)/(2*(gamma - 1));
A_Astar = ((gamma + 1)/2)^(-A_Astar_expo) * (((1 + ((gamma - 1)/2)*(mach_exit^2))^A_Astar_expo)/mach_exit);

%%

% for mach_exit = 1.2:0.1:2.5
%     
%   T_Tfinal = (1 + ((gamma - 1)/2 * mach_exit^2))^-1;
%   P_Pfinal = (T_Tfinal)^(gamma/(gamma-1));
%   A_Astar_expo = (gamma + 1)/(2*(gamma - 1));
%   A_Astar = ((gamma + 1)/2)^(-A_Astar_expo) * (((1 + ((gamma - 1)/2)*(mach_exit^2))^A_Astar_expo)/mach_exit);
%   fprintf('At a Mach of %f:\n T/T0 is %f \n P/P0 is %f \n A/Astar is %f \n\n', mach_exit, T_Tfinal, P_Pfinal, A_Astar)
%   
% end
 %%
% 
%  for A_exit = 2:0.1:4
%      A_throat = A_exit/A_Astar;
%      fprintf('Throat area = %f \n', A_throat)
%  end
% 
%  
 %%
 clc
%  for P_ambient = 13.05:-0.1:11.53
%  P_chamber = P_ambient * (T_Tfinal)^(gamma/(gamma-1));
%  fprintf('The chamber pressure = %f\n', P_chamber)
%  
%  end
 P_atm = 13.05;
 syms P_chamber
 eqn = mach_exit^2 == (2/(gamma - 1))*(((P_chamber/P_atm)^((gamma - 1)/gamma)) - 1); 
 sol = vpasolve(eqn,P_chamber)
 
 