%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           AIAA UCF Ramjet Thrust pressure trade study Script            %
%                                                                         %
%                              Authored by                                %
%               Samer Armaly, Karam Paul, Matthew Aubertin                %
%                           February 11, 2021                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methodology

% Inputs: Thrust Required, Exit area, adiabatic flame temperature,
%         Known: gamma, R, 
% Solving: area ratio, throat area, throat density, stagnation density,
%          throat temperature, exit velocity, and stagnation pressure 
% iterating on mass flow and time within each iteration

% for each mass flow value the exit velocity is found based on thrust
% requirement T = M_dot * V_e . With the exit velocity known the area ratio
% to choke can found, the throat area can be found with this. the density
% at the throat can be solved with m_dot = rho_throat*A_throat*V_throat .
% with density, combustion chamber pressure can be solved for with 
% P_0 = rho_0*R*T_0


% Inputs
Thrust = 780;  % Thrust required <N>
Exit_Area = 10;  % (preliminary) nozzle exit area <m^2>

