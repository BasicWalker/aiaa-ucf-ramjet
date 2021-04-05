%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name:  Nozzle
% 
% Program Description: 
%
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
% 
% File Name: Chamber_pressure_study.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Samer Armaly    02/11/21  Initial Creation 
% --------------------------------------------------------------------- %
clear; clc; close all;
% Iteration parameters
% study variable
AFT_min = 1500;  % lowest adiabatic flame temp to iterate between <kg/s>
AFT_max = 3300;  % highest adiabatic flame temp to iterate between <kg/s>
i_step = 5;  % stepsize to iterate mass flow
% % time
% burntime = 12;  % full burntime of ramjet <s>

% Pre-allocate variable arrays
AFT = [AFT_min:i_step:AFT_max];  % Mass flow array
A_exit = pi()*1.5^2;
% 
[chamber_pressure,rho_stag] = chamberPressure(AFT,500,4,A_exit)

function [P_stag, rho_stag] = chamberPressure(AFT, thrust, mach_exit, A_exit, gamma, R)
% set default values if parameters arent specified
    if ~exist('gamma','var')
        gamma = 1.4;  % specific heat ratio <unitless>
    end
    if ~exist('R','var')
        R = 287;  % gas constant air <J/(kg*K)>
    end
    
    % pre-allocate variable arrays
    P_stag = zeros(1,size(AFT,2));
    T_e = zeros(1,size(AFT,2));
    v_e = zeros(1,size(AFT,2));
    m_dot = zeros(1,size(AFT,2));
    rho_throat = zeros(1,size(AFT,2));
    rho_stag = zeros(1,size(AFT,2));
    T_throat = zeros(1,size(AFT,2));
    v_throat = zeros(1,size(AFT,2));
    P_exit = zeros(1,size(AFT,2));
    %     thrust = zeros(1,size(AFT,2)) + thrust;  % <N>
    %     m_dot = zeros(1,size(AFT,2)) + m_dot;  % mass flow <kg/s>
    %     a_exit = zeros(1,size(AFT,2)) + a_exit;  % nozzle area exit <m^2>
    %     gamma = zeros(1,size(AFT,2)) + gamma;  % specific heat ratio for air <unitless>
    %     R = zeros(1,size(AFT,2)) + R;  % gas constant for air <J/(kg*K)>
    %     v_exit = thrust./m_dot;  % exit velocity <m/s>

    for i=1:size(AFT,2)
        [~, T_rat_exit, P_rat_exit, rho_rat_exit, area_rat_exit] = flowisentropic(gamma,mach_exit);
        T_e(i) = T_rat_exit*AFT(i);  % exit static temperature <K>
        v_e(i) = mach_exit*sqrt(gamma*R*T_e(i));  % exit velocity <m/s>
        m_dot(i) = thrust/v_e(i);  % total mass flow rate <kg/s>
        A_throat = A_exit/area_rat_exit;  % throat area <m^2>
        [~, T_rat_throat, P_rat_throat, rho_rat_throat, ~] = flowisentropic(gamma,mach_exit);
        T_throat(i) = T_rat_throat*AFT(i);  % throat static temperature <K>
        v_throat(i) = 1*sqrt(gamma*R*T_throat(i)); % velocity at throat (mach 1) <m/s>
        rho_throat(i) = m_dot(i)./(A_throat*v_throat(i));  % <kg/m^3>
        rho_stag(i) = rho_throat(i)/rho_rat_throat;  % stagnation density <kg/m^3>
        P_stag(i) = (rho_stag(i)*AFT(i))*R;
        P_exit(i) = P_stag(i)*P_rat_exit;
%         thrust_actual = m_dot*v_exit + ()
    end
    
end



