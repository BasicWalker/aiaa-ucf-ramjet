% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Intake.m 
% 
% File Description: 
% Intake design, oblique shock, normal shock calculation.  This code
% delivers intake flow properties to the combustion chamber inlet.
% 
% Name               Date      SCR  Description
% -----------------  --------  ---  ------------------------------
% Durlak & Aubertin  01/22/21  ---  Initial Creation 
% ---------------------------------------------------------------------- %
    
% Station 1
[~, Ti_To1(n), Pi_Po1(n), Rho_Rho_o(n), ~] = flowisentropic(gamma, Mach_f(n));         % Isentropic flow properties at flight mach
T_o1(n) = (1/Ti_To1(n)) * Temp_a(n);                                                      % Stag Temp 1
Po1(n) = (pressure_atm(n)/Pi_Po1(n));                                                       % Stag Pres 1
Rho_o(n) = Rho_a(n)/Rho_Rho_o(n);                                                         % Stag Density 1
v_1(n) = Mach_f(n)*sqrt(gamma*R*Temp_a(n));                                          % Velocity at station 1
         
% Station 2 (oblique shock) 
[mach_2(n), beta_angle(n)] = obliqueShock(Mach_f(n), def, gamma);                      % Oblique shock properties at flight mach
mach_2_normal(n) = Mach_f(n)*sind(beta_angle(n));                                      % Normal component of mach after oblique shock
[~, T2_To2(n), P2_Po2(n), Rho_Rho_o2(n), A_Astar_2(n)] = flowisentropic(gamma, mach_2(n));  % Isentropic flow properties at station 2 mach
[~, T2_T1(n), ~, ~, ~, Po2_Po1(n), ~] = flownormalshock(gamma, mach_2_normal(n));           % Normal shock properties at mach2 normal component
T_2(n) = T2_T1(n)*Temp_a(n);                                                              % Static Temp 2
Po_2(n) = Po2_Po1(n)*(1/Pi_Po1(n))*pressure_atm(n);                                         % Stag Pres 2
P_2(n) = P2_Po2(n)*Po_2(n);                                                                 % Static Pres 2
Rho_o_2(n) = Rho_o(n)*Po2_Po1(n);                                                           % Stag Density 2
Rho_2(n) = Rho_Rho_o2(n)*Rho_o_2(n);                                                        % Static Density 2
T_o2(n) = T_o1(n);                                                                          % Stag Temp 2
v_2(n) = mach_2(n)*sqrt(gamma*R*T_2(n));                                                    % Velocity at station 2
        
% Station 3 (normal shock)
mach_3(n) = normalShock(mach_2(n), gamma);                                                  % Normal shock calc to solve for mach3
[~, T3_T2(n), ~, ~, ~, Po3_Po2(n), ~] = flownormalshock(gamma, mach_2(n));                  % Normal shock properties at mach2
[~, T3_To3(n), P3_Po3(n), Rho_Rho_o3(n), ~] = flowisentropic(gamma, mach_3(n));             % Isentropic flow properties at station 3 mach
T_3(n) = T3_T2(n)*T2_T1(n)*Ti_To1(n)*T_o1(n);                                               % Static Temp 3
To_3(n) = T_3(n) * (1/T3_To3(n));                                                           % Stag Temp 3
Po_3(n) = Po3_Po2(n)*Po2_Po1(n)*(1/Pi_Po1(n))*pressure_atm(n);                              % Stag Pres 3
P_3(n) = P3_Po3(n)*Po_3(n);                                                                 % Static Pres 3
Rho_o_3(n) = Rho_o_2(n)*Po3_Po2(n);                                                         % Stag Density 3
Rho_3(n) = Rho_Rho_o3(n)*Rho_o_3(n);                                                        % Static Density 3
v_3(n) = mach_3(n)*sqrt(gamma*R*T_3(n));                                                    % Velocity at station 3
MdotAir(n) = Rho_3(n)*Area_3*mach_3(n)*sqrt(gamma*R*T_3(n));                                % Mass flow rate
        
% Station 4 (combustor)
Astar_2(n) = Area_3/A_Astar_2(n);                                                           % Area to choke before normal shock
Astar_3(n) = Astar_2(n)/Po3_Po2(n);                                                         % Area to choke after nomal shock
Acombustor_A3star(n) = Area_combustor/Astar_3(n);                                           % Area ratio of the combustor inlet
[mach_4(n), T4_To4(n), P4_Po4(n), Rho_Rho_o4(n), ~] = flowisentropic(gamma, Acombustor_A3star(n), 'sub');
To_4(n) =  To_3(n);                                                                         % Stag Temp 4
Po_4(n) = Po_3(n);                                                                          % Stag Pres 4
Rho_o_4(n) = Rho_o_3(n);                                                                    % Stag Density 4
T_4(n) = T4_To4(n)*To_4(n);                                                                 % Static Temp 4
P_4(n) = P4_Po4(n)*Po_4(n);                                                                 % Static Pres 4
Rho_4(n) = Rho_Rho_o4(n)*Rho_o_4(n);                                                        % Static Density 4
v_4(n) = mach_4(n)*sqrt(gamma*R*T_4(n));                                                    % Velocity at station 4

% Pass variables to ballistic simulator
InltPres_stag(n) = Po_4(n);                                                                 % Inlet stagnation pressure
InltPres(n) = P_4(n);                                                                       % Inlet static pressure
InltRho(n) = Rho_4(n);                                                                      % Inlet static density
InltTemp(n) = T_4(n);                                                                       % Inlet static temp
InltVel(n) = v_4(n);                                                                        % Inlet velocity
InltMach(n) = mach_4(n);                                                                    % Inlet mach number