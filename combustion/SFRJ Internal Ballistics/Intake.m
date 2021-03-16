% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name:  Intake
% 
% Program Description: 
%
% 
% File Name: NozzleCalculation.m
% 
% File Description: 
% 
% Name                               Date      Description
% ----------------------------      --------  -------------------------
% Jared Durlak & Matt Aubertin      01/22/21  Initial Creation 
% --------------------------------------------------------------------- %

% solving mass flow rate      
% Station 1
[~, Ti_To1(n), Pi_Po1(n), Rho_Rho_o(n), ~] = flowisentropic(gamma, flight_mach(n));
T_o1(n) = (1/Ti_To1(n)) * Temp_atm(n);
Po1(n) = (pressure_atm(n)/Pi_Po1(n));
Rho_o(n) = Rho_atm(n)/Rho_Rho_o(n);
v_1(n) = flight_mach(n)*sqrt(gamma*R*Temp_atm(n));
         
% Station 2 (oblique shock) 
[mach_2(n), beta_angle] = obliqueShock(flight_mach(n), def, gamma);
mach_2_normal(n) = flight_mach(n)*sind(beta_angle);
[~, T2_To2(n), P2_Po2(n), Rho_Rho_o2(n), A_Astar_2(n)] = flowisentropic(gamma, mach_2(n));
[~, T2_T1(n), ~, ~, ~, Po2_Po1(n), ~] = flownormalshock(gamma, mach_2_normal(n));
T_2(n) = T2_T1(n)*Temp_atm(n);
Po_2(n) = Po2_Po1(n)*(1/Pi_Po1(n))*pressure_atm(n);
P_2(n) = P2_Po2(n)*Po_2(n);
Rho_o_2(n) = Rho_o(n)*Po2_Po1(n);
Rho_2(n) = Rho_Rho_o2(n)*Rho_o_2(n);
T_o2(n) = T_o1(n);
v_2(n) = mach_2(n)*sqrt(gamma*R*T_2(n));
        
% Station 3 (normal shock)
mach_3(n) = normalShock(mach_2(n), gamma);
[~, T3_T2(n), ~, ~, ~, Po3_Po2(n), ~] = flownormalshock(gamma, mach_2(n));
[~, T3_To3(n), P3_Po3(n), Rho_Rho_o3(n), ~] = flowisentropic(gamma, mach_3(n));     
T_3(n) = T3_T2(n)*T2_T1(n)*Ti_To1(n)*T_o1(n);
To_3(n) = T_3(n) * (1/T3_To3(n)); 
Po_3(n) = Po3_Po2(n)*Po2_Po1(n)*(1/Pi_Po1(n))*pressure_atm(n);
P_3(n) = P3_Po3(n)*Po_3(n);
Rho_o_3(n) = Rho_o_2(n)*Po3_Po2(n);
Rho_3(n) = Rho_Rho_o3(n)*Rho_o_3(n);
v_3(n) = mach_3(n)*sqrt(gamma*R*T_3(n));
m_dot(n) = Rho_3(n)*Area_3*mach_3(n)*sqrt(gamma*R*T_3(n));
        
% Station 4 (combustor)
Astar_2(n) = Area_3/A_Astar_2(n);
Astar_3(n) = Astar_2(n)/Po3_Po2(n);
Acombustor_A3star(n) = Area_combustor/Astar_3(n);
[mach_4(n), T4_To4(n), P4_Po4(n), Rho_Rho_o4(n), ~] = flowisentropic(gamma, Acombustor_A3star(n), 'sub');
To_4(n) =  To_3(n); 
Po_4(n) = Po_3(n);
Rho_o_4(n) = Rho_o_3(n);
T_4(n) = T4_To4(n)*To_4(n);
P_4(n) = P4_Po4(n)*Po_4(n);
Rho_4(n) = Rho_Rho_o4(n)*Rho_o_4(n);
v_4(n) = mach_4(n)*sqrt(gamma*R*T_4(n));

InltPres_stag(n) = Po_4(n);                     % Inlet stagnation pressure
InltRho(n) = Rho_4(n);                          % Inlet static density
InltTemp(n) = T_4(n);                           % Inlet static temp
InltVel(n) = v_4(n);                            % Inlet velocity
InltMach(n) = mach_4(n);                        % Inlet mach number