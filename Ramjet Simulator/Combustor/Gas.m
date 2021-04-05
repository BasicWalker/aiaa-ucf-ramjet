% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Gas.m 
% 
% File Description: 
% Gas model, calculates combustion chamber pressure and air mass
% properties
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  001  Chamber Pressure Calculation Update
% Ethan Sherlock  03/12/21  002  Added Chemistry Model 
% ---------------------------------------------------------------------- %

% Calculate air mass properties
MairGen(n) = MdotAir(n)*SFRJDt;                                     % Mass of air generated (kg)
MOxdzrGen(n) = MairGen(n) * OxPercent;                              % Mass of oxygen generated by weight % (kg)        
MdotTotal(n) = MdotAir(n) + MdotFuel(n);                            % Total mass flow (kg/s)

f_yield(n) = (MFuelGen(n)/MairGen(n));                              % fuel to air ratio 

% --------------------- Chemistry Model --------------------- %
<<<<<<< HEAD:combustion/SFRJ Internal Ballistics/Gas.m
[phi, T_AFT, gamma_T, R_t] = chem.phiSolver(f_yield(n),T_2(n));     % Call Chemistry Model, need to add T_air before combustion chamber 475
T_stag(n) = T_AFT;                                                  % Stagnation temp in nozzle
gamma_t(n) = gamma_T;                                               % Gamma at the nozzle throat (old:gamma_nzl = 1.3845)
=======
[phi, T_stag(n)] = chem.phiSolver(f_yield(n),475);                      % Call Chemistry Model, need to add T_air before combustion chamber 475
gamma_nzl = 1.3845;                                                 % Temporary
gamma_nzlT(n) = gamma_nzl;                                          % Gamma at the nozzle throat
>>>>>>> f1bb9223acd5578b9b7bb61a399462fc939e2d43:Ramjet Simulator/Combustor/Gas.m
phi_eqv(n) = phi;                                                   % Grab phi value
% ----------------------------------------------------------- %

% Required chamber pressure to choke flow 
[mach,T,P,rho,area] = flowisentropic(gamma_t(n),NzlARatio,'sup');   % Isentropic Relations
PCreq(n) = (1/P)*pressure_atm(n);                                   % Calculates Stag Pres based on Pressure ratio, assumes Pstag = PC
Temp_exit(n) = T*T_AFT;                                             % Temperature at exit plane
Mach_exit(n) = mach;                                                % Mach at exit plane

% Chamber Pressure based on adiabatic flame temp
[mach, T, P, rho, area] = flowisentropic(gamma_t(n), 1 ,'mach');    % Isentropic flow conditions at nozzle throat - BIGGG Assumption
T_static(n)	= T*T_stag(n);                                          % Find static temp (K)
a_Nzl_T(n) = sqrt(gamma_t(n)*R*T_static(n));                        % Speed of sound at nozzle throat (m/s)
V_flowRate(n) = mach * a_Nzl_T(n);                                  % Flow velocity at nozzle throat (m/s)
Rho_static(n) = MdotTotal(n)/(V_flowRate(n) * NzlAT);               % Static density at nozzle throat (kg/m^3)
Rho_stag(n) = (1/rho)*Rho_static(n);                                % Stagnation density at nozzle throat (kg/m^3)
PC_TAFT(n) = Rho_stag(n)*R*T_stag(n)/Pa2kPa;                        % Stagnation pressure = chamber pressure (kPa)

% Isobaric Chamber Pressure 
PC_Isobaric(n) = InltPres_stag(n);                                  % Constant pressure combustion (kPa)

% Fuel Mass Prop & O/F Ratio calculation
OFRatio(n) = MOxdzrGen(n)/MFuelGen(n);                              % O/F Ratio
MassGen(n) = MairGen(n) + MFuelGen(n);                              % Total mass generated (kg)
MassFlow(n) = MassGen(n)/SFRJDt;                                    % Total mass flow (kg/s)