% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Trajectory.m 
% 
% File Description: 
% Trajectory model for the SFRJ
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  005  1DOF trajectory update
% ---------------------------------------------------------------------- %

mass(1)         = dry_mass + FuelMass(1);                                                   % Mass of Vehicle
weight(1)       = gravity*mass(1);                                                          % Weight of Vehicle
acceleration(1) = (Thrustdlvd(1) - drag(1) - weight(1))/ mass(1);                           % Initial Acceleration 

% Interpolated Lookup Tables - Atmospheric values
% pressure_atm(1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);
% pressure_atm(1) = pressure_atm(n)*(1/Pa2kPa);
% Temp_atm(1) = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);
% Rho_atm(1) = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);

if (n > 1)
    
    velocity(n) = velocity(n-1) + acceleration(n-1)*(SFRJDt);                               % Velocity
    altitude(n) = altitude(n-1) + velocity(n-1)*SFRJDt + 0.5*acceleration(n-1)*SFRJDt^2;    % Altitude
    Rho_atm(n) = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(n))/1e3);                     % Interpolated Lookup Table, atm density
    pressure_atm(n) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(n))/1e3);                % Interpolated Lookup Table, atm pressure
    Temp_atm(n) = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(n))/1e3);                       % Interpolated Lookup Table, atm temp
    pressure_atm(n) = pressure_atm(n)*(1/Pa2kPa);                                           % Convert to Kpa
    flight_mach(n) = velocity(n)/sqrt(gamma_atm*R*Temp_atm(n));                             % Vehicle Mach
    drag(n) = c_d*0.5*Rho_atm(n)*velocity(n)^2*S;                                           % Vehicle drag force
    mass(n) = dry_mass + FuelMass(n);                                                       % Vehicle mass
    weight(n) = gravity*mass(n);                                                            % Vehicle weight
    acceleration(n) = (Thrustdlvd(n-1) - drag(n) - weight(n))/ mass(n);                     % Vehicle acceleration
        
end

