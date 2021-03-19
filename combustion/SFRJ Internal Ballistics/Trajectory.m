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

if (n > 1)
    velocity(n) = velocity(n-1) + acceleration(n-1)*(SFRJDt);                               % Velocity (m/s)
    altitude(n) = altitude(n-1) + velocity(n-1)*SFRJDt + 0.5*acceleration(n-1)*SFRJDt^2;    % Altitude (m)
    Rho_atm(n) = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(n))/1e3);                     % Interpolated Lookup Table, atm density 
    pressure_atm(n) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(n))/1e3);                % Interpolated Lookup Table, atm pressure
    Temp_atm(n) = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(n))/1e3);                       % Interpolated Lookup Table, atm temp
    pressure_atm(n) = pressure_atm(n)/Pa2kPa;                                               % Convert to Kpa
    flight_mach(n) = velocity(n)/sqrt(gamma_atm*R*Temp_atm(n));                             % Vehicle Mach
    drag(n) = c_d*0.5*Rho_atm(n)*velocity(n)^2*S;                                           % Vehicle drag force
    mass(n) = dry_mass + FuelMass(n);                                                       % Vehicle mass
    weight(n) = gravity*mass(n);                                                            % Vehicle weight
    acceleration(n) = (Thrustdlvd(n-1) - drag(n-1) - weight(n-1))/ mass(n-1);               % Vehicle acceleration
end

