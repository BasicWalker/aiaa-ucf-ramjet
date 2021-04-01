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


    velocity(n+1) = velocity(n) + acceleration(n)*(SFRJDt);                               % Velocity (m/s)
    altitude(n+1) = altitude(n) + velocity(n)*SFRJDt + 0.5*acceleration(n)*SFRJDt^2;    % Altitude (m)
    acceleration(n+1) = (Thrustdlvd2(n) - drag(n) - weight(n))/ mass(n);              % Vehicle acceleration
                                                           
    Rho_atm(n+1) = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(n+1))/1e3);                     % Interpolated Lookup Table, atm density 
    pressure_atm(n+1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(n+1))/1e3)/Pa2kPa;         % Interpolated Lookup Table, atm pressure
    Temp_atm(n+1) = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(n+1))/1e3);                       % Interpolated Lookup Table, atm temp
    flight_mach(n+1) = velocity(n+1)/sqrt(gamma_atm*R*Temp_atm(n+1));                             % Vehicle Mach
    drag(n+1) = c_d*0.5*Rho_atm(n+1)*velocity(n+1)^2*S;                                           % Vehicle drag force
    mass(n+1) = dry_mass + FuelMass(n);                                                       % Vehicle mass
    weight(n+1) = gravity*mass(n);  % Vehicle weight

