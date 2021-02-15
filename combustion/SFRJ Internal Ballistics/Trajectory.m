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

mass(1) = dry_mass + FuelMass(1);            
weight(1) = gravity*mass(1);
acceleration(1) = (Thrustdlvd(1) - drag(1) - weight(1))/ mass(1);

if (n > 1)
    velocity(n) = velocity(n-1) + acceleration(n-1)*(SFRJDt);            
    altitude(n) = altitude(n-1) + velocity(n-1)*SFRJDt + 0.5*acceleration(n-1)*SFRJDt^2;
    density(n) = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(n))/1e3);
    pressure_atm(n) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(n))/1e3);
    pressure_atm(n) = pressure_atm(n)*(1/Pa2kPa);
    temperature(n) = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(n))/1e3);
    flight_mach(n) = velocity(n)/sqrt(gamma_atm*R*temperature(n));
    drag(n) = c_d*0.5*density(n)*velocity(n)^2*S;
    mass(n) = dry_mass + FuelMass(n);
    weight(n) = gravity*mass(n);
    acceleration(n) = (Thrustdlvd(n) - drag(n) - weight(n))/ mass(n);
        
end

