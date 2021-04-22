% --------- AIAA Internal Ballistic Simulator code for UCF HPR --------- %
% File Name: Thrust.m 
%
% File Description: 
% Thruster model, calculates instantaneous CStar, vehicle.thrust, specific impulse
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  03/23/21  ---  Update Thrust Calculations
% ---------------------------------------------------------------------- %

% % Thrust Calculations - Ramjet Model
% SpeedSound_exit(n) = sqrt(nozzle.gamma(n)*constants.R*Temp_exit(n));
% Velocity_exit(n) = Mach_exit(n) * SpeedSound_exit(n);
% Thrustdlvd2(n) = MdotAir(n) * ((1 + combustion.FAR(n))*Velocity_exit(n) - v_2(n));
vehicle.thrust(n) = nozzle.massFlow(2,n)*nozzle.velocity(2,n) - (intake.massFlow(2,n)*intake.velocity(2,n));

vehicle.Impulse(n) = vehicle.thrust(n)*SFRJDt;         % Impulse delivered
vehicle.Total_Impulse(n) = sum(vehicle.Impulse(:));             % Total Impulse