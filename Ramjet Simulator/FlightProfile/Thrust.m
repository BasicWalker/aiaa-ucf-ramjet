% --------- AIAA Internal Ballistic Simulator code for UCF HPR --------- %
% File Name: Thrust.m 
%
% File Description: 
% Thruster model, calculates instantaneous CStar, thrust, specific impulse
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  03/23/21  ---  Update Thrust Calculations
% ---------------------------------------------------------------------- %

% Calculate CStar based on O/F Ratio
if (OFRatio(n) < 5 && OFRatio(n) > 3)
     CStar(n) = - 16.25 * (OFRatio(n) - 5)^2 + 1585 + abs((7500 - InltPres(n))/7500)*55;
elseif (OFRatio(n) <= 3)
     CStar(n) = 125 * (OFRatio(n) - 1) + 1250;
else
     CStar(n) = (-115/7) * (OFRatio(n) - 5) + 1585;
end        

if(CStar(n) < 0)
    CStar(n) = 0.0;
end 

% Thrust Calculations - Legacy Model
Thrustdlvd(n) = MassGen(n)*CStar(n)/SFRJDt; % Thrust delivered, assuming thrust coefficient = 1.0

% Thrust Calculations - Ramjet Model
SpeedSound_exit(n) = sqrt(gamma_t(n)*R*Temp_exit(n));
Velocity_exit(n) = Mach_exit(n) * SpeedSound_exit(n);
Thrustdlvd2(n) = MdotAir(n) * ((1 + f_yield(n))*Velocity_exit(n) - v_2(n));

Impulse(n) = Thrustdlvd2(n)*SFRJDt;         % Impulse delivered
TotallImp(n) = sum(Impulse(:));             % Total Impulse