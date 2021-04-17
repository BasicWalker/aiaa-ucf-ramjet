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

% % Calculate CStar based on O/F Ratio
% if (combustion.OFRatio(n) < 5 && combustion.OFRatio(n) > 3)
%      CStar(n) = - 16.25 * (combustion.OFRatio(n) - 5)^2 + 1585 + abs((7500 - intake.staticPres(end,n))/7500)*55;
% elseif (combustion.OFRatio(n) <= 3)
%      CStar(n) = 125 * (combustion.OFRatio(n) - 1) + 1250;
% else
%      CStar(n) = (-115/7) * (combustion.OFRatio(n) - 5) + 1585;
% end        
% 
% if(CStar(n) < 0)
%     CStar(n) = 0.0;
% end 

% % Thrust Calculations - Legacy Model
% Thrustdlvd(n) = MassGen(n)*CStar(n)/SFRJDt; % Thrust delivered, assuming thrust coefficient = 1.0

% % Thrust Calculations - Ramjet Model
% SpeedSound_exit(n) = sqrt(nozzle.gamma(n)*R*Temp_exit(n));
% Velocity_exit(n) = Mach_exit(n) * SpeedSound_exit(n);
% Thrustdlvd2(n) = MdotAir(n) * ((1 + combustion.FAR(n))*Velocity_exit(n) - v_2(n));
thrust.F_t(n) = nozzle.massFlow(2,n)*nozzle.velocity(2,n) - (intake.massFlow(2,n)*intake.velocity(2,n));

Impulse(n) = thrust.F_t(n)*SFRJDt;         % Impulse delivered
TotallImp(n) = sum(Impulse(:));             % Total Impulse