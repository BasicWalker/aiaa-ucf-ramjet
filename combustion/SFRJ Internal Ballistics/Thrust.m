% --------- AIAA Internal Ballistic Simulator code for UCF HPR --------- %
% File Name: Thrust.m 
%
% File Description: 
% Thruster model, calculates instantaneous CStar, thrust, specific impulse
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
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

% Thrust Calculations 
Thrustdlvd(n) = MassGen(n)*CStar(n)/SFRJDt; % Thrust delivered, assuming thrust coefficient = 1.0
Impulse(n) = Thrustdlvd(n)*SFRJDt;          % Impulse delivered
TotallImp(n) = sum(Impulse(:));             % Total Impulse