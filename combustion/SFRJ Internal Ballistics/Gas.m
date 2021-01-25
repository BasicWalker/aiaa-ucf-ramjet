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
% ---------------------------------------------------------------------- %

% Set Inlet velocity and pressure
InltVel(n) = InltVel(1);
InltPres(n) = InltPres(1);

f(n) = MFuelGen(n) / SFRJDt;

% Calculate chamber pressure
PC(n) = (3/(121*(NzlAT^2))) * ((20*sqrt(5))*sqrt(121*(NzlAT^2)*(InltArea^2)*InltRho(1)*InltPres(n)...
      - 33*NzlAT*f(n)*(InltArea^2)*InltRho(1) + 4500*(InltArea^4)*InltRho(1)^2) ...
      + 11*NzlAT*f(n) - 3000*(InltArea^2)*InltRho(1));

% Throw error of chamber pressure is greater than inlet pressure  
if PC(n) > InltPres
    fprintf('Error calculating combustion chamber pressure \n')
end

% Calculate air velocity after sudden expansion inside the chamber 
AirVel(n) = InltVel(n) * (InltArea/PortArea(n));% Velocity of air

% Calculate air mass properties
MairGen(n) = AirVel(n) * InltRho(1) * SFRJDt;   % Mass of air generated
MOxdzrGen(n) = MairGen(n) * 0.21;               % Mass of oxygen generated                                    

