% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: CombustionChamber.m 
%
% File Description: 
% Simulates Combustion Chamber processes within the Ramjet to define flow
% properties
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Karam Paul      04/26/21  000  Initial Creation 
% ---------------------------------------------------------------------- %    

% iscentropic expansion for step height
[~, ~, ~, ~, combustion.A_Astar(1,n)] = flowisentropic(constants.gamma, intake.mach(end,n), 'mach');
combustion.Astar(1,n) = combustion.InletArea / combustion.A_Astar(1,n);
combustion.A_Astar(2,n) = fuel.PortArea(n) / combustion.Astar(1,n);

% flow Properties at combustor inlet
[combustion.mach(1,n), combustion.TempRatio(1,n), combustion.PresRatio(1,n), combustion.DensRatio(1,n), ~] = ...
    flowisentropic(constants.gamma, combustion.A_Astar(2,n), 'sub');  % ratios are static over stagnation
combustion.stagPres(1,n) = intake.stagPres(end,n);  % does not change; iscentropic <Pa>
combustion.stagDens(1,n) = intake.stagDens(end,n);  % does not change; iscentropic <kg/m3>
combustion.stagTemp(1,n) = intake.stagTemp(end,n);  % does not change; iscentropi <K>
combustion.staticPres(1,n) = combustion.stagPres(1,n) * combustion.PresRatio(1,n);  % <Pa>
combustion.staticDens(1,n) = combustion.stagDens(1,n) * combustion.DensRatio(1,n);  % <kg/m3>
combustion.staticTemp(1,n) = combustion.stagTemp(1,n) * combustion.TempRatio(1,n);  % <K>
combustion.velocity(1,n) = combustion.mach(1,n)*sqrt(constants.gamma*constants.R*combustion.staticTemp(1,n));  % <m/s>
combustion.massFlow(1,n) = combustion.staticDens(1,n)*combustion.velocity(1,n)*fuel.PortArea(n);

if n == 1 % intial iteration without combustion
    combustion.stagPres(2,n) = combustion.stagPres(1,n);  % does not change; iscentropic <Pa>
    combustion.stagDens(2,n) = combustion.stagDens(1,n);  % does not change; iscentropic <kg/m3>
    combustion.stagTemp(2,n) = combustion.stagTemp(1,n);  % does not change; iscentropi <K>
    combustion.staticPres(2,n) = combustion.staticPres(1,n);  % <Pa>
    combustion.staticDens(2,n) = combustion.staticDens(1,n);  % <kg/m3>
    combustion.staticTemp(2,n) = combustion.staticTemp(1,n);  % <K>
    combustion.velocity(2,n) = combustion.velocity(1,n);  % <m/s>
    combustion.massFlow(2,n) = combustion.massFlow(1,n);
    combustion.stagPresLoss(n) = 0;
else
    % use rayleigh flow to solve for properties at combustion chamber exit
    [combustion.mach(2,n), combustion.staticTemp(2,n), combustion.staticPres(2,n), combustion.stagTemp(2,n), combustion.stagPres(2,n), combustion.stagPresLoss(n)] = ...
        RayleighFlow(constants.gamma, combustion.mach(1,n), combustion.stagTemp(1,n), vehicle.AFT(n), combustion.stagPres(1,n));
    
    combustion.staticDens(2,n) = combustion.staticPres(2,n)/(combustion.staticTemp(2,n)*constants.R);
    combustion.stagDens(2,n) = combustion.stagPres(2,n)/(combustion.stagTemp(2,n)*constants.R);
    combustion.velocity(2,n) = combustion.mach(2,n)*sqrt(constants.gamma*constants.R*combustion.staticTemp(2,n));
    combustion.massFlow(2,n) = combustion.staticDens(2,n)*combustion.velocity(2,n)*fuel.PortArea(n) + fuel.MassFlow(n);  % add in the fuel mass flow to rayleigh heated air mass flow
end
% find pressure needed to choke nozzle
% intake.chokeStagPres(n+1) = pressureToChoke(intake.massFlow(2,n), nozzle.Area_throat, combustion.stagTemp(2,n)) + combustion.stagPresLoss(n);  % choke pressure AFT
intake.chokeStagPres(n+1) = pressureToChoke(combustion.massFlow(2,n), nozzle.Area_throat, combustion.stagTemp(2,n)) + combustion.stagPresLoss(n);  % choke pressure AFT


