% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Intake.m 
% 
% File Description: 
% Intake design, oblique shock, normal shock calculation.  This code
% delivers intake flow properties to the combustion chamber inlet.
% 
% Name               Date      SCR  Description
% -----------------  --------  ---  ------------------------------
% Durlak & Aubertin  01/22/21  ---  Initial Creation 
% ---------------------------------------------------------------------- %

% Station definitions used for indices
% station 1: free stream (ambient)
% station 2: downstream of oblique shock, upstream of normal shock (external ramp)
% station 3: downstream of normal shock, upstream of subsonic iscentropic expansion (throat)
% station 4: downstream of subsonic iscentropic expansion combustor inlet,(combustion chamber opening)

intake.Area_enter =  Area_intake;

% Station 1 properties (free-stream)-----
intake.mach(1,n) = flight_mach(n);
intake.staticPres(1,n) = pressure_atm(n)*1e3;  % <Pa>
intake.staticDens(1,n) = Rho_atm(n);  % <kg/m^3>   
intake.staticTemp(1,n) = Temp_atm(n);  % <K> 
[~, intake.tempRatio(1,n), intake.presRatio(1,n), intake.densRatio(1,n), ~] = flowisentropic(gamma, intake.mach(1,n), 'mach');  % ratios are static over stagnation
intake.stagTemp(1,n) = intake.staticTemp(1,n)/intake.tempRatio(1,n);  % <K>                                                 
intake.stagPres(1,n) = intake.staticPres(1,n)/intake.presRatio(1,n);  % <Pa>                                                 
intake.stagDens(1,n) = intake.staticDens(1,n)/intake.densRatio(1,n);  % <kg/m^3>                                                   
intake.velocity(1,n) = intake.mach(1,n)*sqrt(gamma*R*intake.staticTemp(1,n));  % <m/s>                                          
         
% Station 2 properties (oblique shock)----- 
[intake.mach(2,n), intake.shockAngle(n)] = obliqueShock(intake.mach(1,n), def, gamma);  % uniform mach number of external ramp region
intake.mach1Norm(n) = intake.mach(1,n)*sind(intake.shockAngle(n));  % flow component crossing normal to oblique shock                                  
[~, intake.tempRatio(2,n), intake.presRatio(2,n), intake.densRatio(2,n), ~, intake.stagPresRatio(2,n)] = ...
    flownormalshock(gamma, intake.mach1Norm(n), 'mach');  % ratios are downstream over upstream
intake.stagPres(2,n) = intake.stagPres(1,n) * intake.stagPresRatio(2,n);  % <Pa>
intake.stagDens(2,n) = intake.stagDens(1,n) * intake.stagPresRatio(2,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
intake.stagTemp(2,n) = intake.stagTemp(1,n);  % does not change over normal shock <K>
intake.staticPres(2,n) = intake.staticPres(1,n) * intake.presRatio(2,n);  % <Pa>
intake.staticDens(2,n) = intake.staticDens(1,n) * intake.densRatio(2,n);  % <kg/m3>
intake.staticTemp(2,n) = intake.staticTemp(1,n) * intake.tempRatio(2,n);  % <K>
intake.velocity(2,n) = intake.mach(2,n)*sqrt(gamma*R*intake.staticTemp(2,n));  % <m/s>
[~, ~, ~, ~, intake.A_Astar(2,n)] = flowisentropic(gamma, intake.mach(2,n), 'mach');
intake.Astar(2,n) = intake.Area_enter/intake.A_Astar(2,n);  % reference area before normal shock

intake.massFlow(2,n) = intake.staticDens(2,n)*intake.Area_enter*intake.velocity(2,n);  % air mass flow at intake opening <kg/s>
if n < 2
    intake.chokeStagPres(n) = pressureToChoke(intake.massFlow(2,n), nozzle.Area_throat, intake.stagTemp(2,n));  % no stag temp change from combustion
else
    intake.chokeStagPres(n) = pressureToChoke(intake.massFlow(2,n), nozzle.Area_throat, combustion.stagTemp(2,n-1));  % choke pressure using previous AFT
end

% we need to accelerate until normal shock location to make the mass flow
% at 3 work




% check for choked conditions       
if intake.stagPres(2,n) < intake.chokeStagPres(n)
    error('intake:obliqueStagLoss',...
        'Error. \nOblique shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry decreasing altitude, deflection angle or increasing speed.');
end
[~, ~, ~, ~, ~, intake.r_d_chokeLimit(n)] = flownormalshock(gamma, intake.mach(2,n), 'mach');  % ratios are downstream over upstream
if (intake.stagPres(2,n)*intake.r_d_chokeLimit(n)) < intake.chokeStagPres(n)
    error('intake:normalStagLoss',...
        'Error. \nNormal shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry reducing the strength of the normal shock.');
end

% Station 3 properties (normal shock)-----
% solve normal shock placement
intake.r_d = intake.chokeStagPres(n)/intake.stagPres(2,n);
[intake.machShock(n), intake.tempRatio(3,n), intake.presRatio(3,n), intake.densRatio(3,n), intake.mach(3,n), intake.stagPresRatio(3,n), ~]...
    = flownormalshock(gamma, intake.r_d, 'totalp');
if intake.machShock(n) < intake.mach(2,n)
    error('intake:shockEjected',...
        'Error. \nNormal shock ejected; losing mass flow. \ndiffuser ratio (r_d) too high.');
end

[~, ~, ~, ~, intake.Ashock_A2star] = flowisentropic(gamma, intake.machShock(n), 'mach');
intake.AShock(n) = intake.Ashock_A2star*intake.Astar(2,n);   % area of where the normal shock occurs

intake.stagPres(3,n) = intake.stagPres(2,n) * intake.stagPresRatio(3,n);  % <Pa>
intake.stagDens(3,n) = intake.stagDens(2,n) * intake.stagPresRatio(3,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
intake.stagTemp(3,n) = intake.stagTemp(2,n);  % does not change over normal shock <K>
intake.staticPres(3,n) = intake.staticPres(2,n) * intake.presRatio(3,n);  % <Pa>
intake.staticDens(3,n) = intake.staticDens(2,n) * intake.densRatio(3,n);  % <kg/m3>
intake.staticTemp(3,n) = intake.staticTemp(2,n) * intake.tempRatio(3,n);  % <K>
intake.velocity(3,n) = intake.mach(3,n)*sqrt(gamma*R*intake.staticTemp(3,n));  % <m/s>
[~, ~, ~, ~, intake.A_Astar(3,n)] = flowisentropic(gamma, intake.mach(3,n), 'mach');
intake.Astar(3,n) = intake.Astar(2,n)/intake.stagPresRatio(3,n);  % reference area after normal shock 
% intake.AShock = intake.A_Astar(3,n)*intake.Astar(3,n);  
intake.massFlow(3,n) = intake.staticDens(3,n)*intake.AShock(n)*intake.velocity(3,n);

% Station 4 properties (subsonic iscentropic expansion)-----  
intake.A_Astar(4,n) = Area_combustor/intake.Astar(3,n);  % Area ratio of the combustor inlet
[intake.mach(4,n), intake.tempRatio(4,n), intake.presRatio(4,n), intake.densRatio(4,n), ~] = ...
    flowisentropic(gamma, intake.A_Astar(4,n), 'sub');  % ratios are static over stagnation
intake.stagPres(4,n) = intake.stagPres(3,n);  % does not change; iscentropic <Pa>
intake.stagDens(4,n) = intake.stagDens(3,n);  % does not change; iscentropic <kg/m3>
intake.stagTemp(4,n) = intake.stagTemp(3,n);  % does not change; iscentropi <K>
intake.staticPres(4,n) = intake.stagPres(4,n) * intake.presRatio(4,n);  % <Pa>
intake.staticDens(4,n) = intake.stagDens(4,n) * intake.densRatio(4,n);  % <kg/m3>
intake.staticTemp(4,n) = intake.stagTemp(4,n) * intake.tempRatio(4,n);  % <K>
intake.velocity(4,n) = intake.mach(4,n)*sqrt(gamma*R*intake.staticTemp(4,n));  % <m/s>
intake.massFlow(4,n) = intake.staticDens(4,n)*Area_combustor*intake.velocity(4,n);

% Pass variables to ballistic simulator
InltPres_stag(n) = intake.stagPres(4,n);  % Inlet stagnation pressure
InltTemp_stag(n) = intake.stagTemp(4,n);  % Inlet stagnation temperature
InltTemp_dens(n) = intake.stagDens(4,n);  % Inlet stagnation density
InltPres(n) = intake.staticPres(4,n);  % Inlet static pressure
InltRho(n) = intake.staticDens(4,n);  % Inlet static density
InltTemp(n) = intake.staticTemp(4,n);  % Inlet static temp
InltVel(n) = intake.velocity(4,n);  % Inlet velocity
Inltmach(n) = intake.mach(4,n);  % Inlet mach number