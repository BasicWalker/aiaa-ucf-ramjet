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
% Samer Armaly       --------  ---  Sim Revamp
% ---------------------------------------------------------------------- %

% Station intake.DeflAngleinitions used for indices
% station 1: free stream (ambient)
% station 2: downstream of oblique shock, (external ramp)
% station 3: upstream of normal shock, supersonic iscentropic expansion (internal Ramp)
% station 4: downstream of normal shock, (internal Ramp)
% station 5: combustor inlet, subsonic iscentropic expansion ,(combustion chamber opening)

% Station 1 properties (free-stream)-----
intake.mach(1,n) = vehicle.Mach(n);
intake.staticPres(1,n) = trajectory.pressure_a(n);  % <Pa>
intake.staticDens(1,n) = trajectory.Rho_a(n);  % <kg/m^3>   
intake.staticTemp(1,n) = trajectory.Temp_a(n);  % <K> 
[~, intake.tempRatio(1,n), intake.presRatio(1,n), intake.densRatio(1,n), ~] = flowisentropic(constants.gamma, intake.mach(1,n), 'mach');  % ratios are static over stagnation
intake.stagTemp(1,n) = intake.staticTemp(1,n)/intake.tempRatio(1,n);  % <K>                                                 
intake.stagPres(1,n) = intake.staticPres(1,n)/intake.presRatio(1,n);  % <Pa>                                                 
intake.stagDens(1,n) = intake.staticDens(1,n)/intake.densRatio(1,n);  % <kg/m^3>                                                   
intake.velocity(1,n) = intake.mach(1,n)*sqrt(constants.gamma*constants.R*intake.staticTemp(1,n));  % <m/s>                                          
         
% Station 2 properties (oblique shock)----- 
[intake.mach(2,n), intake.shockAngle(n)] = obliqueShock(intake.mach(1,n), intake.DeflAngle, constants.gamma);  % uniform mach number of external ramp region
intake.mach1Norm(n) = intake.mach(1,n)*sind(intake.shockAngle(n));  % flow component crossing normal to oblique shock                                  
[~, intake.tempRatio(2,n), intake.presRatio(2,n), intake.densRatio(2,n), ~, intake.stagPresRatio(2,n)] = ...
    flownormalshock(constants.gamma, intake.mach1Norm(n), 'mach');  % ratios are downstream over upstream
intake.stagPres(2,n) = intake.stagPres(1,n) * intake.stagPresRatio(2,n);  % <Pa>
intake.stagDens(2,n) = intake.stagDens(1,n) * intake.stagPresRatio(2,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
intake.stagTemp(2,n) = intake.stagTemp(1,n);  % does not change over normal shock <K>
intake.staticPres(2,n) = intake.staticPres(1,n) * intake.presRatio(2,n);  % <Pa>
intake.staticDens(2,n) = intake.staticDens(1,n) * intake.densRatio(2,n);  % <kg/m3>
intake.staticTemp(2,n) = intake.staticTemp(1,n) * intake.tempRatio(2,n);  % <K>
intake.velocity(2,n) = intake.mach(2,n)*sqrt(constants.gamma*constants.R*intake.staticTemp(2,n));  % <m/s>
[~, ~, ~, ~, intake.A_Astar(2,n)] = flowisentropic(constants.gamma, intake.mach(2,n), 'mach');
intake.Astar(2,n) = intake.Area_enter/intake.A_Astar(2,n);  % reference area before normal shock

intake.massFlow(2,n) = intake.staticDens(2,n)*intake.Area_enter*intake.velocity(2,n);  % air mass flow at intake opening <kg/s>
if n == 1  % intial iteration without combustion
    intake.chokeStagPres(n) = pressureToChoke(intake.massFlow(2,n), nozzle.Area_throat, intake.stagTemp(2,n));  % no stag temp change from combustion
end
% chokeStagPres is calculated in combustion when n > 1


% check for choked conditions after oblique shock      
if intake.stagPres(2,n) < intake.chokeStagPres(n)
    warning('intake:obliqueStagLoss',...
        'Oblique shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry decreasing altitude, intake.DeflAnglelection angle or increasing speed.');
end
% check for choked conditions after normal shock at intake enter (least stagnation pressure loss)
[~, ~, ~, ~, ~, intake.r_n_chokeLimit(n)] = flownormalshock(constants.gamma, intake.mach(2,n), 'mach');  % ratios are downstream over upstream
intake.MaxAvailableStag(n) = intake.stagPres(2,n)*intake.r_n_chokeLimit(n);
if (intake.MaxAvailableStag(n)) < intake.chokeStagPres(n)
    warning('intake:normalStagLoss',...
        'Normal shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry reducing the strength of the normal shock.');
end

% Station 3 properties (before normal shock)-----
intake.r_n = intake.chokeStagPres(n)/intake.stagPres(2,n);  % total stagnation pressure loss ratio from normal shock
% normal shock properties
[intake.mach(3,n), intake.tempRatio(4,n), intake.presRatio(4,n), intake.densRatio(4,n), intake.mach(4,n), intake.stagPresRatio(4,n), ~]...
    = flownormalshock(constants.gamma, intake.r_n, 'totalp');
% check if normal shock is ejected
if intake.mach(3,n) < intake.mach(2,n)
    warning('intake:shockEjected',...
        'Normal shock ejected; losing mass flow. \ndiffuser ratio (r_n) too high.');
end





% supersonic expansion to normal shock
intake.Astar(3,n) = intake.Astar(2,n);  % reference area same; iscentropic expansion
[~, intake.tempRatio(3,n), intake.presRatio(3,n), intake.densRatio(3,n), intake.A_Astar(3,n)]...
    = flowisentropic(constants.gamma, intake.mach(3,n), 'mach'); % ratios are static over stagnation
% solve normal shock location
intake.AShock(n) = intake.A_Astar(3,n)*intake.Astar(3,n);   % area of where the normal shock occurs
% supersonic iscentropic expansion
intake.stagPres(3,n) = intake.stagPres(2,n);  % does not change; iscentropic <Pa>
intake.stagDens(3,n) = intake.stagDens(2,n);  % does not change; iscentropic <kg/m3>
intake.stagTemp(3,n) = intake.stagTemp(2,n);  % does not change; iscentropic <K>
intake.staticPres(3,n) = intake.stagPres(3,n) * intake.presRatio(3,n);  % <Pa>
intake.staticDens(3,n) = intake.stagDens(3,n) * intake.densRatio(3,n);  % <kg/m3>
intake.staticTemp(3,n) = intake.stagTemp(3,n) * intake.tempRatio(3,n);  % <K>
intake.velocity(3,n) = intake.mach(3,n)*sqrt(constants.gamma*constants.R*intake.staticTemp(3,n));  % <m/s>
intake.massFlow(3,n) = intake.staticDens(3,n)*intake.AShock(n)*intake.velocity(3,n);  % air mass flow before normal shock <kg/s>


% Station 4 properties (after normal shock)-----
intake.stagPres(4,n) = intake.stagPres(3,n) * intake.stagPresRatio(4,n);  % <Pa>
intake.stagDens(4,n) = intake.stagDens(3,n) * intake.stagPresRatio(4,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
intake.stagTemp(4,n) = intake.stagTemp(3,n);  % does not change over normal shock <K>
intake.staticPres(4,n) = intake.staticPres(3,n) * intake.presRatio(4,n);  % <Pa>
intake.staticDens(4,n) = intake.staticDens(3,n) * intake.densRatio(4,n);  % <kg/m3>
intake.staticTemp(4,n) = intake.staticTemp(3,n) * intake.tempRatio(4,n);  % <K>
intake.velocity(4,n) = intake.mach(4,n)*sqrt(constants.gamma*constants.R*intake.staticTemp(4,n));  % <m/s>
[~, ~, ~, ~, intake.A_Astar(4,n)] = flowisentropic(constants.gamma, intake.mach(4,n), 'mach');
intake.Astar(4,n) = intake.Astar(3,n)/intake.stagPresRatio(4,n);  % reference area after normal shock 
intake.massFlow(4,n) = intake.staticDens(4,n)*intake.AShock(n)*intake.velocity(4,n);


% Station 5 properties (subsonic iscentropic expansion)-----  
intake.A_Astar(5,n) = combustion.InletArea/intake.Astar(4,n);  % Area ratio of the combustor inlet
%check for combustor inlet flow limits
if intake.A_Astar(5,n) <  3.91034275  % A_Astar value for M 0.15
   ideal_combustor_area =  3.91034275*intake.Astar(4,n);
   warning('intake:CombustorInletFast',...
       'Combustor inlet flow too fast, mach number > 0.15\n Combustor inlet too small try %g <m^2>\n',ideal_combustor_area)
end
if intake.A_Astar(5,n) < 1
    warning('intake:CombustorInletChoke',...
        'Combustor inlet too small; combustion chamber flow reaches Mach 1. Flow is choked\n');
end


[intake.mach(5,n), intake.tempRatio(5,n), intake.presRatio(5,n), intake.densRatio(5,n), ~] = ...
    flowisentropic(constants.gamma, intake.A_Astar(5,n), 'sub');  % ratios are static over stagnation
intake.stagPres(5,n) = intake.stagPres(4,n);  % does not change; iscentropic <Pa>
intake.stagDens(5,n) = intake.stagDens(4,n);  % does not change; iscentropic <kg/m3>
intake.stagTemp(5,n) = intake.stagTemp(4,n);  % does not change; iscentropi <K>
intake.staticPres(5,n) = intake.stagPres(5,n) * intake.presRatio(5,n);  % <Pa>
intake.staticDens(5,n) = intake.stagDens(5,n) * intake.densRatio(5,n);  % <kg/m3>
intake.staticTemp(5,n) = intake.stagTemp(5,n) * intake.tempRatio(5,n);  % <K>
intake.velocity(5,n) = intake.mach(5,n)*sqrt(constants.gamma*constants.R*intake.staticTemp(5,n));  % <m/s>
intake.massFlow(5,n) = intake.staticDens(5,n)*combustion.InletArea*intake.velocity(5,n);
