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

% Station 1 properties (free-stream)-----
mach(1,n) = flight_mach(n);
staticPres(1,n) = pressure_atm(n);  % <Pa>
staticDens(1,n) = Rho_atm(n);  % <kg/m^3>   
staticTemp(1,n) = Temp_atm(n);  % <K> 
[~, tempRatio(1,n), presRatio(1,n), densRatio(1,n), ~] = flowisentropic(gamma, mach(1,n), 'mach');  % ratios are static over stagnation
stagTemp(1,n) = staticTemp(1,n)/tempRatio(1,n);  % <K>                                                 
stagPres(1,n) = staticPres(1,n)/presRatio(1,n);  % <Pa>                                                 
stagDens(1,n) = staticDens(1,n)/densRatio(1,n);  % <kg/m^3>                                                   
velocity(1,n) = mach(1,n)*sqrt(gamma*R*staticTemp(1,n));  % <m/s>                                          
         
% Station 2 properties (oblique shock)----- 
[mach(2,n), shockAngle(n)] = obliqueShock(mach(1,n), def, gamma);  % uniform mach number of external ramp region                      % Oblique shock properties at flight mach
Mach1Norm(n) = mach(1,n)*sind(shockAngle(n));  % flow component crossing normal to oblique shock                                      % Normal component of mach after oblique shock
[~, tempRatio(2,n), presRatio(2,n), densRatio(2,n), ~, stagPresRatio(2,n)] = ...
    flownormalshock(gamma, Mach1Norm(n), 'mach');  % ratios are downstream over upstream
stagPres(2,n) = stagPres(1,n) * stagPresRatio(2,n);  % <Pa>
stagDens(2,n) = stagDens(1,n) * stagPresRatio(2,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
stagTemp(2,n) = stagTemp(1,n);  % does not change over normal shock <K>
staticPres(2,n) = staticPres(1,n) * presRatio(2,n);  % <Pa>
staticDens(2,n) = staticDens(1,n) * densRatio(2,n);  % <kg/m3>
staticTemp(2,n) = staticTemp(1,n) * tempRatio(2,n);  % <K>
velocity(2,n) = mach(2,n)*sqrt(gamma*R*staticTemp(2,n));  % <m/s>
[~, ~, ~, ~, A_Astar(2,n)] = flowisentropic(gamma, mach(2,n), 'mach');
Astar(2,n) = Area_intake/A_Astar(2,n);  % Area to choke before normal shock
        
% Station 3 properties (normal shock)-----
[~, tempRatio(3,n), presRatio(3,n), densRatio(3,n), mach(3,n), stagPresRatio(3,n)] = ...
    flownormalshock(gamma, mach(2,n), 'mach');  % ratios are downstream over upstream
stagPres(3,n) = stagPres(2,n) * stagPresRatio(3,n);  % <Pa>
stagDens(3,n) = stagDens(2,n) * stagPresRatio(3,n);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
stagTemp(3,n) = stagTemp(2,n);  % does not change over normal shock <K>
staticPres(3,n) = staticPres(2,n) * presRatio(3,n);  % <Pa>
staticDens(3,n) = staticDens(2,n) * densRatio(3,n);  % <kg/m3>
staticTemp(3,n) = staticTemp(2,n) * tempRatio(3,n);  % <K>
velocity(3,n) = mach(3,n)*sqrt(gamma*R*staticTemp(3,n));  % <m/s>
MdotAir(n) = staticDens(3,n)*Area_intake*velocity(3,n);  % air mass flow at intake opening <kg/s>
Astar(3,n) = Astar(2,n)/stagPresRatio(3,n);  % Area to choke after nomal shock
A_Astar(3,n) = Area_combustor/Astar(3,n);  % Area ratio of the combustor inlet  

% Station 4 properties (subsonic iscentropic expansion)-----  
[mach(4,n), tempRatio(4,n), presRatio(4,n), densRatio(4,n), ~] = ...
    flowisentropic(gamma, A_Astar(3,n), 'sub');  % ratios are static over stagnation
stagPres(4,n) = stagPres(3,n);  % does not change; iscentropic <Pa>
stagDens(4,n) = stagDens(3,n);  % does not change; iscentropic <kg/m3>
stagTemp(4,n) = stagTemp(3,n);  % does not change; iscentropi <K>
staticPres(4,n) = stagPres(4,n) * presRatio(4,n);  % <Pa>
staticDens(4,n) = stagDens(4,n) * densRatio(4,n);  % <kg/m3>
staticTemp(4,n) = stagTemp(4,n) * tempRatio(4,n);  % <K>
velocity(4,n) = mach(4,n)*sqrt(gamma*R*staticTemp(4,n));  % <m/s>

% Pass variables to ballistic simulator
InltPres_stag(n) = stagPres(4,n);                                                                 % Inlet stagnation pressure
InltPres(n) = staticPres(4,n);                                                                       % Inlet static pressure
InltRho(n) = staticDens(4,n);                                                                      % Inlet static density
InltTemp(n) = staticTemp(4,n);                                                                       % Inlet static temp
InltVel(n) = velocity(4,n);                                                                        % Inlet velocity
InltMach(n) = mach(4,n);                                                                    % Inlet mach number