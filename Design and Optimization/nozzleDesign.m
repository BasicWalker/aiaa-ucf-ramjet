% ---------- SFRJ Nozzle Design Tool / UCF CAPSTONE PROJECT ------------ %
% File Name: nozzleDesign.m 
% 
% File Description: 
% The primary objective of this script is to design the nozzle throat area
% such that the chamber (stagnation) pressure present in the combustor is
% ALWAYS capable of providing choked flow.
%
% PARAMETERS:
% Area Throat:              Limits the choking conditions of the ramjet
% Stagnation Temperature:   AFT - Changes primarily with air mass flow
% Chamber Pressure:         Shall be determined as the best operating mode of the SFRJ
% Air Mass Flow:            Determined by intake area. Changes throughout flight 
%
% ASSUMPTIONS:
% This code shall assume that the nozzle is choked. Ideally, there should
% be a warning when the ramjet is operating in an unchoked mode, however,
% this would be best suited as a simulation warning callout.
% 
% 
% Name                      Date        SCR             Description
% ------------------      --------      ---     ------------------------------
% Karam, Jason            01/22/21      000           Initial Creation 
% ---------------------------------------------------------------------- %
clc; clear; close all

intake = Intake_PropertyCalculator();
chem = Chemistry();
thrust = ThrustCalculator();

% Notes
    % Iterate on throat area (limited by expansion ratio)
    % Set an intermediate fuel mass flow
    % Iterate on air mass flow
    % Yield iterative fuel-air ratio
    % Yield iterative phi's (air mass flow related to specific fuel mass flow)
    % Plug in variable phi
    % Yield AFT 
    % Use NASA mass flow equation to solve for stag pressure
    % Find reasonable stag pressure
    % Ensure that sufficient thrust is achieved with air mass flow

% Declare Set Variables 
gamma = 1.4; 
R = 287;
intakeStaticTemp = 475;
f_st = 0.0819;
fuelMassFlow = 0.09;  % <kg/s>
in2m = 0.0254;

% Iterate on altitude
altStep = 1000;
min_alt = 1000;
max_alt = 30000;

% Iterate on Air Mass Flow Rate
airMassFlowStep = 0.01;
min_airMassFlow = 0.5;
max_airMassFlow = 3;

% Iterate on Throat Area <in>
throatDiameterStep = 0.1;
min_throatDiameter = 0.5; 
max_throatDiameter = 2.5;

% Array Sizes 
l = ((max_alt - min_alt) / altStep) + 1;
m = ((max_airMassFlow - min_airMassFlow) / airMassFlowStep) + 1;
n = ((max_throatDiameter - min_throatDiameter) / throatDiameterStep) + 1;

% Pre-allocate Resources (Iteration 1)
altitude = zeros(l,1);
initialStagPres = zeros(l,1);
finalStagPres = zeros(l,1);
stagPresLoss = zeros(l,1);
velocity = zeros(l,1);

% Pre-allocate Resources (Iteration 2)
airMassFlow     = zeros(m,n);
totalMassFlow   = zeros(m,n);
f               = zeros(m,n);
phi             = zeros(m,n);
T_AFT           = zeros(m,n);

throatDiameter_in   = zeros(m,n);
throatArea_in       = zeros(m,n);
throatDiameter_m    = zeros(m,n);
throatArea_m        = zeros(m,n);

chamberPres         = zeros(m,n);

exitArea_in         = zeros(m,n);
exitArea_m          = zeros(m,n);
exitDiameter_in     = zeros(m,n);
exitDiameter_m      = zeros(m,n);

exitPres        = zeros(m,n);
exitTemp        = zeros(m,n);
exitVelocity    = zeros(m,n);
idealThrust     = zeros(m,n);

%% Begin Interation 1
    % Need to iterate on intake designs to understand achievable stagnation pressures
intake_mach         = 2;
intake_deflection   = 15;
for i = 1:l
    altitude(i) = min_alt + altStep*(i-1);
    
    [initialStagPres(i), finalStagPres(i), stagPresLoss(i)] = ...
        intake.StagnationLoss(intake_mach,altitude(i),intake_deflection);
    
    velocity(i) = intake.Velocity(intake_mach,altitude(i));
end

%% Begin Iteration 2
    % Iterate on throat area 
    % Set an intermediate fuel mass flow
    % Iterate on air mass flow
    % Yield iterative fuel-air ratio
    % Yield iterative phi's (air mass flow related to specific fuel mass flow)
    % Plug in variable phi
    % Yield AFT 
    % Use NASA mass flow equation to solve for stag pressure
    % Find reasonable stag pressure

% Enter Expansion Ratio
expansionRatio = 2;
    % expansionRatio = input("What expansion ratio did you want to study?\n");

counter = 1;
for i = 1:m
    for j = 1:n

        % Calcuate AFT for each air mass flow rate 
        airMassFlow(i,j) = min_airMassFlow + airMassFlowStep*(i-1);
        totalMassFlow(i,j) = airMassFlow(i,j) + fuelMassFlow;
        f(i,j) = fuelMassFlow / airMassFlow(i,j);
        [phi(i,j), T_AFT(i,j)] = chem.phiSolver(f(i,j), intakeStaticTemp);

        % Calculate throat area in m2
        throatDiameter_in(i,j) = min_throatDiameter + throatDiameterStep*(j-1);
        throatArea_in(i,j) = (pi/4) * throatDiameter_in(i,j)^2;
        throatDiameter_m(i,j) = throatDiameter_in(i,j) * in2m;
        throatArea_m(i,j) = (pi/4) * throatDiameter_m(i,j)^2;

        % Calculate stagnation pressure to choke nozzle with this mass flow
        chamberPres(i,j) = (totalMassFlow(i,j) * sqrt(T_AFT(i,j))) / ( throatArea_m(i,j) * (sqrt(gamma/R)) ...
                            * ((gamma+1)/2)^-((gamma+1)/(2*(gamma-1))));  

        % Calculate Thrust Given an Area Ratio
        exitArea_in(i,j) = throatArea_in(i,j) * expansionRatio; 
        exitDiameter_in(i,j) = sqrt((4/pi)*exitArea_in(i,j));
        exitArea_m(i,j) = throatArea_m(i,j) * expansionRatio;
        exitDiameter_m(i,j) = sqrt((4/pi)*exitArea_m(i,j));

        [exitMach, tempRatio, presRatio, ~,~] = flowisentropic(gamma,expansionRatio,'sup');
        exitPres(i,j) = chamberPres(i,j) * presRatio;
        exitTemp(i,j) = T_AFT(i,j) * tempRatio;
        exitVelocity(i,j) = exitMach * sqrt(gamma*R*exitTemp(i,j));
    end

    fprintf("Computing... %.3f%% Complete\n",(counter/m)*100);
    counter = counter + 1;
end

%% Plot Data
plotFlag = 1;

if plotFlag == 1
    
    figure(1)
    plot(altitude, finalStagPres);
    title('Stagnation Pressure after Normal Shock vs. Altitude');
    xlabel('Altitude < m >')
    ylabel('Stagnation Pressure after Shock < Pa >')
    
    a = airMassFlow(:,1);
    b = throatDiameter_in(1,:);
    [X,Y] = ndgrid(a,b);
    figure(2)
    mesh(X, Y, chamberPres)
    xlabel('Air Mass Flow Rate <kg/s>')
    ylabel('Throat Diameter <in>')
    zlabel('Minimum Chamber Pressure to Choke < Pa >');
    
    c = phi(:,1);
    d = throatDiameter_in(1,:);
    [X,Y] = ndgrid(c,d);
    figure(3)
    mesh(X, Y, chamberPres)
    xlabel('Equivalence Ratio')
    ylabel('Throat Diameter <in>')
    zlabel('Minimum Chamber Pressure to Choke < Pa >');
    
    
end


%% Call Thrust Function
% Input Parameters for Thrust Function
%     airMassFlow     = 0.5;          % <kg/s>
%     f               = 0.16;         % <unitless>
%     exitVelocity    = 1740.5;       % <m/s>
%     flightVelocity  = 1006;         % <m/s>
%     Pe              = 5.997e5;      % <Pa>
%     Pa              = 80325;        % <Pa>
%     Ae              = 2.5335e-4;    % <m2>

% If you want Ideal Thrust
    % Ideal_Thrust = thrust.ideal_Thrust(airMassFlow, f, exitVelocity, flightVelocity)

% If you want Non-Ideal Thrust
    % Non_Thrust = thrust.nonideal_Thrust(airMassFlow, f, exitVelocity, flightVelocity, Pe, Pa, Ae)

% Target a phi or AFT and then find a mass flow for that
% Iterate on altitude to find a reasonable chamber pressure

% Choose a max altitude at a mach number --> final stag pressure
% Set final stag = chamber pressure
% unknowns: mass flow rate and AFT and A_intake and A_throat
% equations: 
    %   nozzle mass flow 
    %   equivalence ratio equation
    %   intake mass flow
    %   mass flow rate at throat
    
% one altitude, one mach number (sets stag pressures)
% plot of areas (x --> intake area, y--> area of throat)
% Equating stagnation pressure or intake higher than throat
% One plot: areas and chamber pressure...
%       overlay a surface that is shock stag pressure
% Two plot: area and thrust (one study per expansion) 


% Define an air mass flow rate requirement to proivde lean conditions 
% Minimum constraint is whatever stoichiometric (0.6-0.9)
% ------------ % 
% For max fuel flow rate, air mass flow can never drop below fst (requirement) 
% We need to iterate on altitude at smallest mach number to determine possible min intake stags
% For those air mass flow rates, we need to change throat areas to remain 
%   equal to or below intake stag pressure









