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
chem = Chemistry();

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

% Iterate on Air Mass Flow Rate
airMassFlowStep = 0.001;
min_airMassFlow = 0.5;
max_airMassFlow = 2;

% Iterate on Throat Area
throatRadiusStep = 0.001;
min_throatRadius = 0.1; % < 0.1 in >
max_throatRadius = 1;  % < 1 in >

% Array Sizes 
m = ((max_airMassFlow - min_airMassFlow) / airMassFlowStep) + 1;
n = ((max_throatRadius - min_throatRadius) / throatRadiusStep) + 1;

% Pre-allocate Resources 
airMassFlow     = zeros(m,n);
totalMassFlow   = zeros(m,n);
f               = zeros(m,n);
phi             = zeros(m,n);
T_AFT           = zeros(m,n);

throatRadius    = zeros(m,n);
throatArea      = zeros(m,n);

stagPres     = zeros(m,n);

counter = 1;
% Begin Iterations
for i = 1:m
    
    for j = 1:n
        
        % Calcuate AFT for each air mass flow rate 
        airMassFlow(i,j) = min_airMassFlow + airMassFlowStep*(i-1);
        totalMassFlow(i,j) = airMassFlow(i,j) + fuelMassFlow;
        f(i,j) = fuelMassFlow / airMassFlow(i,j);
        [phi(i,j), T_AFT(i,j)] = chem.phiSolver(f(i,j), intakeStaticTemp);
        
        % Calculate throat area in m2
        throatRadius(i,j) = min_throatRadius + throatRadiusStep*(j-1);
        throatRadius(i,j) = throatRadius(i,j) * in2m;
        throatArea(i,j) = pi * throatRadius(i,j)^2;
        
        % Calculate stagnation pressure to choke nozzle with this mass flow
        stagPres(i,j) = (totalMassFlow(i,j) * sqrt(T_AFT(i,j))) / ( throatArea(i,j) * (sqrt(gamma/R)) ...
                            * ((gamma+1)/2)^-((gamma+1)/(2*(gamma-1))));  
    end
    
    disp(counter)
    counter = counter + 1;
end
toc





