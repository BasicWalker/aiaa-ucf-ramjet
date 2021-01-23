%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name:  Ramjet Nozzle with the Method of Characteristics
% 
% Program Description: 
%
% 
% File Name: ramjetMOC.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Karam Paul      01/17/21  Initial Creation 
% --------------------------------------------------------------------- %

%% Setup
close all; clear; clc

%% User Defined Inputs
gamma       = 1.4;                % Cp/Cv
exitMach    = 2.4;                % Mach Number Leaving the Nozzle
n           = 7;                  % Number of Characteristic Lines
nodes       = nodeCalculator(n);   % Calculates Number of Nodes

% Calculate Total Number of Nodes Created by the Number of Characteristics


%% Initialize Variables
theta   = zeros(1, nodes);
PM      = zeros (1, nodes);
KL      = zeros(1, nodes);
KR      = zeros(1, nodes);

%% Calculate Angles at the Throat 
exitPM          = PrandtlMeyer(gamma, exitMach);
thetaMax        = exitPM / 2;
thetaChange     = thetaMax / n;

%% Begin MOC Procedure

% Store values for each point along starting characteristics 
for i = 1:n
    theta(i)    = (i)*thetaChange;
    PM(i)       = theta(i);
    KL(i)       = theta(i) - PM(i); 
    KR(i)       = theta(i) + PM(i);
end

% Store Values for the point in contact with the nozzle wall
theta(i+1)  = theta(i);
PM(i+1)     = PM(i);
KL(i+1)     = KL(i);
KR(i+1)     = KR(i);

    



    







