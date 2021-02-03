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
n           = 5;                  % Number of Characteristic Lines
nodes       = nodeCalculator(n);  % Calculates Number of Nodes

% Known Nozzle Geometry Parameters
exitRadius      = 1.5;                                % Radius of nozzle exit <in>
exitArea        = pi * exitRadius^2;                % Area of nozzle exit <in2>
areaRatio       = isentropicFlow(gamma, exitMach);  % Area ratio given exit mach
throatArea      = exitArea / areaRatio;             % Area of nozzle throat <in2>
throatRadius    = sqrt(throatArea / pi);            % Radius of nozzle throat <in>

%% Initialize Variables
theta       = zeros(nodes, 1);
PM          = zeros(nodes, 1);
KL          = zeros(nodes, 1);
KR          = zeros(nodes, 1);

mach        = zeros(nodes, 1);
machAngle   = zeros(nodes, 1);

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

% Initialize at 1st right running characteristic
p = 2; 
q = n + 2;
for iter = 1:n-1 
    j = p; %2 %3 %4 %5
    h = q; %7 %12 %16 %19
    
    % Set values along centerline
    theta(h)    = 0;
    KR(h)       = KR(j);
    PM(h)       = KR(h) - theta(h);
    KL(h)       = theta(h) - PM(h);
    j           = j + 1; %3 %4 %5

        % Loop through left running characteristic except for wall
        for i = h+1:n-p+q   % from 8 to 10 % from 13 to 14 % from 17 to 17
            KR(i)    = KR(j); 
            KL(i)    = KL(i-1);
            theta(i) = 0.5 * (KR(i) + KL(i));
            PM(i)    = 0.5 * (KR(i) - KL(i));
            j = j + 1;
        end
    
    % Increment to wall. Last iteration will not go through above loop
    if i == n-p+q
        h = i + 1; %11 %15 %18
    else
        h = h + 1; %20
    end
    
    % Set values along the nozzle wall
    theta(h)    = theta(h-1);
    PM(h)       = PM(h-1);
    KL(h)       = KL(h-1);
    KR(h)       = KR(h-1);
    
    % Set up for next right running characteristic line
    p = p + 1;
    q = h + 1; 
end

% Determine Mach Number and Mach Angle at Each Node
for i = 1:nodes
    [mach(i), machAngle(i)] = inversePrandtlMeyer(gamma, PM(i));
end

% At this point, all of the value should be set, including along the nozzle
% wall. Tabulate results
T = table(theta, PM, KL, KR, mach, machAngle);

%% Plot Nozzle Wall
x = [nodes, 1];
y = [nodes, 1];
centerline_x = [n, 1];
centerline_y = [n, 1];

xlim = 4;
ylim = exitRadius; 
line

% Plots Centerline
plot([0 xlim], [0 0], '--k', 'LineWidth', 2)
hold on

%
for i = 1:n
    centerline_x(i) = throatRadius * tand(90 - PM(i) - theta(i));
    centerline_y(i) = 0;
    plot( [0 centerline_x(i)] , [throatRadius centerline_y(i)] )
    hold on
end

    
    











