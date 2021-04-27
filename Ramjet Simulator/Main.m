% ---------- SFRJ Numerical Simulation / UCF CAPSTONE PROJECT ---------- %
% Program Name:  SFRJ Internal Ballistic Simulator
%
% Program Description: 
% This program models and simulates critical performance parameters of a 
% Solid-Fuel Ramjet. 
% 
% File Description:
% Main executive file.  
%
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial creation 
% Ethan Sherlock  02/14/21  005  1DOF Trajectory Initialization
% Samer & karam   --------  ---  Sim Revamp
% Ethan Sherlock  04/14/21  ---  2DOF Trajectory Update
% ---------------------------------------------------------------------- %
clear; clc; close all;
addpath(genpath(pwd))
addpath('..\Common Resources')
tic
% Import data
load RamjetDimensions.mat  % load in the ramjet design
load GRAM_Model.mat  % GRAM atmospheric model
load Constants.mat  % load in constants and conversions

% Initialize Simulation Variables
SFRJDt              = 1/100;                                % Simulation rate (Hz) 
n                   = 1;                                    % Initialize counter
time                = 0.0;                                  % Initialize time (s)
StopBurn            = false;                                % Burn status flag (boolean) 
Burnout             = false;                                % Burn out status flag (boolean)

% User Defined Parameters 
vehicle.Mach(1)         = 2.5;  % Ramjet Initial Mach number
trajectory.Z_pos(1)     = 1000.0;  % Initial altitude for ramjet start (m)
alpha                   = 0.0;  % Launch Angle (deg) - in reference to horizon
trajectory.LiftOnOff    = 1;  % 0.0 = off; % 1.0 = on

% -------------------- Commence Ramjet Simulation  --------------------- %
Initialize  % set initial conditions
while StopBurn == 0
    BurnTime(n) = time;                                     % Simulation Time 
    if Burnout == 0
        RegressionRate                                      % Call Regression Rate Model
        GrainGeometry                                       % Call Instantaneous Grain Geometry Model
        Intake                                              % Call Intake Model
        Gas                                                 % Call Gas Model (And Chemistry Model)  
        CombustionChamber                                   % Call Combustion Chamber Model
        Nozzle                                              % Call Nozzle Model
        Thrust                                              % Call Thrust Model 
        if StopBurn == 1                                    % do not update trajectory since fuel is depleted
            [~, nozzle.idealIndex] = min(abs(nozzle.exitPressureDifference)); % Mach @ Ideal Expansion
            nozzle.idealExitMach = intake.mach(1,nozzle.idealIndex);
            break
        end
        Trajectory                                          % Call Trajectory Model
        MessageFile                                         % Call Message Outputs       
    end
    time = time + SFRJDt;                                   % Step through simulation time
    n = n + 1;                                              % Increase Index
end
clear gamma R
PlotData                                                    % Plot Simulation Results