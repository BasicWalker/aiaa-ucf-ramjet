% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% Program Name:  SFRJ Internal Ballistic Simulatior
% 
% Program Description: 
% First order approximation of the performance 
% capabilities of a solid fuel ramjet
% 
% File Name: StartUp.m 
% 
% File Description: 
% Main executable file for the SFRJ Internal Ballistic Simulator. Defines
% user parameters and initializes key variables
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %
close all
clear
clc
% ------------------------- Import Tables ------------------------- %

load ('N2OProperties.mat')
load ('AltitudeData.mat')

% --------------------- Initialize Variables ---------------------- %

SFRJDt          = 1/100;                    % Simulation rate (Hz) 
n               = 1;                        % Initialize counter
time            = 0.0;                      % Initialize time (s)
BurnTime(1)     = 0.0;                      % Initialize BurnTime variable (s)
StopBurn        = false;                    % Burn status flag (boolean) 
gravity         = 9.81;                     % gravitation acceleration constant (m/s^2)
In2Mtr          = 39.3701;                  % Inch to meter converstion 
Bar2kPa         = 100.0;                    % Psi to Pa conversion
C2K             = 273.15;                   % Celcius to Kelvin conversion
R               = 287;                      % Universal Gas Constant for air
AltitudeTbl     = table2array(AltitudeData);% Altitude Table

%% User Defined Parameters 
% --------------- Environmental User Defined Parameters --------------- %

Height(1)       = 1100;                     % Altitude (m) 

% --------------- Fuel Grain User Defined Parameters --------------- %

GrainOD         =  2.75 /In2Mtr;            % Grain OD (m)
GrainID(1)      =  1.00 /In2Mtr;            % Grain ID (m)
GrainL          = 15.00 /In2Mtr;            % Grain Length (m)
FuelRho         = 1020;                     % Grain Density (kg/m^3)

% ----------------- Intake User Defined Parameters ----------------- %

InltD           = 0.75 / In2Mtr;            % Diameter of inlet (m)
InltArea        = pi*InltD^2*(1/4);         % Area of inlet (m)
InltPres(1)     = 6.1493 * Bar2kPa;         % Pressure of inlet (kPa)
InltRho         = 4.5122;                   % Density of air at the inlet (kg/m^3)
InltTemp        = 250.05 + C2K;             % Temp of air at the inlet (K)
gamma           = 1.3845;                   % Specific heat ratio of air 
InltSpeedSnd    = sqrt(gamma*R*InltTemp);   % Speed of Sound at inlet (m/s)
InltMach        = 0.2;                      % Mach number at the inlet
InltVel(1)      = InltMach*InltSpeedSnd;    % Velocity of air at the inlet
InltMassFlw     = InltRho*InltVel*InltArea; % Mass flow rate of air at the inlet

% ----------------------------- Nozzle ----------------------------- %

NzlThrtDia      = 0.985 /In2Mtr;            % Throat Diameter, assuming exit area is 1.6 in diameter
NzlAT           = pi*(NzlThrtDia/2)^2;      % Throat area

% --------------- Chemistry User Defined Parameters ---------------- %




%% Main Code
Main