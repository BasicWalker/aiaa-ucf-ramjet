% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Main.m 
%
% File Description: 
% Main executive model. controls logic of the program
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial creation 
% Ethan Sherlock  02/14/21  001  Chamber pressure calculation update
% Ethan Sherlock  02/14/21  005  1DOF trajectory update
% Ethan Sherlock  03/17/21  ---  Code Clean Up
% ---------------------------------------------------------------------- %

% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% Program Name:  SFRJ Internal Ballistic Simulatior
% 
% Program Description: 
% First order approximation of the performance 
% capabilities of a solid fuel ramjet
% 
% File Name: StartUp_Ramjet.m 
% 
% File Description: 
% Main executable file for the SFRJ Internal Ballistic Simulator. Defines
% user parameters and initializes key variables
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  005  Add 1DOF Trajectory Initialization
% Ethan Sherlock  03/12/21  ---  Add Intake Initialization  
% 
% ---------------------------------------------------------------------- %
close all; clear; clc;
% ------------------------- Import Tables ------------------------- %

if exist('T','var')==0                      % GRAM atmosphereic model
    load GRAM_Model.mat                     
end
GRAM = T; 

% --------------------- Initialize Variables ---------------------- %

SFRJDt          = 1/100;                    % Simulation rate (Hz) 
n               = 1;                        % Initialize counter
time            = 0.0;                      % Initialize time (s)
BurnTime(1)     = 0.0;                      % Initialize BurnTime variable (s)
StopBurn        = false;                    % Burn status flag (boolean) 
gravity         = 9.81;                     % gravitation acceleration constant (m/s^2)
In2Mtr          = 39.3701;                  % Inch to meter converstion 
Bar2kPa         = 100.0;                    % Bar to kPa conversion
Pa2kPa          = 1000.0;                   % Pa to Kpa
C2K             = 273.15;                   % Celcius to Kelvin conversion
R               = 287.05;                   % Universal Gas Constant for air
OxPercent       = 0.2314;                   % Density percentage of oxygen in air by mass

%% User Defined Parameters 
% --------------- Environmental User Defined Parameters --------------- %

flight_mach(1)  = 2.5;                      % Booster max mach
altitude(1)     = 3100;                     % Initial altitude for ramjet start (m)
c_d             = 0.23;                     % Drag coefficient (0.35)
S               = 0.008119;                 % Frontal surface area (m^2)
gamma_atm       = 1.4;                      % Specific heat ratio
dry_mass        = 6.80389;                  % Mass of ramjet without fuelgrain (kg)

% --------------- Fuel Grain User Defined Parameters --------------- %

GrainOD         =  2.75 /In2Mtr;                        % Grain OD (m)
GrainID(1)      =  1.50 /In2Mtr;                        % Grain ID (m)
GrainL          = 15.00 /In2Mtr;                        % Grain Length (m)
FuelRho         = 1020;                                 % Grain Density (kg/m^3)
PortArea(1)     = pi*(GrainID(1)^2)*(1/4);              % Fuel Port Area (m^2)
FuelCS(1)       = pi*(GrainOD^2)*(1/4) - PortArea(1);   % Fuel Grain Crossectional Area (m^2)
FuelVol(1)      = FuelCS(1) * GrainL;                   % Fuel Grain Volume (m^3)
FuelSA(1)       = GrainID(1)* pi * GrainL;              % Fuel Grain Surface Area (m^2)
FuelMass(1)     = FuelRho*FuelVol(1);                   % Grain fuel mass, instantaneous (kg)

% ----------------- Intake User Defined Parameters ----------------- %

InltD           = 1.00 / In2Mtr;            % Diameter of Combustor inlet (m)
InltArea        = pi*InltD^2*(1/4);         % Area of inlet (m)
gamma_Inlt      = 1.3845;                   % Specific heat ratio of air 
Area_3          = 5.5061e-04;               % Area of throat (m^2) - Drives mass flow rate through intake
radius_combustor= InltD/2;                  % Radius of the combustor inlet (m)
Area_combustor  = pi*radius_combustor^2;    % Area of the combustor inlet (m^2)
def             = 10;                       % Deflection angle (deg)
gamma           = 1.4;                      % Specific heat ratio (atm)

% ----------------------------- Nozzle ----------------------------- %

NzlThrtDia      = 1.8 /In2Mtr;              % Throat Diameter, assuming exit area is 1.6 in diameter (from HPR), 0.985
NzlAT           = pi*(NzlThrtDia/2)^2;      % Throat area (m^2)
NzlARatio       = 1.6;                      % Nozzle expansion ratio

% --------------- Chemistry User Defined Parameters ---------------- %

chem = Chemistry();

% ----------------- Trajectory Initial Conditions ------------------ %
                    
Rho_atm(1)      = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);    % Atmospheric Density (kg/m^3)
pressure_atm(1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);    % Atmospheric Pressure (Pa)
pressure_atm(1) = pressure_atm(1)*(1/Pa2kPa);                               % Atmospheric Pressure (kPa)
Temp_atm(1)     = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);       % Atmospheric Temperature (K)
velocity(1)     = flight_mach(1)*sqrt(gamma_atm*R*Temp_atm(1));             % Atmospheric Velocity (m/s)
drag(1)         = c_d*0.5*Rho_atm(1)*velocity(1)^2*S;                       % Induced Drag (N)
Thrustdlvd(1)   = 0.0;                                                      % Initialize First Thrust Value 
mass(1)         = dry_mass + FuelMass(1);                                   % Mass of Vehicle
weight(1)       = gravity*mass(1);                                          % Weight of Vehicle
acceleration(1) = (Thrustdlvd(1) - drag(1) - weight(1))/ mass(1);           % Initial Acceleration

% ----------------- Commence Ramjet Simulation  -------------------- %

while StopBurn == 0
    BurnTime(n) = time;                         % Simulation Time
    
    RegressionRate                              % Call Regression Rate Model
    GrainGeometry                               % Call Instantaneous Grain Geometry Model
    Trajectory                                  % Call Trajectory Model
    Intake                                      % Call Intake Model
    Gas                                         % Call Gas Model (And Chemistry Model)
    BoundaryLayer                               % Call Boundary Layer Model
    Thrust                                      % Call Thrust Model
   
    time = time + SFRJDt;                       % Step through simulation time
    n = n + 1;
end
Nozzle                                          % Call Nozzle Model
PlotData                                        % Plot Simulation Results