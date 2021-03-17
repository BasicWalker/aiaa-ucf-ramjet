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

%% User Defined Parameters 
% --------------- Environmental User Defined Parameters --------------- %

flight_mach(1)  = 2.0;                      % Booster max mach
altitude(1)     = 1100;                     % Initial altitude for ramjet start (m)
c_d             = 0.20;                     % Drag coefficient (0.35)
S               = 0.008119;                 % Frontal surface area (m^2)
gamma_atm       = 1.4;                      % Specific heat ratio
R               = 287.05;                   % Ideal gas constant (J/kg*K)
dry_mass        = 6.80389;                  % Mass of ramjet without fuelgrain (kg)

% --------------- Fuel Grain User Defined Parameters --------------- %

GrainOD         =  2.75 /In2Mtr;            % Grain OD (m)
GrainID(1)      =  1.25 /In2Mtr;            % Grain ID (m)
GrainL          = 15.00 /In2Mtr;            % Grain Length (m)
FuelRho         = 1020;                     % Grain Density (kg/m^3)

% ----------------- Intake User Defined Parameters ----------------- %

InltD           = 0.85 / In2Mtr;            % Diameter of inlet (m)
InltArea        = pi*InltD^2*(1/4);         % Area of inlet (m)
gamma_Inlt      = 1.3845;                   % Specific heat ratio of air 
Area_3          = 2.5459e-04;               % Area of throat (m^2)
radius_combustor= InltD/2;                  % Radius of the combustor inlet (m)
Area_combustor  = pi*radius_combustor^2;    % Area of the combustor inlet (m^2)
def             = 17;                       % Deflection angle (deg)
gamma           = 1.4;                      % Specific heat ratio (atm)
% mach_i          = 2.1;                      % Initial Vehicle Mach 

% ----------------------------- Nozzle ----------------------------- %

NzlThrtDia      = 0.985 /In2Mtr;            % Throat Diameter, assuming exit area is 1.6 in diameter
NzlAT           = pi*(NzlThrtDia/2)^2;      % Throat area (m^2)

% --------------- Chemistry User Defined Parameters ---------------- %

chem = Chemistry();

% ----------------- Trajectory Initial Conditions ------------------ %
                    
density(1)      = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);    % Atmosphereic Density (kg/m^3)
pressure_atm(1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);    % Atmosphereic Pressure (Pa)
pressure_atm(1) = pressure_atm(1)*(1/Pa2kPa);                               % Atmosphereic Pressure (kPa)
Temp_atm(1)     = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);       % Atmosphereic Temperature (K)
velocity(1)     = flight_mach(1)*sqrt(gamma_atm*R*Temp_atm(1));             % Atmosphereic Velocity (m/s)
drag(1)         = c_d*0.5*density(1)*velocity(1)^2*S;                       % Induced Drag (N)
Rho_atm(1)      = density(1);                                               % Density value for Intake Model - Clean up later
Thrustdlvd(1)   = 400;                                                      % Fake first Thrust value 



%% Main Code
Main