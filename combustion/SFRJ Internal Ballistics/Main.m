% ---------- SFRJ Numerical Simulation / UCF CAPSTONE PROJECT ---------- %
% Program Name:  SFRJ Internal Ballistic Simulatior
%
% Program Description: 
% This program numerically simulates the performance capabilities of a
% supersonic airbreathing solid fuel ramjet. First order approximation. 
% 
% File Description:
% Main executive file.  
%
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial creation 
% Ethan Sherlock  02/14/21  005  1DOF Trajectory Initialization
% Ethan Sherlock  04/03/21  ---  2DOF Trajectory Update
% ---------------------------------------------------------------------- %
clear; clc; close all;
tic
% ------------------------- Import Tables ------------------------- %

if exist('T','var')==0                      % GRAM atmosphereic model
    load GRAM_Model.mat
    GRAM = T; 
end

% --------------------- Initialize Variables ---------------------- %

SFRJDt          = 1/100;                    % Simulation rate (Hz) 
n               = 1;                        % Initialize counter
time            = 0.0;                      % Initialize time (s)
BurnTime        = zeros(1,1500);            % Preallocate BurnTime array (s)
StopBurn        = false;                    % Burn status flag (boolean) 
Burnout         = false;                    % Burn out status flag (boolean)
gravity         = 9.81;                     % gravitation acceleration constant (m/s^2)
In2Mtr          = 39.3701;                  % Inch to meter converstion 
Bar2kPa         = 100.0;                    % Bar to kPa conversion
Pa2kPa          = 1000.0;                   % Pa to Kpa
C2K             = 273.15;                   % Celcius to Kelvin conversion
R               = 287.05;                   % Universal Gas Constant for air
OxPercent       = 0.2314;                   % Density percentage of oxygen in air by mass

%% User Defined Parameters 
% --------------- Environmental User Defined Parameters --------------- %

Mach_f(1)       = 2.0;                      % Booster max mach
altitude(1)     = 1000.0;                   % Initial altitude for ramjet start (m)
c_d             = 0.23;                     % Drag coefficient (0.35)
S               = 0.008119;                 % Frontal surface area (m^2)
k               = 1.4;                      % Specific heat ratio (air)
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

InltD           = 1.00 / In2Mtr;                        % Diameter of Combustor inlet (m)
InltArea        = pi*InltD^2*(1/4);                     % Area of inlet (m)
gamma_Inlt      = 1.3845;                               % Specific heat ratio of air 
Area_3          = 5.5061e-04;                           % Area of throat (m^2) - Drives mass flow rate through intake
radius_combustor= InltD/2;                              % Radius of the combustor inlet (m)
Area_combustor  = pi*radius_combustor^2;                % Area of the combustor inlet (m^2)
def             = 10;                                   % Deflection angle (deg)
gamma           = 1.4;                                  % Specific heat ratio (atm)

% ----------------------------- Nozzle ----------------------------- %

NzlThrtDia      = 1.8 /In2Mtr;                          % Throat Diameter, assuming exit area is 1.6 in diameter (from HPR), 0.985
NzlAT           = pi*(NzlThrtDia/2)^2;                  % Throat area (m^2)
NzlARatio       = 1.6;                                  % Nozzle expansion ratio

% --------------- Chemistry User Defined Parameters ---------------- %

chem = Chemistry();                                     % Initialize Chemistry Model

% ----------------- Trajectory Initial Conditions ------------------ %

alpha           = 0.0;                                                      % Launch Angle (deg) - in reference to horizon
Rho_a(1)        = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);    % Atmospheric Density (kg/m^3)
pressure_atm(1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);    % Atmospheric Pressure (Pa)
pressure_atm(1) = pressure_atm(1)/Pa2kPa;                                   % Atmospheric Pressure (kPa)
Temp_a(1)       = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);       % Atmospheric Temperature (K)
Mass(1)         = dry_mass + FuelMass(1);                                   % Mass of Vehicle (Kg)
Vel(1)          = Mach_f(1)*sqrt(k*R*Temp_a(1));                            % Velocity (m/s)
Vel_x(1)        = Vel(1)*cosd(alpha);                                       % Velocity X (m/s)
Vel_z(1)        = Vel(1)*sind(alpha);                                       % Velocity Z (m/s)
F_d(1)          = 0.5*Rho_a(1)*c_d*S*Vel(1)^2;                              % Drag (N)
F_dx(1)         = 0.5*Rho_a(1)*c_d*S*Vel(1)*Vel_x(1);                       % Drag X (N)
F_dz(1)         = 0.5*Rho_a(1)*c_d*S*Vel(1)*Vel_z(1);                       % Drag Z (N)
F_t(1)          = 0.0;                                                      % Thrust (N) 
F_tx(1)         = F_t(1)*cosd(alpha);                                       % Thrust X (N)
F_tz(1)         = F_t(1)*sind(alpha);                                       % Thrust Z (N)
F_x(1)          = F_tx(1) - F_dx(1);                                        % Force X (N)
F_z(1)          = F_tz(1) - F_dz(1);% - Mass(1)*gravity;                      % Force Z (N)
F_net(1)        = sqrt(F_x(n)^2 + F_z(n)^2);                                % Force (N)
Acc(1)          = F_net(1)/Mass(1);                                         % Acceleration (m/s/s)
Acc_x(1)        = F_x(1)/Mass(1);                                           % Acceleration X (m/s/s)
Acc_z(1)        = F_z(1)/Mass(1);                                           % Acceleration Z (m/s/s)
X_pos(1)        = 0.0;                                                      % Initial X Position (m)
Z_pos(1)        = altitude(1);                                              % Initial Z Position (m)
Vel_t(1)        = sqrt((2*Mass(1)*gravity)/(c_d*Rho_a(1)*S));               % Terminal Velocity (m/s)


while StopBurn == 0
    BurnTime(n) = time;     % Simulation Time
    if Burnout == 0
        RegressionRate      % Call Regression Rate Model
        GrainGeometry       % Call Instantaneous Grain Geometry Model
        Trajectory          % Call Trajectory Model
        Intake              % Call Intake Model
        Gas                 % Call Gas Model (And Chemistry Model)
        BoundaryLayer       % Call Boundary Layer Model
        Thrust              % Call Thrust Model
        MessageFile         % Call Simulation Message Outputs
    else
        Trajectory          % Call Trajectory Model
        MessageFile         % Call Simulation Message Outputs
    end
    time = time + SFRJDt;   % Step through simulation time
    n = n + 1;
end
Nozzle                      % Call Nozzle Model
PlotData                    % Plot Simulation Results