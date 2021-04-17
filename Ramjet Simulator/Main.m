% ---------- SFRJ Numerical Simulation / UCF CAPSTONE PROJECT ---------- %
% Program Name:  SFRJ Internal Ballistic Simulatior
%
% Program Description: 
% This program models and simulates critical performance parameters of a 
% Solif-Fuel Ramjet. 
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
tic
% ------------------------- Import Tables ------------------------- %

if exist('GRAM','var')==0                                   % GRAM atmospheric model
    load GRAM_Model.mat                     
end

% --------------------- Initialize Variables ---------------------- %

SFRJDt              = 1/100;                                % Simulation rate (Hz) 
n                   = 1;                                    % Initialize counter
time                = 0.0;                                  % Initialize time (s)
BurnTime(1)         = 0.0;                                  % Initialize BurnTime variable (s)
StopBurn            = false;                                % Burn status flag (boolean) 
Burnout             = false;                                % Burn out status flag (boolean)
% ------------------        Define Constants     ------------------ %
gravity             = 9.81;                                 % gravitation acceleration constant (m/s^2)
In2Mtr              = 39.3701;                              % Inch to meter converstion 
Bar2kPa             = 100.0;                                % Bar to kPa conversion
Pa2kPa              = 1000.0;                               % Pa to Kpa
C2K                 = 273.15;                               % Celcius to Kelvin conversion
R                   = 287.05;                               % Universal Gas Constant for air
OxPercent           = 0.2314;                               % Density percentage of oxygen in air by mass
gamma               = 1.4;                                  % Specific heat ratio (atm)
k                   = 1.4;                                  % Specific heat ratio (air)

% User Defined Parameters 
% ------------------ Environmental Initialization ------------------ %
Mach_f(1)           = 2.0;                                  % Booster max mach
altitude(1)         = 1000.0;                               % Initial altitude for ramjet start (m)
c_d                 = 0.23;                                 % Drag coefficient (0.35)
S                   = 0.008119;                             % Frontal surface area (m^2)
dry_mass            = 6.80389;                              % Mass of ramjet without fuelgrain (kg)

% -------------------- Combustor  Initialization ------------------- %
combustion.InletArea = 0.0064693;%pi*combustion.InletDia^2*(1/4);  % (m)
combustion.InletDia = sqrt(4*combustion.InletArea/pi);%  <m>   
combustion.InletDiaINCH = sqrt(4*combustion.InletArea/pi)* In2Mtr;  % <in>    %1.4 / In2Mtr;  % (m)
combustion.InletGamma = 1.3845;  % Specific heat ratio of air 
combustion.ChamberArea = combustion.InletArea*1.5;  % 50% larger than inlet area     %pi*combustion.ChamberRadius^2/4; %pi*combustion.ChamberRadius^2;  % Area of combustion chamber (m^2)
combustion.ChamberDia = sqrt(4*combustion.ChamberArea/pi); %<m>            %2.75 /In2Mtr / 2;  % Radius of the combustion chamber (m)
combustion.ChamberDiaINCH = sqrt(4*combustion.ChamberArea/pi)* In2Mtr; %<in>
% -------------------- Fuel Grain Initialization ------------------- %
fuel.DiaOuter             =  combustion.ChamberDia;%2.75 /In2Mtr;                        % Grain OD (m)
fuel.DiaOuterINCH             =  combustion.ChamberDia*In2Mtr;
fuel.StepHeight(1) = 1/8 /In2Mtr;
fuel.DiaInner(1)          = 2*fuel.StepHeight(1) + combustion.InletDia;    %1.50 /In2Mtr;   
fuel.DiaInnerINCH(1)          = (2*fuel.StepHeight(1) + combustion.InletDia)* In2Mtr;    %1.50 /In2Mtr;  % Grain ID (m)
fuel.Length              = 15.00 /In2Mtr;                        % Grain Length (m)
fuel.Density             = 1020;                                 % Grain Density (kg/m^3)
fuel.PortArea(1)         = pi*(fuel.DiaInner(1)^2)*(1/4);              % Fuel Port Area (m^2)
fuel.CsxArea(1)           = pi*(fuel.DiaOuter^2)*(1/4) - fuel.PortArea(1);   % Fuel Grain Crossectional Area (m^2)
fuel.Volume(1)          = fuel.CsxArea(1) * fuel.Length;                   % Fuel Grain Volume (m^3)
fuel.SurfArea(1)           = fuel.DiaInner(1)* pi * fuel.Length;              % Fuel Grain Surface Area (m^2)
fuel.Mass(1)         = fuel.Density*fuel.Volume(1);                   % Grain fuel mass, instantaneous (kg)
Mass(1)             = dry_mass + fuel.Mass(1);               % Mass of Vehicle (Kg)

% ---------------------- Intake Initialization --------------------- %
intake.Area_enter         = 0.0007;                               % Area of throat (m^2) - Drives mass flow rate through intake
intake.DeflAngle                 = 7;                                    % Deflection angle (deg)

% ---------------------- Nozzle Initialization --------------------- %
nozzle.DiaThroat          = 1.8 /In2Mtr;                          % Throat Diameter, assuming exit area is 1.6 in diameter (from HPR), 0.985
nozzle.Area_throat  = pi*(nozzle.DiaThroat)^2/4;                  % Throat area (m^2)
nozzle.Area_exit    = pi*(2.75/In2Mtr)^2/4;                 % nozzle exit area (m^2)

% --------------- Chemistry User Defined Parameters ---------------- %
chem = Chemistry();                                         % Initialize Chemistry Model

% ------------------- Trajectory Initialization -------------------- %
thrust(1)               = 0.0;                                                                  % Thrust (N) 
alpha                   = 0.0;                                                                  % Launch Angle (deg) - in reference to horizon
trajectory.LiftOnOff    = 0.0;                                                                  % 0.0 = off; % 1.0 = on
trajectory.Lift(1)      = Mass(1)*gravity*trajectory.LiftOnOff;                                 % Lift (N)
trajectory.Rho_a(1)     = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);                % Atmospheric Density (kg/m^3)
trajectory.pressure_a(1)= interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);                % Atmospheric Pressure (Pa)
trajectory.Temp_a(1)    = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);                   % Atmospheric Temperature (K)
trajectory.Vel(1)       = Mach_f(1)*sqrt(k*R*trajectory.Temp_a(1));                             % Velocity (m/s)
trajectory.Vel_x(1)     = trajectory.Vel(1)*cosd(alpha);                                        % Velocity X (m/s)
trajectory.Vel_z(1)     = trajectory.Vel(1)*sind(alpha);                                        % Velocity Z (m/s)
trajectory.F_d(1)       = 0.5*trajectory.Rho_a(1)*c_d*S*trajectory.Vel(1)^2;                    % Drag (N)
trajectory.F_dx(1)      = 0.5*trajectory.Rho_a(1)*c_d*S*trajectory.Vel(1)*trajectory.Vel_x(1);  % Drag X (N)
trajectory.F_dz(1)      = 0.5*trajectory.Rho_a(1)*c_d*S*trajectory.Vel(1)*trajectory.Vel_z(1);  % Drag Z (N)
trajectory.F_tx(1)      = thrust(1)*cosd(alpha);                                            % Thrust X (N)
trajectory.F_tz(1)      = thrust(1)*sind(alpha);                                            % Thrust Z (N)
trajectory.F_x(1)       = trajectory.F_tx(1)-trajectory.F_dx(1);                                    % Force X (N)
trajectory.F_z(1)       = trajectory.F_tz(1)-trajectory.F_dz(1)-Mass(1)*gravity+trajectory.Lift(1); % Force Z (N)
trajectory.F_net(1)     = sqrt(trajectory.F_x(n)^2+trajectory.F_z(n)^2);                        % Force (N)
trajectory.Acc(1)       = trajectory.F_net(1)/Mass(1);                                          % Acceleration (m/s/s)
trajectory.Acc_x(1)     = trajectory.F_x(1)/Mass(1);                                            % Acceleration X (m/s/s)
trajectory.Acc_z(1)     = trajectory.F_z(1)/Mass(1);                                            % Acceleration Z (m/s/s)
trajectory.X_pos(1)     = 0.0;                                                                  % Initial X Position (m)
trajectory.Z_pos(1)     = altitude(1);                                                          % Initial Z Position (m)
trajectory.Vel_t(1)     = sqrt((2*Mass(1)*gravity)/(c_d*trajectory.Rho_a(1)*S));                % Terminal Velocity (m/s)

% -------------------- Commence Ramjet Simulation  --------------------- %

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
        if StopBurn == 1  % do not update trajectory since fuel is depleted
            break
        end
        Trajectory                                          % Call Trajectory Model
        MessageFile                                         % Call Message Outputs       
    end
    time = time + SFRJDt;                                   % Step through simulation time
    n = n + 1;                                              % Increase Index
end
PlotData                                                    % Plot Simulation Results