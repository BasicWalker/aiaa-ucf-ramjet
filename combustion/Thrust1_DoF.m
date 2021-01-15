%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  AIAA UCF Ramjet Thrust 1-DOF Script                    %
%                                                                         %
%                              Authored by                                %
%               Samer Armaly, Karam Paul, Matthew Aubertin                %
%                           January 15, 2021                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methodology
% using kinematics and Newtons 2nd law solve for thrust and acceleration of
% the ramjet.
% initial mach to design mach: using burntime find acceleration using
% kinematics V_f = V_0 + a*t then find the thrust value with the solved
% acceleration using newtons 2nd law.
% design mach onwards: using thrust value solve for acceleration then solve
% for V_f for each step using kinematics.

%% ---------- define variables ----------

% flight properties
initial_mach = 1.5;  % booster max mach
initial_altitude = 1100;  %  initial altitude for ramjet start <m>
design_mach = 2;  % mach number for criticl flight operations
burntime = 5;  % combustion burn to reach design mach <sec>


% vehicle properties
dry_mass = 4.536;  % mass of ramjet without fuelgrain <kg>
fuel_mass = 1.134;  % mass of fuel grain <kg>
wet_mass = dry_mass + fuel_mass;  % mass of ramjet without fuelgrain <kg>
fuel_mass_flow = 0.2;  % <kg/s>
c_d = 0.12;  % drag coefficient

% environment properties
g = 9.81;  % gravitaional constant <m/s^2>
initial_temp = 292.91;  % temperature at initial altitude <K>
initial_pressure = 89695;  % pressure at initial altitude <pa>
initial_density = 1.0589;  % density at initial altitude <kg/m^3>
gamma = 1.4;  % specific heat ratio
R = 287;  % <J/kg*K>

% simulation properties
step_size = 0.1;

%% ---------- simulation ----------

t = 0:step_size:burntime;  % time iteration array

% find initial acceleration required

V_0 = initial_mach*sqrt(gamma*R*initial_temp);
V_f = design_mach*sqrt(gamma*R*initial_temp);
initial_acceleration = (V_f - V_0)/burntime;

% pre-allocate array variables
thrust = zeros(1, size(t,2));
drag = zeros(1, size(t,2));
mass = zeros(1, size(t,2));
weight = zeros(1, size(t,2));

mass(1) = wet_mass;
weight(1) = g*mass(1);

for i = 2:size(t,2)
    mass() = mass






