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
close all; clear; clc
if exist('T','var')==0
    load GRAM_Model.mat
end

% flight properties
initial_mach = 1.5;  % booster max mach
initial_altitude = 1100;  %  initial altitude for ramjet start <m>
design_mach = 2;  % mach number for criticl flight operations
design_altitude = 4051.8;  % <m>
burntime = 5;  % combustion burn to reach design mach <sec>


% vehicle properties
dry_mass = 4.536;  % mass of ramjet without fuelgrain <kg>
fuel_mass = 1.134;  % mass of fuel grain <kg>
wet_mass = dry_mass + fuel_mass;  % mass of ramjet without fuelgrain <kg>
fuel_mass_flow = 0.0033;  % <kg/s>
c_d = 0.012;  % drag coefficient
S = 0.008119;  % frontal surface area <m^2>


% environment properties
g = 9.81;  % gravitaional constant <m/s^2>
initial_temperature = 292.91;  % temperature at initial altitude <K>
initial_pressure = 89695;  % pressure at initial altitude <pa>
initial_density = 1.0589;  % density at initial altitude <kg/m^3>
design_temperature = interp1(T.Hgtkm, T.Tmean, (design_altitude)/1e3);
design_pressure = interp1(T.Hgtkm, T.PresMean, (design_altitude)/1e3);
design_density = interp1(T.Hgtkm, T.DensMean, (design_altitude)/1e3);
gamma = 1.4;  % specific heat ratio
R = 287;  % <J/kg*K>


% simulation properties
step_size = 0.1;

%% ---------- simulation ----------

t = 0:step_size:burntime;  % time iteration array

% find initial acceleration required

initial_acceleration = (design_mach*sqrt(gamma*R*design_temperature) - ...
    initial_mach*sqrt(gamma*R*initial_temperature))/burntime;

% pre-allocate array variables
thrust = zeros(1, size(t,2));
drag = zeros(1, size(t,2));
mass = zeros(1, size(t,2));
weight = zeros(1, size(t,2));
velocity = zeros(1, size(t,2));
mach = zeros(1, size(t,2));
altitude = zeros(1, size(t,2));
density = zeros(1, size(t,2));
pressure = zeros(1, size(t,2));
temperature = zeros(1, size(t,2));

velocity(1) = initial_mach*sqrt(gamma*R*initial_temperature);
mach(1) = initial_mach;
altitude(1) = initial_altitude;
density(1) = initial_density;
pressure(1) = initial_pressure;
temperature(1) = initial_temperature;
mass(1) = wet_mass;
weight(1) = g*mass(1);
drag(1) = c_d*S*0.5*design_density*velocity(1)^2;
thrust(1) = mass(1)*initial_acceleration + drag(1) + weight(1);

for i = 2:size(t,2)
    velocity(i) = velocity(i-1) + initial_acceleration*(step_size);
    altitude(i) = altitude(i-1) + velocity(i)*step_size + 0.5*initial_acceleration*step_size^2;
    density(i) = interp1(T.Hgtkm, T.DensMean, (altitude(i))/1e3);
    pressure(i) = interp1(T.Hgtkm, T.PresMean, (altitude(i))/1e3);
    temperature(i) = interp1(T.Hgtkm, T.Tmean, (altitude(i))/1e3);
    mach(i) = velocity(i)/sqrt(gamma*R*temperature(i));
    if mass(i-1) > dry_mass
        mass(i) = mass(i-1) - fuel_mass_flow*step_size;
    else
        mass(i) = dry_mass;
        fprintf("out of fuel")
        break
    end
    weight(i) = mass(i)*g;
    drag(i) = c_d*S*0.5*density(i)*velocity(i)^2;
    thrust(i) = mass(i)*initial_acceleration + drag(i) + weight(i);
end

%% Plotting
figure('Name','Thrust & Drag Vs. Time');
plot(t,thrust); hold on;
plot(t,drag);
ylabel('<N>');
xlabel('<s>');
legend('Thrust','Drag');
hold off;

figure('Name','Velocity & Mach Number Vs. Time')
ax1 = subplot(2,1,1);
plot(ax1,t,velocity);
title(ax1,'Velocity');
ylabel(ax1,'<m/s>');
xlabel(ax1,'<s>');
ax1 = subplot(2,1,2);
plot(ax1,t,mach);
title(ax1,'Mach Number');
ylabel(ax1,'<unitless>');

figure('Name','Environment Vs. Time');
ax1 = subplot(4,1,1);
plot(ax1,t,altitude);
title(ax1,'Altitude');
ylabel(ax1,'<m>');
xlabel(ax1,'<s>');
ax2 = subplot(4,1,2);
plot(ax2,t,temperature);
title(ax2,'Temperature');
ylabel(ax2,'<K>');
ax3 = subplot(4,1,3);
plot(ax3,t,pressure);
title(ax3,'Pressure');
ylabel(ax3,'<Pa>');
ax4 = subplot(4,1,4);
plot(ax4,t,density);
title(ax4,'Density');
ylabel(ax4,'<kg/m^3>');



% %% functions
% function 
% final_temp = initial_temp;
% res = 10;
% while res > eps('single')
%     V_0 = initial_mach*sqrt(gamma*R*initial_temp);
%     V_f = design_mach*sqrt(gamma*R*final_temp);
%     initial_acceleration = (V_f - V_0)/burntime;
% 
%     delta_x = V_0*burntime + 0.5*initial_acceleration*burntime^2;
%     design_altitude = delta_x + initial_altitude;
% 
%     new_temp = interp1(T.Hgtkm, T.Tmean, (design_altitude)/1e3);
%     res = abs(final_temp - new_temp);
%     final_temp = new_temp;
%     fprintf('%f\n', res);
% end
% 
% 
% 
% 
