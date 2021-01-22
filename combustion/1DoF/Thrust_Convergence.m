%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  AIAA UCF Ramjet Thrust 1-DOF Script                    %
%                                                                         %
%                              Authored by                                %
%               Samer Armaly, Karam Paul, Matthew Aubertin                %
%                           January 15, 2021                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methodology

%% ---------- define variables ----------
close all; clear; clc
if exist('T','var')==0
    load GRAM_Model.mat
end

% flight properties
initial_mach = 1.2;  % booster max mach
initial_altitude = 1100;  %  initial altitude for ramjet start <m>
design_mach = 2;  % mach number for criticl flight operations
burntime = 3;  % burntime to reach design mach <sec>

% vehicle properties
dry_mass = 4.536;  % mass of ramjet without fuelgrain <kg>
fuel_mass = 1.134;  % mass of fuel grain <kg>
wet_mass = dry_mass + fuel_mass;  % mass of ramjet without fuelgrain <kg>
fuel_mass_flow = 0.0033;  % <kg/s>
c_d = 0.012;  % drag coefficient
S = 0.008119;  % frontal surface area <m^2>

% environment properties
g = 9.81;  % gravitaional constant <m/s^2>
gamma = 1.4;  % specific heat ratio
R = 287;  % <J/kg*K>

% simulation properties
step_size = 0.01;

%% ---------- simulation ----------

t = 0:step_size:burntime;  % time iteration array

% pre-allocate array variables
thrust = zeros(1, size(t,2)) + 200;
drag = zeros(1, size(t,2));
mass = zeros(1, size(t,2));
weight = zeros(1, size(t,2));
velocity = zeros(1, size(t,2));
mach = zeros(1, size(t,2));
altitude = zeros(1, size(t,2));
density = zeros(1, size(t,2));
pressure = zeros(1, size(t,2));
temperature = zeros(1, size(t,2));
acceleration = zeros(1, size(t,2));

% iteration parameters
j = 1;

while true
    % initial timestep parameters
    altitude(j,1) = initial_altitude;
    density(j,1) = interp1(T.Hgtkm, T.DensMean, (altitude(j,1))/1e3);
    pressure(j,1) = interp1(T.Hgtkm, T.PresMean, (altitude(j,1))/1e3);
    temperature(j,1) = interp1(T.Hgtkm, T.Tmean, (altitude(j,1))/1e3);
    mach(j,1) = initial_mach;
    velocity(j,1) = mach(j,1)*sqrt(gamma*R*temperature(j,1));
    drag(j,1) = c_d*0.5*density(j,1)*velocity(j,1)^2*S;
    mass(j,1) = wet_mass;
    weight(j,1) = g*mass(j,1);
    acceleration(j,1) = (thrust(j,1) + drag(j,1) + weight(j,1))/ mass(j,1);

    for i = 2:size(t,2)
        velocity(j,i) = velocity(j,i-1) + acceleration(j,i-1)*(step_size);
        altitude(j,i) = altitude(j,i-1) + velocity(j,i-1)*step_size + 0.5*acceleration(j,i-1)*step_size^2;
        density(j,i) = interp1(T.Hgtkm, T.DensMean, (altitude(j,i))/1e3);
        pressure(j,i) = interp1(T.Hgtkm, T.PresMean, (altitude(j,i))/1e3);
        temperature(j,i) = interp1(T.Hgtkm, T.Tmean, (altitude(j,i))/1e3);
        mach(j,i) = velocity(j,i)/sqrt(gamma*R*temperature(j,i));
        drag(j,i) = c_d*0.5*density(j,i)*velocity(j,i)^2*S;
        if mass(j,i-1) > dry_mass
            mass(j,i) = mass(j,i-1) - fuel_mass_flow*step_size;
        else
            mass(j,i) = dry_mass;
            fprintf("out of fuel")
            break
        end
        weight(j,i) = g*mass(j,i);
        acceleration(j,i) = (thrust(j,i) + drag(j,i) + weight(j,i))/ mass(j,i);
    end
    final_mach(j) = mach(end);
    resid(j) = abs(design_mach-final_mach(j));
    if resid(j) < eps('single')
        break
    end
    if j == 1
        chng(j) = 1;
        thrust(j+1,:) = 201;
    else
        chng(j) = (resid(j-1) - resid(j))/(thrust(j-1,1)-thrust(j,1));
        thrust(j+1,:) = thrust(j,:) - resid(j)/chng(j);
    end
    j = j+1;
    if j > 1000
        break
    end 
end

fprintf("%4.2f N of thrust required for %i seconds to achieve mach %4.2f", thrust(end,end), burntime, final_mach(end))

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

