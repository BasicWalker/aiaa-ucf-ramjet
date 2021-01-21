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
burntime = 5;  % burntime to reach design mach <sec>

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
step_size = 0.1;

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
resid(1) = 10;
j = 1;
chng(1) = 1;


while resid(j) > eps('single')
    
    % initial timestep parameters
    altitude(1) = initial_altitude;
    density(1) = interp1(T.Hgtkm, T.DensMean, (altitude(1))/1e3);
    pressure(1) = interp1(T.Hgtkm, T.PresMean, (altitude(1))/1e3);
    temperature(1) = interp1(T.Hgtkm, T.Tmean, (altitude(1))/1e3);
    mach(1) = initial_mach;
    velocity(1) = mach(1)*sqrt(gamma*R*temperature(1));
    drag(1) = c_d*0.5*density(1)*velocity(1)^2*S;
    mass(1) = wet_mass;
    weight(1) = g*mass(1);
    acceleration(1) = (thrust(j,1) + drag(1) + weight(1))/ mass(1);

    for i = 2:size(t,2)
        velocity(i) = velocity(i-1) + acceleration(i-1)*(step_size);
        altitude(i) = altitude(i-1) + velocity(i-1)*step_size + 0.5*acceleration(i-1)*step_size^2;
        density(i) = interp1(T.Hgtkm, T.DensMean, (altitude(i))/1e3);
        pressure(i) = interp1(T.Hgtkm, T.PresMean, (altitude(i))/1e3);
        temperature(i) = interp1(T.Hgtkm, T.Tmean, (altitude(i))/1e3);
        mach(i) = velocity(i)/sqrt(gamma*R*temperature(i));
        drag(i) = c_d*0.5*density(i)*velocity(i)^2*S;
        if mass(i-1) > dry_mass
            mass(i) = mass(i-1) - fuel_mass_flow*step_size;
        else
            mass(i) = dry_mass;
            fprintf("out of fuel")
            break
        end
        weight(i) = g*mass(i);
        acceleration(i) = (thrust(i) + drag(i) + weight(i))/ mass(i);
    end
    final_mach(j) = mach(end);
    if j == 1
        chng(j) = 1;
    else
        chng(j) = (final_mach(j) - final_mach(j-1))/(thrust(j,1)-thrust(j-1,1));
    end
    thrust(j,:) = thrust(j-1,:) - final_mach(j-1)/chng(j-1);
    resid(j) = abs(design_mach-final_mach(j));
    
   j = j+1;
   if j > 1000
       break
   end
end



% velocity(1) = initial_mach*sqrt(gamma*R*initial_temperature);
% mach(1) = initial_mach;
% altitude(1) = initial_altitude;
% density(1) = initial_density;
% pressure(1) = initial_pressure;
% temperature(1) = initial_temperature;
% mass(1) = wet_mass;
% weight(1) = g*mass(1);
% drag(1) = c_d*0.5*design_density*velocity(1)^2;
% % thrust(1) = mass(1)*initial_acceleration + drag(1) + weight(1);
% for i = 2:size(t,2)
%     velocity(i) = velocity(i-1) + initial_acceleration*(step_size);
%     altitude(i) = altitude(i-1) + velocity(i)*step_size + 0.5*initial_acceleration*step_size^2;
%     density(i) = interp1(T.Hgtkm, T.DensMean, (altitude(i))/1e3);
%     pressure(i) = interp1(T.Hgtkm, T.PresMean, (altitude(i))/1e3);
%     temperature(i) = interp1(T.Hgtkm, T.Tmean, (altitude(i))/1e3);
%     mach(i) = velocity(i)/sqrt(gamma*R*temperature(i));
%     if mass(i-1) > dry_mass
%         mass(i) = mass(i-1) - fuel_mass_flow*step_size;
%     else
%         mass(i) = dry_mass;
%         fprintf("out of fuel")
%         break
%     end
%     weight(i) = mass(i)*g;
%     drag(i) = c_d*0.5*density(i)*velocity(i)^2;
%     thrust(i) = mass(i)*initial_acceleration + drag(i) + weight(i);
% end

% %% Plotting
% figure('Name','Thrust & Drag Vs. Time');
% plot(t,thrust); hold on;
% plot(t,drag);
% ylabel('<N>');
% xlabel('<s>');
% legend('Thrust','Drag');
% hold off;
% 
% figure('Name','Velocity & Mach Number Vs. Time')
% ax1 = subplot(2,1,1);
% plot(ax1,t,velocity);
% title(ax1,'Velocity');
% ylabel(ax1,'<m/s>');
% xlabel(ax1,'<s>');
% ax1 = subplot(2,1,2);
% plot(ax1,t,mach);
% title(ax1,'Mach Number');
% ylabel(ax1,'<unitless>');
% 
% figure('Name','Environment Vs. Time');
% ax1 = subplot(4,1,1);
% plot(ax1,t,altitude);
% title(ax1,'Altitude');
% ylabel(ax1,'<m>');
% xlabel(ax1,'<s>');
% ax2 = subplot(4,1,2);
% plot(ax2,t,temperature);
% title(ax2,'Temperature');
% ylabel(ax2,'<K>');
% ax3 = subplot(4,1,3);
% plot(ax3,t,pressure);
% title(ax3,'Pressure');
% ylabel(ax3,'<Pa>');
% ax4 = subplot(4,1,4);
% plot(ax4,t,density);
% title(ax4,'Density');
% ylabel(ax4,'<kg/m^3>');



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
