close all;
clc;
clearvars -except airProp

if(exist('airProp') == 0)
    airProp = [];
end
if(isempty(airProp))
    airProp = xlsread("Atmosphere_100km.xlsx");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Initial Conditions                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial Altitude;      (m), altitude at which the ramjet engine starts
start_alt = 1100; 
%Initial Velocity;      (m/s), upwards velocity at which the ramjet engine starts
start_vel = 408;
%Initial Air Density;   (kg/m^3), Air Density at the specified intial altitude
start_den = interp1(airProp(:, 1), airProp(:, 5), start_alt);
%Initial Air Pressure;  (Pa), Air Pressure at the specified initial altitude
start_pres = interp1(airProp(:, 1), airProp(:, 6), start_alt);
%Initial Air Temp;      (C), Air Temperatire at the specified altitude
start_temp = interp1(airProp(:, 1), airProp(:, 7), start_alt);

%Air Properties
gamma = 1.4;
%Gas Constant of Air;   (m^2/s^2/K)
air_r = 287.1;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%       VEHICLE PROPERTIES                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MASS PROPERTIES
%Vehicle Mass;          (kg), initial mass of the ramjet vehicle
veh_mass_wet = 5;

%FUEL PROPERTIES
%Fuel Mass;             (kg), mass of the fuel onboard the vehicle that is burned
veh_fuel_mass = 0.5;
%Burn Time;             (s), total time for which the fuel is burning
veh_burn_time = 14;
%Fuel Flow Rate;        (kg/s), flow rate of the fuel out of the vehicle
veh_fuel_rate = veh_fuel_mass / veh_burn_time;

%DIMENSIONAL PROPERTIES
%Vehicle Diameter;      (in), max diameter of the ramjet vehicle in inches
veh_dia_in = 3;
veh_dia = veh_dia_in * 2.54/100;    %convert Vehicle Diameter to meters
veh_A_frontal = pi * veh_dia^2 / 4;
veh_C_d = 0.12;
veh_tumblingDragFactor = 16.3;

%PROPULSION SYSTEM PROPERTIES
%Thrust;                (N), Thrust of the Ramjet Vehicle
veh_thrust = 200;
%Mass Flow Rate         (kg/s), Mass Flow Rate of Propellants through Nozzle
prop_massFlow_design = 1;
%Exit Velocity;         (m/s), Velocity of Combustion Products at Nozzle Exit
prop_exitVelocity_design = 800;
%Nozzle Exit Area       (m^2), Cross-Sectional Area of the Nozzle Exit
prop_exitArea_design = 4;
%Nozzle Exit Pressure   (Pa), The Pressure of the Combustion Products at Nozzle Exit
prop_exitPres_design = 49000;
%Combustion Design Pres;(Pa), The Design Pressure of the Combustion Chamber during Firing
prop_combPres_design = 1500e3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%       SIMULATION PROPERTIES                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ADVANCED_SIM = true;

i = 1;
dt = 0.05;
height(i) = start_alt;
vel(i) = start_vel;
vel_mach(i) = start_vel/sqrt(gamma*air_r*start_temp);
time_flight = (i-1)*dt;

%Max Vals:
max_height = 0;
max_vel = 0;

while height(i) > 0
    
    %Air Properties Interpolation
    air_dens(i) = airProp(floor(height(i)), 5);  %interp1(airProp(:, 1), airProp(:, 5), height(i));
    air_pres(i) = airProp(floor(height(i)), 6);  %interp1(airProp(:, 1), airProp(:, 6), height(i));
    air_temp(i) = airProp(floor(height(i)), 7);  %interp1(airProp(:, 1), airProp(:, 7), height(i));
    
    %Free-Stream Air Mass Flow Rate
    air_freeStream_flowRate(i) = air_dens(i) * vel(i) * veh_A_frontal;
    
    
    %Mass Caluclations
    veh_mass(i) = veh_mass_wet-veh_fuel_mass/veh_burn_time * time_flight;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   BURNING CALCULATIONS    %
    if(time_flight < veh_burn_time)
        
        %Mass Flow Rate Calculations
        prop_prodFlowRate(i) = air_freeStream_flowRate(i) + veh_fuel_rate;
        
        
        %Inlet Calculations
        %Section 1: Free-Stream/Inlet Tip
        
        eng_1_a = sqrt(gamma*air_r*air_temp(i));
        eng_1_vel = vel(i);
        
        [m_i, t_i, p_i, r_i, a_i] = flowisentropic(1.4, vel_mach(i), 'mach');
        
        eng_1_Pstag(i) = air_pres(i) / p_i;
        eng_1_Tstag(i) = air_temp(i) / t_i;
        
        
        %Combustion Chamber Calculations
        %Sections 4 to 5
        %Section 4: Combustor Inlet
        %Section 5: Combustor Outlet
        
        %All of these values below are assumed
        eng_4_mach = 0.2;      % Mach Number of Incoming Air at 4
        eng_4_temp = 420;      %(K); Temperature of Incoming Air at 4
        eng_4_pres = 211e3;    %(Pa); Pressure of Incoming Air at 4
        eng_4_airFlowRate = prop_prodFlowRate(i);
        %Area-Mach Relation at Section 4
        eng_4_areaMach = areaMach(eng_4_mach, gamma);
        
        
        %Convective Heat Transfer from Fuel Grain to the Air to Determine T5
        %Q_dot = h * A * (T_burning - T_air)
        h = 423;
        %A = pi * grain_r_inner^2 * grain_length;
        
        
        eng_5_temp = 1000;
        eng_5_Tstag = 1260;
        eng_5_Pstag = 168e3;
        
        
        [m_i, t_i, p_i, r_i, a_i] = flowisentropic(1.4, eng_4_mach, 'mach');
        eng_4_Tstag(i) = eng_4_temp / t_i;
        %Rayleigh Flow Ratios
        [m_r, t_r, p_r, r_i, v_i, T0_r, P0_r] = flowrayleigh(gamma, eng_4_mach, 'mach');
        
        
        %Converging Diverging Nozzle Calculations
        %Section 5, 6, 7
        %Section 5: Combustor Outlet
        %Section 6: Nozzle Throat
        %Section 7: Nozzle Exit
        
        nozz_PresRatio = prop_exitPres_design/eng_5_Pstag;        %Estimation for right now
        [m_i, t_i, p_i, r_i, a_i] = flowisentropic(gamma, nozz_PresRatio, 'pres');
        
        
        %Thrust and Drag Calculations
        Thrust(i) = veh_thrust;
        veh_mass(i) = veh_mass_wet-veh_fuel_rate * time_flight;
        burnout_i = i;
        burnout_time = time_flight+dt;
        
    else
        Thrust(i) = 0;
        veh_mass(i) = veh_mass_wet-veh_fuel_mass;
    end
    
    %C_d(i) = veh_C_d;
    
    Drag(i) = veh_C_d * air_dens(i) * vel(i)^2 * veh_A_frontal * 0.5;
    if(vel(i) < 0)
        Drag(i)  = Drag(i) * -veh_tumblingDragFactor;
    end
    
    %New Height & Velocity Calculations
    Force_y(i) = Thrust(i) - Drag(i) - 9.81*veh_mass(i);
    accel_y(i) = Force_y(i)/veh_mass(i);
    vel(i+1) = accel_y(i)*dt + vel(i);
    height(i+1) = 0.5*accel_y(i)*dt^2 + vel(i) * dt + height(i);
    
    %Mach Number Calculations
    air_temp_Kelvin(i) = air_temp(i);
    air_mach(i) = sqrt(gamma * air_r * air_temp_Kelvin(i));
    vel_mach(i+1) = vel(i)/air_mach(i);
    
    %Exit Pressure Optimization Calculation - %%%%%%%DEFUNCT%%%%%
    if(time_flight < veh_burn_time)
        air_presWeighted(i) = air_pres(i) * (vel(i)*dt);
    end
    
    
    %Maximum Values Saving
    if(height(i) > max_height)
        max_height = height(i);
    end
    if(vel(i) > max_vel)
        max_vel = vel(i);
    end
    
    
    %Iterating to next
    i = i + 1;
    time_flight = (i-1)*dt;
    time(i) = time_flight;
    i
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%       RESULTS VISUALIZATIONS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Height & Velocity vs. Time (Full Duration)')
clf(1)
t1 = title('Height & Velocity vs. Time (Full Duration)');
x1 = xlabel('Time (s)');
yyaxis left;
y1 = ylabel('Height (m)');
hold on;
plot(time(:), height(:));
yyaxis right;
y1 = ylabel('Velocity (m/s');
plot(time(:), vel(:));
hold off;

figure('Name', 'Height & Velocity vs. Time (Thrust Duration)')
clf(2)
t1 = title('Height & Velocity vs. Time (Thrust Duration)');
x1 = xlabel('Time (s)');
yyaxis left;
y1 = ylabel('Height (m)');
hold on;
plot(time(1:burnout_i), height(1:burnout_i));
yyaxis right;
y1 = ylabel('Velocity (m/s)');
plot(time(1:burnout_i), vel(1:burnout_i));
hold off;

figure('Name', 'Mach vs. Time (Thrust Duration)')
clf(3)
plot(time(1:burnout_i), vel_mach(1:burnout_i));
t1 = title('Velocity (Mach) vs. Time (Thrust Duration)');
x1 = xlabel('Time (s)');
y1 = ylabel('Mach Number');

figure('Name', 'Air Mass Flow Rate vs. Time (Thrust Duration)')
clf(4)
hold on;
yyaxis left;
y1 = ylabel('Mass Flow Rate (kg/s)');
plot(time(1:burnout_i), air_freeStream_flowRate(1:burnout_i));
yyaxis right;
y1 = ylabel('Air Density (kg/m^3)');
plot(time(1:burnout_i), air_dens(1:burnout_i));
hold off;
t1 = title('Air Mass Flow Rate vs. Time (Thrust Duration)');

figure('Name', 'Air Pressure & Density vs. Time (Thrust Duration)')
clf(5)
x1 = xlabel('Height (m)');
yyaxis left;
y1 = ylabel('Pressure (Pa)');
hold on
plot(height(1:burnout_i), air_pres(1:burnout_i));
yyaxis right;
y1 = ylabel('Density (kg/m^3)');
plot(height(1:burnout_i), air_dens(1:burnout_i));
hold off;
t1 = title('Air Pressure & Density vs. Time (Thrust Duration)');

figure('Name', 'Free-Stream Stagnation Pressure vs. Time (Thrust Duration)')
clf(6)
x1 = xlabel('Time (s)');
hold on;
yyaxis left;
y1 = ylabel('Stagnation Pressure (Pa)');
plot(time(1:burnout_i), eng_1_Pstag(1:burnout_i));
yyaxis right;
y1 = ylabel('Ambient Pressure (Pa)');
plot(time(1:burnout_i), air_pres(1:burnout_i));
t1 = title('Free-Stream Stagnation Pressure vs. Time (Thrust Duration)');

figure('Name', 'Free-Stream Stagnation Temperature vs. Time (Thrust Duration)')
clf(7)
plot(time(1:burnout_i), eng_1_Tstag(1: burnout_i));
t1 = title('Free-Stream Stagnation Temperature vs. Time (Thrust Duration)');
x1 = xlabel('Time (s)');
y1 = ylabel('Stagnation Temperature (K)');

figure('Name', 'Aerodynamic Drag vs. Time (Thrust Duration)')
clf(8)
hold on;
xlabel('Time (s)');
yyaxis left;
y1 = ylabel('Aerodynamic Drag (N)');
plot(time(1:burnout_i), Drag(1:burnout_i));
yyaxis right;
y1 = ylabel('Mach Number');
plot(time(1:burnout_i), vel_mach(1:burnout_i));
hold off;
t1 = title('Aerodynamic Drag vs. Time (Thrust Duration)');

figure('Name', 'Aerodynamic Drag vs. Time (Full Duration)')
clf(9)
plot(time(1:length(Drag)), Drag(:));
t1 = title('Aerodynamic Drag vs. Time (Full Duration)');

figure('Name', 'Local Speed of Sound vs. Altitude')
clf(10)
plot(height(1:burnout_i), air_mach(1:burnout_i));
t1 = title('Local Speed of Sound vs. Altitude');
x1 = xlabel('Height (m)');
y1 = ylabel('Speed of Sound (m/s)');

figure('Name', 'Ambient Air Temperature vs. Height');
clf(11);
plot(height(1:burnout_i), air_temp(1:burnout_i));
x1 = xlabel('Heihgt (m)');
y1 = ylabel('Temperature (K)');
t1 = title('Ambient Air Temperature vs. Height');

figure('Name', 'Drag vs. Velocity');
clf(12);
plot(vel_mach(1:burnout_i), Drag(1:burnout_i));

format longG

height(i-1)
max_height
max_height_ft = max_height*3.28
max_vel
avg_machBurning = mean(vel_mach(1,1:burnout_i))
avg_airPresBurning = mean(air_pres(1, 1:burnout_i))
max_propFlow = max(prop_prodFlowRate)
avg_propFlow = mean(prop_prodFlowRate)

%Exit Pressure Optimization Calculation





function r = areaMach(M, gamma)
    %Area Mach Function
    %Provides connection between local channel cross-sectional area and Mach number
    
    %a^b * M/(c^b)
    
    a = (gamma + 1) / 2;
    b = (gamma + 1) / (2 * (gamma - 1));
    c = 1 + (gamma - 1)/2 * M^2;
    
    r = a^b * M/(c^b);
end



