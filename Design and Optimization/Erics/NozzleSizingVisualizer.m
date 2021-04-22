clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%
%      Constraints
%%%%%%%%%%%%%%%%%%%%%%%
constrain = 1;
constrain_FlowSep = 0.4;
constraint_ThroatTemp = 3200; %Kelvin

chamberArea = 0.008107319665560;

%%%%%%%%%%%%%%%%%%%%%%%
%   Prop Properties
%%%%%%%%%%%%%%%%%%%%%%%
gamma = 1.4;
mw = 420;

Ae_input = input('Enter Nozzle Exit Diameter (in):   ');
exitArea_max = (Ae_input*2.54/1e2)^2 * pi / 4;

%Minimum Factor Exit Pressure can be before flow seperation occurs
constrain_FlowSep = 0.4;
%Ambient Pressure
ambientPressure = 89.7e3;

%Pascal
chamberPres_start = 200e3;
chamberPres_end = 600e3;
chamberPres_step = 10e3;

%Kelvin
chamberTemp_start = 1500;
chamberTemp_end = 3600;
chamberTemp_step = 50;

%Mach
chamberVel_start = 100;
chamberVel_end = 600;
chamberVel_step = 10;

i = int64((chamberVel_end - chamberVel_start) / chamberVel_step);
j = (chamberTemp_end - chamberTemp_start) / chamberTemp_step;
k = (chamberPres_end - chamberPres_start) / chamberPres_step;

data = zeros(i, j, k, 11);

%%%%%%%%%%%%%%%%%%%%%%%
%   Data Structure
%%%%%%%%%%%%%%%%%%%%%%%

%   data(x, y, z, p)
%   
%   x = Flow Velocity at (5)
%   y = Static Flow Temperature at (5)
%   z = Static Pressure at (5)
%
%   Parameters:
%
%   1   |   Throat Area
%   2   |   Throat Pressure (Static)
%   3   |   Throat Temperature (Static)
%   4   |   Inlet Flow Velocity
%   5   |   Inlet Flow Temperature (Static)
%   6   |   Expansion Ratio (Supersonic)
%   7   |   Exit Mach
%   8   |   Exit Pressure (Static)
%   9   |   Exit Velocity
%   10  |   Exit Temperature (Static)
%   11  |   
%

throatArea = zeros(i, j, k);
throatPres = zeros(i, j, k);
throatTemp = zeros(i, j, k);

i = 1;
j = 1;
k = 1;
c = 1;


for m = chamberVel_start:chamberVel_step:chamberVel_end
    m
    for t = chamberTemp_start:chamberTemp_step:chamberTemp_end
        
        for p = chamberPres_start:chamberPres_step:chamberPres_end
            
            %[throatArea(i, j, k),throatPres(i, j, k),throatTemp(i, j, k)] = nozzleThroatArea(p,m,t,chamberArea);
            %[Throat Area, Throat Pres, Throat Temp]
            
            [data(i, j, k, 1), data(i, j, k, 2), data(i, j, k, 3), data(i, j, k, 4)] = nozzleThroatArea(p, m, t, chamberArea);
            data(i, j, k, 5) = t;
           
            
            
            %Expansion Ratio
            data(i, j, k, 6) = exitArea_max / data(i, j, k, 1);
            
            if(data(i, j, k, 6) < 1)
                data(i, j, k, 8) = NaN;
                data(i, j, k, 9) = NaN;
                data(i, j, k, 1) = NaN;
                data(i, j, k, 6) = NaN;
                data(i, j, k, 7) = NaN;
                data(i, j, k, 10) = NaN;
            else
                [data(i, j, k, 7), data(i, j, k, 8), data(i, j, k, 9), data(i, j, k, 10)] = divNozzle(data(i, j, k, 6), data(i, j, k, 2), data(i, j, k, 3));
            end
            
            %Compressible Flow Calculations:
            
            
            %Reject DataPoints if the potential Exit Pressure is below flow seperation criteria
            if(constrain == 1)
                
                %Flow Seperation Constraint
                if(data(i, j, k , 8) < constrain_FlowSep * ambientPressure)
                    data(i, j, k, 8) = NaN;
                    data(i, j, k, 9) = NaN;
                    data(i, j, k, 1) = NaN;
                    data(i, j, k, 6) = NaN;
                    data(i, j, k, 7) = NaN;
                    data(i, j, k, 10) = NaN;
       
                end
                   
                %Throat Temperature Constraint
                if(data(i, j, k, 3) > constraint_ThroatTemp)
                    data(i, j, k, 8) = NaN;
                    data(i, j, k, 9) = NaN;
                    data(i, j, k, 1) = NaN;
                    data(i, j, k, 6) = NaN;
                    data(i, j, k, 7) = NaN;
                    data(i, j, k, 10) = NaN;
                end
                
            end
            
            
            k = k + 1;
        end
        j = j + 1;
        k = 1;
        
    end
    i = i + 1;
    j = 1;
end

%testPres = throatPres(1, :, :);

fprintf("Enter a Combustion Chamber Pressure Value in kPa (Steps of 10 kPa, up to " +chamberPres_end/1000 +" kPa)   ")
p = input(' ');
p = p * 1e3;
p_i = (p - chamberPres_start) / chamberPres_step;

throatPres_p = squeeze(data(1, :, :, 2));
throatTemp_p = squeeze(data(1, :, :, 3));
throatArea_p = squeeze(data(1, :, :, 1));

%For constant Pressure:
%Axis Values
temp_Vals = squeeze(data(:, :, p_i, 5));
vel_Vals = squeeze(data(:, :, p_i, 4));

%Calculated Values
area_Vals = squeeze(data(:, :, p_i, 1));
area_Vals = sqrt(area_Vals * (4/pi)) * 100 / 2.54;
expansion_Vals = squeeze(data(:, :, p_i, 6));
mach_Vals = squeeze(data(:, :, p_i, 7));
exitTemp_Vals = squeeze(data(:, :, p_i, 10));
exitVel_Vals = squeeze(data(:, :, p_i, 9));

%Exit Pressure Mesh Plot Data Sanitizing



exitPres_Vals = squeeze(data(:, :, p_i, 8));


figure(1)
cp1 = contourf(temp_Vals(:, :), vel_Vals(:, :), area_Vals(:, :), 50, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Throat Diameter to Choke vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
%legend;
c = colorbar;
c.Label.String = 'Throat Diameter (in)';

figure(2)
cp2 = contourf(temp_Vals(:, :), vel_Vals(:, :), exitVel_Vals(:, :), 500, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Exhaust Velocity vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
c = colorbar;
colormap(jet);
c.Label.String = 'Exhaust Velocity (m/s)';

figure(3)
cp3 = contourf(temp_Vals(:, :), vel_Vals(:, :), exitPres_Vals(:, :)/1e3, 500, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Exit Pressure vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
c = colorbar;
c.Label.String = 'Exit Pressure (kPa)';

figure(4)
cp3 = contourf(temp_Vals(:, :), vel_Vals(:, :), expansion_Vals(:, :), 500, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Expansion Ratio vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
c = colorbar;
c.Label.String = 'Expansion Ratio (\epsilon)';


figure(5)
cp3 = contourf(temp_Vals(:, :), vel_Vals(:, :), mach_Vals(:, :), 500, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Exit Mach vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
c = colorbar;
c.Label.String = 'Exit Mach';

figure(6)
cp3 = contourf(temp_Vals(:, :), vel_Vals(:, :), exitTemp_Vals(:, :), 500, 'edgecolor','none');
x1 = xlabel('Combustion Chamber Temperature (\circK)');
y1 = ylabel('Combustion Chamber Flow Velocity (m/s)');
t1 = title("Exit Temperature vs. Chamber Velocity and Temperature at Chamber Pressure of "+ p/1e3 +" kPa");
c = colorbar;
c.Label.String = 'Temperature Static (\circK)';

%mp1 = mesh(chamberVel_start:chamberVel_step:chamberVel_end, chamberTemp_start:chamberTemp_step:chamberTemp_end, chamberPres_start:chamberPres_step:chamberPres_end, data(i, j, k, 8));



%contourf(throatPres_p, throatTemp_p, throatArea_p)

%contourf(throatPres(1, :, :), throatTemp(1, :, :), throatArea(1, :, :));
