%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            AIAA UCF Ramjet Intake Geometry Design Script                %
%                                                                        %
%                              Authored by                                %
%        Samer Armaly, Karam Paul, Jared Durlak, Matthew Aubertin         %
%                           October 28, 2020                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

%% define variables

defl = [5 8 12];  % deflection angle array
n = size(defl,2);  % number of oblique shocks

% pre-allocate arrays
mach = zeros(n+2,1);  % array for mach numbers
pressure = zeros(n+2,1);  % array for static pressure
temperature = zeros(n+2,1);  % array for temperature
density = zeros(n+2,1);  % array for density
ramp = zeros(n+1,2);  % array for ramp coordinates
theta = zeros(n,1);  % array for each shock angle formed
mach_normal = zeros(n,1);  % array for the normal components crossing the oblique wave


%input known terms 
mach(1) = 1.5;  % freestream mach
pressure(1) = 26.5e3;  % freestream static pressure (Pa)

cowl_inner_dia = 2.5;  % inner diameter for cowl (inches)
gamma = 1.4;  % specific heat ratio for flow

%% oblique shocks
% solving for weak oblique shocks (imaginary components signify detached shock)
for i = 1:n
    if i == 1
        [mach(i+1), theta(i)] = obliqueShock(mach(i), defl(i));
    else
        [mach(i+1), theta(i)] = obliqueShock(mach(i), (defl(i)-defl(i-1)) );
    end
end

% %% normal shock
% %solve for normal shock
% [mach(end)] = normalShock(mach(n-1));

%% flow properties
% find normal components of flow as it approaches oblique
for i=1:n
    mach_normal(i) = mach(i)*sind(theta(i));
end
% since last shock is normal shock there are no components
mach_normal(n+1) = mach(n+1);


[~, t_r, p_r, d_r, M, ~, ~] = flownormalshock(gamma, mach_normal, 'mach');

mach(end) = M(end);

for i=2:n+2
    pressure(i) = pressure(i-1) * p_r(i-1);
end



%% geometry propogation
% find spike length by intersecting initial shock with cowl height
spike_length = cowl_inner_dia/tand(theta(1));
ramp(end,1) = spike_length;  % ramp termination x value

% loop to find additional ramp geometry by intersecting n-shock with ramp
if n >1  % if there are more than one oblique shock solve for the rest
    for i = 2:n
        % solve for ramp x locations
        ramp(i,1) = ( cowl_inner_dia - tand(theta(i))*spike_length )/...
            ( tand(defl(i-1)) - tand(theta(i)) );
        % solve for ramp x locations
        ramp(i,2) = tand(defl(i-1))*ramp(i,1);
    end
end

% ramp termination y value
ramp(end,2) = tand(defl(n))*spike_length;


%% plotting
f1 = figure('Name','ramp shock structure');
hold on
% plot ramp
plot(ramp(:,1),ramp(:,2), 'k', 'LineWidth',2);
plot([0,spike_length],[0,0],':k', 'LineWidth',2);
plot([spike_length,spike_length],[0,ramp(end,2)],':k', 'LineWidth',2);

% plot shocks
for i= 1:n
    % check if shocks are attached
    if isreal(mach(i+1))
        plot([ramp(i,1),spike_length],[ramp(i,2),cowl_inner_dia], '--c');
    % show detached shock origination
    else
        plot([ramp(i,1),spike_length],[ramp(i,2),cowl_inner_dia], 'r*')
        t = text(ramp(i,1),ramp(i,2)+(cowl_inner_dia - ramp(i,2))/5,'detached');
        t.Color = 'red';
    end
end
hold off
f2 = figure('Name','Pressure throughout intake');
hold on
plot([0; 1; 2; 3; 4], pressure)


