%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            AIAA UCF Ramjet Intake Geometry Design Script                %
%                                                                        %
%                              Authored by                                %
%        Samer Armaly, Karam Paul, Jared Durlak, Matthew Aubertin         %
%                           October 28, 2020                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;

%% define variables

delta = [1 2 3 4 5];  % deflection angle array
n = size(delta,2);  % number of oblique shocks
gamma = 1.4;  % specific heat ratio for flow

% pre-allocate arrays
mach = zeros(n+2,1);  % array for mach numbers
pressure = zeros(n+2,1);  % array for static pressure
stagPressure = zeros(n+2,1);
temperature = zeros(n+2,1);  % array for temperature
density = zeros(n+2,1);  % array for density
ramp = zeros(n+1,2);  % array for ramp coordinates
theta = zeros(n,1);  % array for each shock angle formed
lambda = zeros(n,1);
phi = zeros(n,1);
mach_normal = zeros(n,1);  % array for the normal components crossing the oblique wave

%input known terms 
mach(1) = 3;  % freestream mach
pressure(1) = 26.5e3;  % freestream static pressure (Pa)
[~, tempRatio, presRatio, ~, ~] = flowisentropic(gamma, mach(1),'mach');
stagPressure(1) = pressure(1) * presRatio; 

cowl_height = 1.25;  % inner radius for cowl <inches>

%% oblique shocks
% solving for weak oblique shocks (imaginary components signify detached shock)
for i = 1:n+1
    try
        if i < n+1
            [~, theta(i)] = obliqueShock(mach(i), delta(i));
            mach_normal(i) = mach(i)*sind(theta(i));
        else
            mach_normal(i) = mach(i);
        end

        % solve for normal shock relations
        [~, t_r, p_r, d_r, M, po_r, ~] = flownormalshock(gamma, mach_normal(i), 'mach');

        if i < n+1
            mach(i+1) = M / sind(theta(i) - delta(i));
        else
            mach(i+1) = M;
        end

        pressure(i+1) = pressure(i)*p_r;  % multiplies upstream pressure by normal shock ratio to get downstream value
        temperature(i+1) = temperature(i)*t_r;  % multiplies upstream temperature by normal shock ratio to get downstream value
        density(i+1) = density(i)*d_r;  % multiplies upstream density by normal shock ratio to get downstream value
        stagPressure(i+1) = po_r * stagPressure(i);
    catch warning('detached shock due to imaginary solution');
        break
    end
end

%% geometry propogation
% find spike length by intersecting initial shock with cowl height
cowl_x = cowl_height / tand(theta(1));


% define delta and theta wrt the horizontal
phi(1) = delta(1);  % deflection angle wrt to the horizontal
lambda(1) = theta(1);  
for i = 2:n
    phi(i) = delta(i) + phi(i-1);
    lambda(i) = theta(i) + phi(i-1);
end

% loop to find additional ramp geometry by intersecting n-shock with ramp
if n >1  % if there are more than one oblique shock solve for the rest
    for i = 2:n
        % solve for ramp x locations
        ramp(i,1) = (tand(lambda(i))*cowl_x - cowl_height - tand(phi(i-1))*ramp(i-1,1) + ramp(i-1,2))/(tand(lambda(i))-tand(phi(i-1)));
        % solve for ramp y locations
%         ramp(i,2) = tan(lambda(i))*(ramp(i,1)-cowl_x)+cowl_height;
        ramp(i,2) = tand(phi(i-1))*(ramp(i,1)-ramp(i-1,1)) + ramp(i-1,2);
    end
end

% solve for ramp termination points
ramp(n+1,1) = (cowl_x/tand(90+phi(n))-cowl_height-tand(phi(n))*ramp(n,1)+ramp(n,2))/...
    (1/tand(90+phi(n))-tand(phi(n)));

ramp(n+1,2) = 1/tand(90+phi(n))*(ramp(n+1,1)-cowl_x) + cowl_height;
% ramp(n+1,2) = tand(90 + phi(n))*(ramp(n+1,1)-ramp(n,1));
centerbody_height = ramp(n+1,2);
isolator_height = cowl_height - centerbody_height;
lip_length = isolator_height * tand(phi(n));
spike_length = cowl_x + lip_length;


% 
%% plotting
f1 = figure('Name','ramp shock structure');
hold on
% plot ramp
plot(ramp(:,1),ramp(:,2), 'k', 'LineWidth',2);
plot([0,spike_length],[0,0],':k', 'LineWidth',2);
plot([spike_length,spike_length],[0,ramp(end,2)],':k', 'LineWidth',2);
% 
% plot shocks
for i= 1:n+1
    % check if shocks are attached
    if isreal(mach(i+1))
        plot([ramp(i,1),cowl_x],[ramp(i,2),cowl_height], '--c');
    % show detached shock origination
    else
        plot([ramp(i,1),spike_length],[ramp(i,2),cowl_inner_dia], 'r*')
        t = text(ramp(i,1),ramp(i,2)+(cowl_height - ramp(i,2))/5,'detached');
        t.Color = 'red';
    end
end
% hold off
% f2 = figure('Name','Pressure throughout intake');
% hold on
% plot([0; 1; 2; 3; 4], pressure)


