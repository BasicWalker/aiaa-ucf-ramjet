%% Set up 
clc; clear all

if exist('T','var')==0
    load GRAM_Model.mat
end

gamma               = 1.4;
R                   = 287;

initial_altitude    = 1100;
final_altitude      = 10000;    % Arbitrary 
n                   = ((final_altitude - initial_altitude) / 100) + 1;  % Increment altitude by 100m

starting_exitMach   = 2; 
final_exitMach      = 4; 
m                   = ((final_exitMach - starting_exitMach)/0.05) + 1;  % Increment exit mach by 0.05

%% Preallocate
altitude        = zeros(m,n);
exitTemp        = zeros(m,n);
backPres        = zeros(m,n);
exitMach        = zeros(m,n);
exitVelocity    = zeros(m,n);
areaRatio       = zeros(m,n);
exitArea        = zeros(m,n); 
throatArea      = zeros(m,n);
presRatio       = zeros(m,n);
stagPres        = zeros(m,n);
thrust          = zeros(m,n);
mdot_exit       = zeros(m,n);

areaTable       = zeros(n);

%% Begin Calculations Iterating Exit Mach Number and Altitude
for i = 1:m
    for j= 1:n
        altitude(i,j)                       = initial_altitude + 100*(j-1);
        exitTemp(i,j)                       = interp1(T.Hgtkm, T.Tmean, altitude(i,j)*0.001);
        backPres(i,j)                       = interp1(T.Hgtkm, T.PresMean, altitude(i,j)*0.001);

        exitMach(i,j)                       = 2 + 0.05*(i-1);
        exitVelocity(i,j)                   = exitMach(i,j) * sqrt(gamma * R * exitTemp(i,j));

        [presRatio(i,j), areaRatio(i,j)]    = isentropicFlow(gamma, exitMach(i,j));

        exitArea(i,j)                       = 2;
        throatArea(i,j)                     = exitArea(i,j) / areaRatio(i,j);

        stagPres(i,j)                       = backPres(i,j) / presRatio(i,j);

        thrust(i,j)                         = 758.74;
        mdot_exit(i,j)                      = thrust(i,j) / exitVelocity(i,j);
    end  
end

%% Tabulate Results in a 
for j = 1:2
    areaTable = table( exitMach(:,j), areaRatio(:,j), exitArea(:,j), throatArea(:,j) );
    file = sprintf('areaTable_%d', j);
    writetable(areaTable, file);
end


% T = table(exitMach, areaRatio, exitArea, throatArea);
% U = table(thrust, exitMach, exitVelocity, mdot_exit);

