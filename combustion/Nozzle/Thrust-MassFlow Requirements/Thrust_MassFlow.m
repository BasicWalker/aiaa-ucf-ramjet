%% Set up 
clc; clear
tic

if exist('T','var')==0
    load GRAM_Model.mat
end

gamma               = 1.4;
R                   = 287;

alt_step            = 100;
initial_altitude    = 1100;
final_altitude      = 10000;    % Arbitrary 
n                   = ((final_altitude - initial_altitude) / alt_step) + 1;  % Increment altitude by 100m

mach_step           = 0.1;
starting_exitMach   = 2; 
final_exitMach      = 3; 
m                   = ((final_exitMach - starting_exitMach)/mach_step) + 1;  % Increment exit mach by 0.05

%% Preallocate
altitude        = zeros(m,n);
exitTemp        = zeros(m,n);
backPres        = zeros(m,n);
exitMach        = zeros(m,n);
exitVelocity    = zeros(m,n);
areaRatio       = zeros(m,n);
exitArea        = zeros(m,n); 
exitDiameter    = zeros(m,n); 
throatArea      = zeros(m,n);
throatDiameter  = zeros(m,n);
presRatio       = zeros(m,n);
stagPres        = zeros(m,n);
thrust          = zeros(m,n);
mdot_exit       = zeros(m,n);

mdotData        = zeros(m,4);

%% Begin Calculations Iterating Exit Mach Number and Altitude
for i = 1:m
    for j= 1:n
        altitude(i,j)                       = initial_altitude + alt_step*(j-1);
        exitTemp(i,j)                       = interp1(T.Hgtkm, T.Tmean, altitude(i,j)*0.001);
        backPres(i,j)                       = interp1(T.Hgtkm, T.PresMean, altitude(i,j)*0.001);

        exitMach(i,j)                       = starting_exitMach + mach_step*(i-1);
        exitVelocity(i,j)                   = exitMach(i,j) * sqrt(gamma * R * exitTemp(i,j));

        [presRatio(i,j), areaRatio(i,j)]    = isentropicFlow(gamma, exitMach(i,j));

        exitArea(i,j)                       = 2;
        exitDiameter(i,j)                   = sqrt( (4/pi) * exitArea(i,j) );
        
        throatArea(i,j)                     = exitArea(i,j) / areaRatio(i,j);
        throatDiameter(i,j)                 = sqrt( (4/pi) * throatArea(i,j) );


        stagPres(i,j)                       = backPres(i,j) / presRatio(i,j);

        thrust(i,j)                         = 200;
        mdot_exit(i,j)                      = thrust(i,j) / exitVelocity(i,j);
    end  
end

%% Tabulate Results
warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;

% Tabulate Exit and Throat Areas with Changing Exit Mach Number
areaVarNames        = {'Exit Mach', 'Area Ratio', 'Exit Area', 'Exit Diameter', ...
                        'Throat Area', 'Throat Diameter'};
areaTable           = table( exitMach(:,1) , areaRatio(:,1), exitArea(:,1), exitDiameter(:,1), ...
                        throatArea(:,1), throatDiameter(:,1), 'VariableNames', areaVarNames );
writetable(areaTable, 'mach_vs_area.xlsx' );

% Tabulate varying mdot,exit requirements based on varying exit mach and altitudes
mdotVarNames_altSheets        = {'Thrust', 'Exit Mach', 'Exit Velocity', 'Exit Mass Flow'};
for j = 1:n
    sheet       = j;
    mdotTable_altSheets   = table ( thrust(:,j), exitMach(:,j), exitVelocity(:,j), mdot_exit(:,j), ...
        'VariableNames', mdotVarNames_altSheets );     
    writetable( mdotTable_altSheets, 'exitMassFlow_altSheets.xlsx', 'sheet', sheet );
end

mdotVarNames_machSheets        = {'Thrust', 'Altitude', 'Exit Velocity', 'Exit Mass Flow'};
for i = 1:m
    sheet       = i;
    mdotTable_machSheets   = table ( thrust(i,:)', altitude(i,:)', exitVelocity(i,:)', mdot_exit(i,:)', ...
        'VariableNames', mdotVarNames_machSheets );     
    writetable( mdotTable_machSheets, 'exitMassFlow_machSheets.xlsx', 'sheet', sheet );
end
toc