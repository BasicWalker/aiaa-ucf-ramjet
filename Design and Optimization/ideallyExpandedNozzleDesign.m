clear; clc; close all;

addpath(genpath(pwd))
addpath('..\Common Resources')

% Import data
load RamjetDimensions.mat  % load in the ramjet design
load GRAM_Model.mat  % GRAM atmospheric model
load Constants.mat  % load in constants and conversions

% Run Simulation to Get All Current Values
Main
clc; close all
fprintf("Simulation Complete. Results are saved in the Workspace\n");

% Pre-allocate Resources
error = zeros(n,1);

% Flight Mach Number Where Nozzle Should Perform Ideally
interestMach = 3;

% Find Index at which Flight Mach Number Occurs
for i = 1:n
    error(i) = intake.mach(1,i) - interestMach;
end
[minError,index] = min(abs(error(1:end)));

% Store Necessary Properties at that Flight Mach Number
gamma = nozzle.gamma(index);
ambientPres = intake.staticPres(1,index);
nozzleStagPres = combustion.stagPres(2,index);
throatStaticPres = nozzle.staticPres(1,index);
exitPresRatio = ambientPres / nozzleStagPres;

% Find Mach Number and Expansion Ratio Required for Ideal Expansion
[exitMach,~,~,~,expansionRatio] = flowisentropic(gamma,exitPresRatio,'pres');

% Calculate Exit Dimensions
nozzleExitArea_m = nozzle.Area_throat * expansionRatio; % <m2>
nozzleExitArea_in = nozzleExitArea_m * (constants.In2Mtr)^2; % <in2>
nozzleExitDiameter_m = sqrt( (4/pi)*nozzleExitArea_m); % <m>
nozzleExitDiameter_in = sqrt( (4/pi)*nozzleExitArea_in) % <in>

