%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name: Aerodynamic Property Calculator 
% 
% Program Description: 
% The purpose of this calculator is to solve for velocity, pressure, density  
% and temperature across our shock structure
% 
% File Name: propertyCalculator.m
% 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Karam Paul      01/124/21  Initial Creation 
% --------------------------------------------------------------------- %

%% Setup
clc ; clear ; close all 

if exist('T','var')==0
    load GRAM_Model.mat
end

%% Define Initial Conditions
mach1 = 2;
gamma = 1.4;
defl = 15;  % <deg>
initial_altitude = 1.1;    % <km>

staticDens1 = interp1(T.Hgtkm, T.DensMean, initial_altitude);  % <kg/m3>
staticPres1 = interp1(T.Hgtkm, T.PresMean, initial_altitude);   % <Pa>
staticTemp1 = interp1(T.Hgtkm, T.Tmean, initial_altitude);       % <K>

%% Calculate Stagnation Properties 
[~, tempRatio1, presRatio1, densRatio1, ~] = flowisentropic(gamma, mach1, 'mach'); 
stagDens1 = staticDens1 / densRatio1;
stagPres1 = staticPres1 / presRatio1;
stagTemp1 = staticTemp1 / tempRatio1;

%% Oblique Shock Procedure
[mach2, shockAngle] = obliqueShock(mach1, defl, gamma); 
normalMach1 = mach1 * sind(shockAngle); 

[~, tempRatio2, presRatio2, densRatio2, normalMach2, stagPresRatio2] = ...
    flownormalshock(gamma, normalMach1, 'mach'); 

%% Terminating Normal Shock Procedure 
[~, tempRatio3, presRatio3, densRatio3, normalMach3, stagPresRatio3] = ...
    flownormalshock(gamma, mach2, 'mach'); 
mach3 = normalMach3;

%% Diffuser Procedure
mach4 = 0.2; 
[~, tempRatio4, presRatio4, densRatio4, ~] = flowisentropic(gamma, mach4, 'mach'); 

%% Calculate Final Flow Properties 
staticTemp4 = tempRatio4 * stagTemp1;
staticPres4 = presRatio4 * stagPresRatio3 * stagPresRatio2 * (1/presRatio1) * staticPres1;
staticDens4 = densRatio4 * stagDens1;
stagPres4   = (1/presRatio4) * staticPres4;

%% Convert Property Units
staticTemp4_celcius = staticTemp4 - 273;
staticPres4_bar = staticPres4 * 1e-5;
stagPres4_bar = stagPres4 * 1e-5;

%% Create Table
T = table(staticTemp4_celcius, staticDens4, staticPres4_bar, stagPres4_bar);
T



