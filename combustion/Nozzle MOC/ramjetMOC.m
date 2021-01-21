%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name:  Ramjet Nozzle with the Method of Characteristics
% 
% Program Description: 
%
% 
% File Name: ramjetMOC.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Karam Paul      01/17/21  Initial Creation 
% --------------------------------------------------------------------- %

%% Setup
close all; clear; clc

% Load .mat file and save into table
if exist('T','var')==0
    load GRAM_Model.mat
end

%% Propulsion Propterties
chamberPres = 2.27e6;   % Chamber Pressure [Pa]
chamberTemp = 1200;     % Chamber Temperature [K]
thrust      = 1200;     % Design Thrust [N]
altitude    = 7500;     % Altitude [m]
gamma       = 1.4;      % Coefficient of Heat
gasConstant = 355;      % Gas Constant [J/kgK]

%% Atmospheric Properties
backPressure = interp1(T.Hgtkm, T.PresMean, (altitude)/1e3);

%% Calculations 
presRatio   = backPressure / chamberPressure;           % Pressure Ratio
tempRatio   = presRatio ^((gamma - 1) / gamma) ;        % Temperature Ratio   
throatTemp  = (2 * gamma * chamberTemp) / (gamma - 1);  % Throat Temp [K]
