function [StagPressure] = pressureToChoke(Massflow, Area_throat, StagTemp, gamma, R)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('gamma','var')
    gamma = 1.4;
end
if ~exist('R','var')
    R = 287;
end
StagPressure = Massflow*sqrt(StagTemp)/Area_throat*sqrt(R/gamma)*((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));
end

