% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: RegressionRate.m 
%
% File Description: 
% Regression rate model, calculates fuel regression rate
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation  
% ---------------------------------------------------------------------- %

% Simplified linear regression rate assumption
 fuel.BurnRate     = 0.001;         % Burn Rate (m/s)
 fuel.Regression = fuel.BurnRate*SFRJDt; % Regression Per Step (m)