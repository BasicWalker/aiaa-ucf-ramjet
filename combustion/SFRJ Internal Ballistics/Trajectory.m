% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Trajectory.m 
% 
% File Description: 
% Trajectory model for the SFRJ
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %

% Assume constant height till this model can be developed
Height(n) = Height(1);

% Interpolate to get back pressure
BackPres(n) = interp1(AltitudeTbl(:,2),AltitudeTbl(:,5),Height(n),'spline');