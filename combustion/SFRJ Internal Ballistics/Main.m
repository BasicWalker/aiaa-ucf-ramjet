% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Main.m 
%
% File Description: 
% Main executive model. controls logic of the program
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial creation 
% Ethan Sherlock  02/14/21  001  Chamber pressure calculation update
% Ethan Sherlock  02/14/21  005  1DOF trajectory update
% Ethan Sherlock  03/17/21  ---  Code Clean Up
% ---------------------------------------------------------------------- %

while StopBurn == 0
    BurnTime(n) = time;                         % Simulation Time
    RegressionRate                              % Call Regression Rate Model
    GrainGeometry                               % Call Instantaneous Grain Geometry Model
    Trajectory                                  % Call Trajectory Model
    Intake                                      % Call Intake Model
    Gas                                         % Call Gas Model (And Chemistry Model)
    BoundaryLayer                               % Call Boundary Layer Model
    Thrust                                      % Call Thrust Model
    time = time + SFRJDt;                       % Step through simulation time
    n = n + 1;
end
Nozzle                                          % Call Nozzle Model
PlotData                                        % Plot Simulation Results