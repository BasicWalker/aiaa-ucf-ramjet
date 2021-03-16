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
% ---------------------------------------------------------------------- %

while StopBurn == 0
    BurnTime(n) = time;                         % Simulation Time
    
    RegressionRate                              % Call Regression Rate Model
    GrainGeometry                               % Call Instantaneous Grain Geometry Model
    Gas                                         % Call Gas Model
    BoundaryLayer                               % Call Boundary Layer Model
    Nozzle                                      % Call Nozzle Model

    % Fuel Mass Prop continued & O/F Ratio calculation
    OFRatio(n) = MOxdzrGen(n)/MFuelGen(n);      % O/F Ratio
    MassGen(n) = MairGen(n) + MFuelGen(n);      % Total mass generated (kg)
    MassFlow(n) = MassGen(n)/SFRJDt;            % Total mass flow (kg/s)
    AFRatio(n) = MairGen(n)/MFuelGen(n);        % Fuel to air ratio 
     
    Thrust                                      % Call Thrust Model
    Trajectory                                  % Call Trajectory Model
   
    time = time + SFRJDt;                       % Step through simulation time
    n = n + 1;
end
PlotData