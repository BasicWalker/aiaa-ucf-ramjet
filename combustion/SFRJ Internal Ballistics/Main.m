% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Main.m 
%
% File Description: 
% Main executive model. controls logic of the program
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %

while StopBurn == 0
    BurnTime(n) = time;                         % Simulation Time
    RegressionRate                              % Call Regression Rate Model
    GrainGeometry                               % Call Instantaneous Grain Geometry Model
    
    % Calculate Fuel Mass Properties
    MFuelGen(n) = RgrsPerStp*FuelRho*FuelSA(n); % Fuel mass generated
    FuelMass(n) = FuelRho*FuelVol(n);           % Fuel mass 
    
    Gas                                         % Call Gas Model
    Chemistry                                   % Call Chemistry Model
    BoundaryLayer                               % Call Boundary Layer Model
    Nozzle                                      % Call Nozzle Model

    % Fuel Mass Prop continued & O/F Ratio calculation
    OFRatio(n) = MOxdzrGen(n)/MFuelGen(n);      % O/F Ratio
    MassGen(n) = MairGen(n) + MFuelGen(n);      % Total mass generated (kg)
    MassFlow(n) = MassGen(n)/SFRJDt;            % Total mass flow (kg/s)
    AFRatio(n) = MairGen(n)/MFuelGen(n);        % Fuel to air ratio 
    
    % Print O/F Ratio warning
    if OFRatio(n) > 10
        fprintf('WARNING: O/F Ratio too high \n')
    end
    
    if OFRatio(n) < 3
        fprintf('WARNING: O/F Ratio too low \n')
    end
     
    fprintf('Running... \n')                    % Running Simulator indicator
    
    Thrust                                      % Call Thrust Model
    Trajectory                                  % Call Trajectory Model
    
    % Step through simulation time
    time = time + SFRJDt;
    n = n + 1;
end
PlotData