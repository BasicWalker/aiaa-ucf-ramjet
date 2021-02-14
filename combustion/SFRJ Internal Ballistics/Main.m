% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Main.m 
%
% File Description: 
% Main executive model. controls logic of the program
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  001  Chamber Pressure Calculation Update
% ---------------------------------------------------------------------- %

while StopBurn == 0
    BurnTime(n) = time;                         % Simulation Time
    Trajectory                                  % Call Trajectory Model
    RegressionRate                              % Call Regression Rate Model
    GrainGeometry                               % Call Instantaneous Grain Geometry Model
    
    % Calculate Fuel Mass Properties
    MFuelGen(n) = RgrsPerStp*FuelRho*FuelSA(n); % Fuel mass generated
    MdotFuel(n) = MFuelGen(n)/SFRJDt;           % Fuel mass flow rate
    FuelMass(n) = FuelRho*FuelVol(n);           % Fuel mass 
    
    [AFRst, gamma_nzl, T_AFT] = Chemistry(Phi); % Call Chemistry Model
    T_stag(n) = T_AFT;                          % Stagnation temp in nozzle
    gamma_nzlT(n) = gamma_nzl;                  % Gamma at the nozzle throat
    Gas                                         % Call Gas Model

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
   
    % Step through simulation time
    time = time + SFRJDt;
    n = n + 1;
end
PlotData