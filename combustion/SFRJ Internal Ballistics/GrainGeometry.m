% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: GrainGeometry.m 
%
% File Description: 
% Geometry model, calculates instantaneous fuel grain geometry
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %

if n > 1
    GrainID(n) = GrainID(n-1) + 2 * RgrsPerStp;     % Fuel Grain Inner Diameter (m)
    PortArea(n) = pi*(GrainID(n)^2)*(1/4);          % Fuel Port Area (m^2)
    FuelCS(n) = pi*(GrainOD^2)*(1/4) - PortArea(n); % Fuel Grain Crossectional Area (m^2)
    FuelVol(n) = FuelCS(n) * GrainL;                % Fuel Grain Volume (m^3)
    FuelSA(n) = GrainID(n)* pi * GrainL;            % Fuel Grain Surface Area (m^2)
end

StepHeight(n) = (GrainID(n) - InltD)/2;             % Rearward Step Height

% Fuel Mass Properties
MFuelGen(n) = RgrsPerStp*FuelRho*FuelSA(n);         % Fuel mass generated every time step (kg)
MdotFuel(n) = MFuelGen(n)/SFRJDt;                   % Fuel mass flow rate (kg/s)
FuelMass(n) = FuelRho*FuelVol(n);                   % Grain fuel mass, instantaneous (kg)


% Stop Simulation Flag
if GrainID(n) > GrainOD
    StopBurn = true;
    Burnout = true;
    %fprintf('Motor Burnout Detected\n')
    index = n;
end