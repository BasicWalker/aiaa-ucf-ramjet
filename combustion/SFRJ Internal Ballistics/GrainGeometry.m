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

% Initialize Grain Geometry Model
PortArea(1) = pi*(GrainID(1)^2)*(1/4);              % Fuel Port Area
FuelCS(1) = pi*(GrainOD^2)*(1/4) - PortArea(1);     % Fuel Grain Crossectional Area
FuelVol(1) = FuelCS(1) * GrainL;                    % Fuel Grain Volume
FuelSA(1) = GrainID(1)* pi * GrainL;                % Fuel Grain Surface Area

% After Initialization 
if n > 1
    GrainID(n) = GrainID(n-1) + 2 * RgrsPerStp;     % Fuel Grain Inner Diameter
    PortArea(n) = pi*(GrainID(n)^2)*(1/4);          % Fuel Port Area
    FuelCS(n) = pi*(GrainOD^2)*(1/4) - PortArea(n); % Fuel Grain Crossectional Area
    FuelVol(n) = FuelCS(n) * GrainL;                % Fuel Grain Volume
    FuelSA(n) = GrainID(n)* pi * GrainL;            % Fuel Grain Surface Area
end

StepHeight(n) = (GrainID(n) - InltD)/2;             % Rearward Step Height

if GrainID(n) > GrainOD
    StopBurn = true;
end


