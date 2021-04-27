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


fuel.DiaInner(n+1) = fuel.DiaInner(n) + 2 * fuel.Regression;     % Fuel Grain Inner Diameter (m)

fuel.PortArea(n) = pi*(fuel.DiaInner(n)^2)*(1/4);          % Fuel Port Area (m^2)
fuel.CsxArea(n) = pi*(fuel.DiaOuter^2)*(1/4) - fuel.PortArea(n); % Fuel Grain Crossectional Area (m^2)
fuel.Volume(n) = fuel.CsxArea(n) * fuel.Length;                % Fuel Grain Volume (m^3)
fuel.SurfArea(n) = fuel.DiaInner(n)* pi * fuel.Length;            % Fuel Grain Surface Area (m^2)


fuel.StepHeight(n) = (fuel.DiaInner(n) - combustion.InletDia)/2;             % Rearward Step Height

% Fuel Mass Properties
fuel.MassGen(n) = fuel.Regression*fuel.Density*fuel.SurfArea(n);         % Fuel mass generated every time step (kg)
fuel.MassFlow(n) = fuel.MassGen(n)/SFRJDt;                   % Fuel mass flow rate (kg/s)
fuel.Mass(n) = fuel.Density*fuel.Volume(n);                   % Grain fuel mass, instantaneous (kg)

% Stop Simulation Flag
if fuel.DiaInner(n) > fuel.DiaOuter
    StopBurn = true;
    Burnout = true;
    index = n;
end