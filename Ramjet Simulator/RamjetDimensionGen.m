clear;close all;clc
constants.In2Mtr = 39.3701;  % Inch to meter converstion 




% combustion -----
combustion.InletArea = 0.0064693;  %pi*combustion.InletDia^2*(1/4);  % (m)
combustion.InletDia = sqrt(4*combustion.InletArea/pi);  %  <m>   
combustion.InletDiaINCH = sqrt(4*combustion.InletArea/pi)* constants.In2Mtr;  % <in>    %1.4 / constants.In2Mtr;  % (m)
combustion.ChamberArea = combustion.InletArea*1.5;  % 50% larger than inlet area     %pi*combustion.ChamberRadius^2/4; %pi*combustion.ChamberRadius^2;  % Area of combustion chamber (m^2)
combustion.ChamberDia = sqrt(4*combustion.ChamberArea/pi); %<m>            %2.75 /constants.In2Mtr / 2;  % Radius of the combustion chamber (m)
combustion.ChamberDiaINCH = sqrt(4*combustion.ChamberArea/pi)* constants.In2Mtr; %<in>

% intake -----
intake.Area_enter = 0.0007;  % Area of throat (m^2) - Drives mass flow rate through intake
intake.CowlDia = combustion.ChamberDia;  % <m>
intake.CowlDiaINCH = intake.CowlDia*constants.In2Mtr;  % <in>
intake.EnterDia = sqrt(4/pi*(pi*intake.CowlDia^2/4 - intake.Area_enter));  % <m>
intake.EnterDiaINCH = intake.EnterDia*constants.In2Mtr;  % <in>
intake.DeflAngle = 15;   % Deflection angle (deg)

% fuel -----
fuel.DiaOuter =  combustion.ChamberDia;%2.75 /constants.In2Mtr;  % Grain OD (m)
fuel.DiaOuterINCH =  combustion.ChamberDia*constants.In2Mtr;
fuel.StepHeight(1) = 1/8 /constants.In2Mtr;
fuel.DiaInner(1) = 2*fuel.StepHeight(1) + combustion.InletDia;    %1.50 /constants.In2Mtr;   
fuel.DiaInnerINCH(1) = (2*fuel.StepHeight(1) + combustion.InletDia)* constants.In2Mtr;    %1.50 /constants.In2Mtr;  % Grain ID (m)
fuel.Length = 15.00 /constants.In2Mtr;  % Grain Length (m)
fuel.Density = 1020;  % Grain Density (kg/m^3)

% nozzle -----
nozzle.DiaThroat = 1.8 /constants.In2Mtr;  % Throat Diameter, assuming exit area is 1.6 in diameter (from HPR), 0.985
nozzle.Area_throat = (pi/4)*(nozzle.DiaThroat)^2;  % Throat area (m^2)
nozzle.DiaExit = 2.75/constants.In2Mtr;
nozzle.Area_exit = (pi/4)*(nozzle.DiaExit)^2;  % nozzle exit area (m^2)

% Ramjet Dimensions
vehicle.DragCoeff = 0.23;  % Drag coefficient (0.35)
vehicle.FrontSurfArea = combustion.ChamberArea; % 0.008119;  % Frontal surface area (m^2)
vehicle.DryMass = 6.80389;  % Mass of ramjet without fuelgrain (kg)


% save workspace variables to mat file
clear constants.In2Mtr
save('RamjetDimensions.mat')


