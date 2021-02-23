%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name: Normal Shock Placement 

% Program Description: 
%
% 
% File Name: intake_solver.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Karam, Samer    02/22/21  Initial Creation 
% --------------------------------------------------------------------- %
clc; 

if exist('T','var')==0
    load GRAM_Model.mat
end






% define variables
gamma = 1.4;
R = 287;  %<J/(kg*K)>
mach(1) = 2;                  % free stream mach
deflectionAngle = 5;
altitude = 5000e-3;     % altitude at mach 2 (erich's code)    
radiusCombustor = 0.5;  % in
areaCombustor =  pi*radiusCombustor^2;
pressureCombustor = 400e3;  % kPa
machCombustorLim = 0.2;

%free stream values
staticDens(1) = interp1(T.Hgtkm, T.DensMean, altitude);       % <kg/m3>
staticPres(1) = interp1(T.Hgtkm, T.PresMean, altitude);       % <Pa>
staticTemp(1) = interp1(T.Hgtkm, T.Tmean, altitude);          % <K>
velocity(1) = mach(1)*sqrt(gamma*R*staticTemp(1));
[~, tempRatio(1), presRatio(1), densRatio(1), ~] = flowisentropic(gamma, mach(1), 'mach'); 
stagDens(1) = staticDens(1) / densRatio(1);
stagPres(1) = staticPres(1) / presRatio(1);
stagTemp(1) = staticTemp(1) / tempRatio(1);

% Oblique Shock Procedure
[mach(2), shockAngle] = obliqueShock(mach(1), deflectionAngle, gamma); 
normalCompMach1 = mach(1) * sind(shockAngle); 

[~, tempRatio(2), presRatio(2), densRatio(2), ~, stagPresRatio(2)] = ...
    flownormalshock(gamma, normalCompMach1, 'mach'); 
staticDens(2) = staticDens(1) * densRatio(2);
staticPres(2) = staticPres(1) * presRatio(2);
staticTemp(2) = staticTemp(1) * tempRatio(2);
[~, ~, ~, ~, aThroat_a2star] = flowisentropic(gamma, mach(2), 'mach'); 

while true
    radiusThroat = radiusCombustor
    areaThroat =  pi*radiusThroat^2;
    area2star = areaThroat / aThroat_a2star;

    % Terminating Normal Shock Procedure 
    [~, tempRatio(3), presRatio(3), densRatio(3), mach(3), stagPresRatio(3)] = ...
        flownormalshock(gamma, mach(2), 'mach'); 
    stagDens(3) = staticDens(2) * densRatio(3);
    stagPres(3) = staticPres(2) * presRatio(3);
    stagTemp(3) = staticTemp(2) * tempRatio(3);
    
    area3star = area2star / stagPresRatio(3);
  
   [mach(4), tempRatio(4), presRatio(4), densRatio(4), ~] = flowisentropic(gamma, (areaCombustor/area3star), 'sub'); 
    staticDens(4) = stagDens(3)*densRatio(4);
    staticPres(4) = stagPres(3)*presRatio(4);
    staticTemp(4) = stagTemp(3)*tempRatio(4);
    
    change = (staticPres(4)/pressureCombustor - 1)
    radiusThroat = radiusThroat + change 
    if change < eps('single')
        break
    end
end



