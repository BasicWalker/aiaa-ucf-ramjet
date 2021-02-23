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
stagPres(2) = stagPres(1) * stagPresRatio(2);
stagDens(2) = stagDens(1);
stagTemp(2) = stagTemp(1);
staticDens(2) = staticDens(1) * densRatio(2);
staticPres(2) = staticPres(1) * presRatio(2);
staticTemp(2) = staticTemp(1) * tempRatio(2);
[~, ~, ~, ~, aThroat_a2star] = flowisentropic(gamma, mach(2), 'mach'); 

% Terminating Normal Shock Procedure 
[~, tempRatio(3), presRatio(3), densRatio(3), mach(3), stagPresRatio(3)] = ...
    flownormalshock(gamma, mach(2), 'mach'); 
stagPres(3) = stagPres(2) * stagPresRatio(3);
stagDens(3) = stagDens(2);
stagTemp(3) = stagTemp(2);
staticDens(3) = staticDens(2) * densRatio(3);
staticPres(3) = staticPres(2) * presRatio(3);
staticTemp(3) = staticTemp(2) * tempRatio(3);

radiusThroat(1) = radiusCombustor;  % initial guess
i = 1;

while true
    areaThroat(i) =  pi*radiusThroat(i)^2;
    area2star(i) = areaThroat(i) / aThroat_a2star;
    area3star(i) = area2star(i) / stagPresRatio(3);
  
   [mach4(i), tempRatio4(i), presRatio4(i), densRatio4(i), ~] = ...
       flowisentropic(gamma, (areaCombustor/area3star(i)), 'sub'); 
    staticDens4(i) = stagDens(3)*densRatio4(i);
    staticPres4(i) = stagPres(3)*presRatio4(i);
    staticTemp4(i) = stagTemp(3)*tempRatio4(i);
    
    change(i) = (staticPres4(i)/pressureCombustor - 1);
    error(i) = pressureCombustor - staticPres4(i);
    
    if abs(change(i)) < eps('single')
        break
    end
    radiusThroat(i+1) = radiusThroat(i) + change(i)/10;
    i = i+1;
end

plot(areaThroat,abs(error))

% [x,ind] = min(abs(error));
% disp(x)
% disp(mach1(ind))

