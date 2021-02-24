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
clc; clear; close all

% Station definitions used for indices
% station 1: free stream (ambient)
% station 2: downstream of oblique shock, upstream of normal shock (external ramp)
% station 3: downstream of normal shock, upstream of subsonic iscentropic expansion (throat)
% station 4: downstream of subsonic iscentropic expansion combustor inlet,(combustion chamber opening)

% load in atmosphere model based on altitude in km
if exist('T','var')==0
    load GRAM_Model.mat
end

% define variables
gamma = 1.4;  % specific heat ratio air
R = 287;  % gas constant air <J/(kg*K)>
mach1 = 2.1;  % free stream mach
deflectionAngle = 17;  % single intake ramp deflection angle <deg>
altitude = 5;     % altitude at mach 2 (erich's code) <km>    
radiusCowl = 1.25*0.0254;  % radius of cowl <m>
radiusCombustor = 0.5*0.0254;  % radius of combustion chamber inlet opening <m>
areaCombustor =  pi*radiusCombustor^2;  % area of combustor inlet <m^2>
pressureCombustor = 400e3;  % combustion chamber pressure <kPa>
% machCombustorLim = 0.2;  

stepLim = 1000;  % convergence iteration step limit
tolerance = 0.0001;  % tolerance for residual

% Station 1 properties (free-stream)-----
staticDens1 = interp1(T.Hgtkm, T.DensMean, altitude);  % <kg/m3>
staticPres1 = interp1(T.Hgtkm, T.PresMean, altitude);  % <Pa>
staticTemp1 = interp1(T.Hgtkm, T.Tmean, altitude);  % <K>
[~, tempRatio1, presRatio1, densRatio1, ~] = flowisentropic(gamma, mach1, 'mach');  % ratios are static over stagnation
stagDens1 = staticDens1 / densRatio1;  % <kg/m3>
stagPres1 = staticPres1 / presRatio1;  % <Pa>
stagTemp1 = staticTemp1 / tempRatio1;  % <K>
velocity1 = mach1*sqrt(gamma*R*staticTemp1);  % <m/s>
if stagPres1 < pressureCombustor
    error('intake:ambientStagLow',...
        'Error. \nStagnation pressure too small; cannot reach combustion chamber pressure. \nTry increasing speed or ambient static pressure.');
end

% Station 2 properties (oblique shock)-----
[mach2, shockAngle] = obliqueShock(mach1, deflectionAngle, gamma);  % uniform mach number of external ramp region
normalCompMach1 = mach1 * sind(shockAngle);  % flow component crossing normal to oblique shock
[~, tempRatio2, presRatio2, densRatio2, ~, stagPresRatio2] = ...
    flownormalshock(gamma, normalCompMach1, 'mach');  % ratios are downstream over upstream
stagPres2 = stagPres1 * stagPresRatio2;  % <Pa>
stagDens2 = stagDens1;  % does not change over normal shock <kg/m3>
stagTemp2 = stagTemp1;  % does not change over normal shock <K>
staticDens2 = staticDens1 * densRatio2;  % <kg/m3>
staticPres2 = staticPres1 * presRatio2;  % <Pa>
staticTemp2 = staticTemp1 * tempRatio2;  % <K>
velocity2 = mach2*sqrt(gamma*R*staticTemp2);  % <m/s>
[~, ~, ~, ~, aThroat_a2star] = flowisentropic(gamma, mach2, 'mach'); 
if stagPres2 < pressureCombustor
    error('intake:obliqueStagLoss',...
        'Error. \nOblique shock stagnation pressure loss too great; cannot reach combustion chamber pressure. \nTry increasing deflection angle or speed.');
end

% Station 3 properties (normal shock)-----
[~, tempRatio3, presRatio3, densRatio3, mach3, stagPresRatio3] = ...
    flownormalshock(gamma, mach2, 'mach');  % ratios are downstream over upstream
stagPres3 = stagPres2 * stagPresRatio3;  % <Pa>
stagDens3 = stagDens2;  % does not change over normal shock <kg/m3>
stagTemp3 = stagTemp2;  % does not change over normal shock <K>
staticDens3 = staticDens2 * densRatio3;  % <kg/m3>
staticPres3 = staticPres2 * presRatio3;  % <Pa>
staticTemp3 = staticTemp2 * tempRatio3;  % <K>
velocity3 = mach3*sqrt(gamma*R*staticTemp3);  % <m/s>
if stagPres3 < pressureCombustor
    error('intake:normalStagLoss',...
        'Error. \nNormal shock stagnation pressure loss too great; cannot reach combustion chamber pressure. \nTry reducing the strength of the normal shock.');
else
    % Station 4 properties (subsonic iscentropic expansion)-----
    stagPres4 = stagPres3;  % does not change; iscentropic <Pa>
    stagDens4 = stagDens3;  % does not change; iscentropic <kg/m3>
    stagTemp4 = stagTemp3;  % does not change; iscentropic <K>
    
%     radiusThroat(1) = radiusCowl/2;  % first guess
%     radiusThroat(2) = radiusThroat(1) - 0.01;  % second guess


    % pre-allocate arrays
    areaThroat = zeros(stepLim,1);
    area2star = zeros(stepLim,1);
    area3star = zeros(stepLim,1);
    aCombustor_a3star = zeros(stepLim,1);
    mach4 = zeros(stepLim,1);
    tempRatio4 = zeros(stepLim,1);
    presRatio4 = zeros(stepLim,1);
    densRatio4 = zeros(stepLim,1);
    staticDens4 = zeros(stepLim,1);
    staticPres4 = zeros(stepLim,1);
    staticTemp4 = zeros(stepLim,1);
    residual = zeros(stepLim,1);
    convergeFlag = 0;
    ind = 0;
    areaThroat(1) = areaCombustor/5;
    areaThroat(2) = areaThroat(1)-0.00001;
    % secant method to converge on radius of throat
    for i=1:stepLim
        if i < 3
            % ignore first 2 steps (guesses)
        else
            change = areaThroat(i-1) - areaThroat(i-2);
            if abs(change) < eps()
                convergeFlag = 1;  % flag for iterative convergence
            end

            areaDelta = (residual(i-1)*change) /(residual(i-1)-residual(i-2));
            areaThroat(i) = areaThroat(i-1) - (areaDelta/10);
        end

        area2star(i) = areaThroat(i) / aThroat_a2star;  % sonic area before normal shock
        area3star(i) = area2star(i) / stagPresRatio3;  % sonic area after normal shock
        aCombustor_a3star(i) = areaCombustor/area3star(i);
        [mach4(i), tempRatio4(i), presRatio4(i), densRatio4(i), ~] = ...
            flowisentropic(gamma, (aCombustor_a3star(i)), 'sub');  % ratios are static over stagnation
        staticDens4(i) = stagDens4*densRatio4(i);
        staticPres4(i) = stagPres4*presRatio4(i);
        staticTemp4(i) = stagTemp4*tempRatio4(i);
        residual(i) = pressureCombustor - staticPres4(i);
        if convergeFlag == 1
            ind = i;
            break
        end
        if abs(residual(i)) < tolerance
            ind = i;
            break
        end
        ind = i;
        i = i+1;
    end
    radiusThroat = sqrt((pi*radiusCowl^2 - areaThroat(ind))/pi);
    massflow3 = staticDens3*velocity3*areaThroat(ind);  % <kg/s> *************not the same as massflow4*********
    velocity4 = mach4(ind)*sqrt(gamma*R*staticTemp4(ind));  % <m/s>
    massflow4 = staticDens4(ind)*velocity4*areaCombustor;  % <kg/s>
    
    % plotting
    
    stagPres = [stagPres1 stagPres2 stagPres3 stagPres4];
    figure
    plot(stagPres)
    title('Stagnation Pressure')
    
    staticPres = [staticPres1 staticPres2 staticPres3 staticPres4(ind)];
    figure
    plot(staticPres)
    title('Static Pressure')
    
    stagDens = [stagDens1 stagDens2 stagDens3 stagDens4];
    figure
    plot(stagDens)
    title('Stagnation Density')    
    
    staticDens = [staticDens1 staticDens2 staticDens3 staticDens4(ind)];
    figure
    plot(staticDens)
    title('staticDens')
    
    stagTemp = [stagTemp1 stagTemp2 stagTemp3 stagTemp4];
    figure
    plot(stagTemp)
    title('stagTemp')
        
    staticTemp = [staticTemp1 staticTemp2 staticTemp3 staticTemp4(ind)];
    figure
    plot(staticTemp)
    title('staticTemp')
    
    mach = [mach1 mach2 mach3 mach4(ind)];
    figure
    plot(mach)
    title('Mach')
        
    velocity = [velocity1 velocity2 velocity3 velocity4];
    figure
    plot(velocity)
    title('Velocity')
        
    % iteration plots
%     plot(1:ind, areaThroat(1:ind))
%     title('areaThroat')
%     figure
%     plot(1:ind, area2star(1:ind))
%     title('area2star')
%     figure
%     plot(1:ind, area3star(1:ind))
%     title('area3star')
%     figure
%     plot(1:ind, mach4(1:ind))
%     title('mach4')
%     figure
%     plot(1:ind, tempRatio4(1:ind))
%     title('tempRatio4')
%     figure
%     plot(1:ind, presRatio4(1:ind))
%     title('presRatio4')
%     figure
%     plot(1:ind, densRatio4(1:ind))
%     title('densRatio4')
%     figure
%     plot(1:ind, staticDens4(1:ind))
%     title('staticDens4')
%     figure
%     plot(1:ind, staticPres4(1:ind))
%     title('staticPres4')
%     figure
%     plot(1:ind, staticTemp4(1:ind))
%     title('staticTemp4')
%     figure
%     plot(1:ind, abs(residual(1:ind)))
%     title('residual')
    
end

