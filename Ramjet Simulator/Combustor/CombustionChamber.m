% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: CombustionChamber.m 
%
% File Description: 
% Simulates Combustion Chamber processes within the Ramjet to define flow
% properties
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Karam Paul      04/26/21  000  Initial Creation 
% ---------------------------------------------------------------------- %    

RegressionRate                              % Call Regression Rate Model
GrainGeometry                               % Call Instantaneous Grain Geometry Model
Gas                                         % Call Gas Model (And Chemistry Model)
BoundaryLayer                               % Call Boundary Layer Model

% calculate the expansion for step height
% cmbMach2*GrainID(n) = InltD


[~, ~, ~, ~, ccA_Astar(1,n)] = flowisentropic(gamma, InltMach(n), 'mach');
Astar(1,n) = InltArea / ccA_Astar(1,n);
ccA_Astar(2,n) = PortArea(n) / Astar(1,n);
% 
% 
[ccMach(1,n), ccTempRatio(1,n), ccPresRatio(1,n), ccDensRatio(1,n), ~] = ...
    flowisentropic(gamma, ccA_Astar(2,n), 'sub');  % ratios are static over stagnation
ccStagPres(1,n) = InltPres_stag(n);  % does not change; iscentropic <Pa>
ccStagDens(1,n) = InltTemp_dens(n);  % does not change; iscentropic <kg/m3>
ccStagTemp(1,n) = InltTemp_stag(n);  % does not change; iscentropi <K>
ccStaticPres(1,n) = ccStagPres(1,n) * ccPresRatio(1,n);  % <Pa>
ccStaticDens(1,n) = ccStagDens(1,n) * ccDensRatio(1,n);  % <kg/m3>
ccStaticTemp(1,n) = ccStagTemp(1,n) * ccTempRatio(1,n);  % <K>
ccVelocity(1,n) = ccMach(1,n)*sqrt(gamma*R*ccStaticTemp(1,n));  % <m/s>
ccMdot(1,n) = ccStaticDens(1,n)*ccVelocity(1,n)*PortArea(n);
% 
% use rayleigh flow to solve for properties at combustion chamber exit
[ccMach(2,n), ccStaticTemp(2,n), ccStaticPres(2,n), ccStagTemp(2,n), ccStagPres(2,n), ccStagPresLoss(n)] = ...
    RayleighFlow(1.2, ccMach(1,n), InltTemp_stag(n), T_stag(n), InltPres_stag(n));

ccStaticDens(2,n) = ccStaticPres(2,n)/(ccStaticTemp(2,n)*R);
ccStagDens(2,n) = ccStagPres(2,n)/(ccStagTemp(2,n)*R);
ccVelocity(2,n) = ccMach(2,n)*sqrt(1.2*600*ccStaticTemp(2,n));
ccMdot(2,n) = ccStaticDens(2,n)*ccVelocity(2,n)*PortArea(n);



