% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: RayleighFlow.m 
%
% File Description: 
% Simulates Combustion Process as Rayleigh Flow and Solves for Properties
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Karam Paul      04/26/21  000  Initial Creation 
% ---------------------------------------------------------------------- %

function [mach2, staticTemp2, staticPres2, stagTemp2, stagPres2, stagPresLoss] = ...
    RayleighFlow(gamma, Cp, mach1, staticTemp1, staticPres1, q)

    [~, isenTempRatio1, isenPresRatio1, ~, ~] = flowisentropic(gamma, mach1, 'mach');
    
    stagTemp1 = staticTemp1 / isenTempRatio1;
    
    stagPres1 = staticPres1 / isenPresRatio1;                             
        
    stagTemp2 = stagTemp1 + q/Cp;
    
    [~, rayStaticTempRatio1, rayStaticPresRatio1, ~, ~, rayStagTempRatio1, ~] ...  
        = flowrayleigh(gamma, mach1, 'mach');
    
    rayStagTempRatio2 = (stagTemp2/stagTemp1) * rayStagTempRatio1;
    
    [mach2, rayStaticTempRatio2, rayStaticPresRatio2, ~, ~, ~, ~] ...
        = flowrayleigh(gamma, rayStagTempRatio2, 'totaltsub');
    
    [~, ~, isenPresRatio2, ~, ~] = flowisentropic(gamma, mach2, 'mach');

    staticTemp2 = rayStaticTempRatio2 * (1/rayStaticTempRatio1) * staticTemp1;
    
    staticPres2 = rayStaticPresRatio2 * (1/rayStaticPresRatio1) * staticPres1;
    
    stagPres2   = (1/isenPresRatio2) * staticPres2; 
    
    stagPresLoss = stagPres1 - stagPres2; 
    
end
    
    
    
    
    
    
    
    
    
    
    

    

