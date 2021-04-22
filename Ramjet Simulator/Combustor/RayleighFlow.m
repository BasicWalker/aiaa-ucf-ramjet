% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: RayleighFlow.m 
%
% File Description: 
% Simulates Combustion Process as Rayleigh Flow and Solves for Properties
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Karam Paul      03/26/21  000  Initial Creation 
% ---------------------------------------------------------------------- %

function [mach2, staticTemp2, staticPres2, stagTemp2, stagPres2, stagPresLoss] = ...
    RayleighFlow(g, mach1, stagTemp1, AFT, stagPres1)

    [~, ~, ~, ~, ~, stagTempRatRef1, stagPresRatRef1] = flowrayleigh(g, mach1, 'mach');
    
    stagTempRef1 = stagTemp1 / stagTempRatRef1;
    stagPresRef1 = stagPres1 / stagPresRatRef1;
    
    stagTemp2 = AFT; 
    
    stagTempRatRef2 = stagTemp2 / stagTempRef1;
    
    [mach2, ~, ~, ~, ~, ~, stagPresRatRef2] = flowrayleigh(g, stagTempRatRef2, 'totaltsub');
    
    stagPres2 = stagPresRef1 * stagPresRatRef2; 
    
    [~,tempRatio2, presRatio2, ~,~] = flowisentropic(g, mach2, 'mach');
    
    staticPres2 = stagPres2 * presRatio2; 
    
    staticTemp2 = stagTemp2 * tempRatio2; 
    
    stagPresLoss = stagPres1 - stagPres2; 
    
end

    
    
    
    
    
    
    
    
    
    
    

    

