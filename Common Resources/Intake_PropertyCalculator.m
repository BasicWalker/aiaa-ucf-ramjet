classdef Intake_PropertyCalculator
    
    properties
        gamma
        R
    end
    
    methods
        function init = Intake_PropertyCalculator()
            init.gamma  = 1.4;
            init.R      = 297;
        end
        
        function [initialStagPres, finalStagPres, stagPresLoss] = StagnationLoss(self, mach1, altitude, defl)
            if exist('T','var')==0
                load GRAM_Model.mat
                T = GRAM;
            end
            
            % Declare Freestream Properties
            staticPres1 = interp1(T.Hgtkm, T.PresMean, altitude/1000);  % <Pa>
            
            [~, ~, presRatio1, ~, ~] = flowisentropic(self.gamma, mach1, 'mach'); 
            stagPres1 = staticPres1 / presRatio1;
            initialStagPres = stagPres1;
            
            % Oblique Shock Procedure
            [mach2, shockAngle] = obliqueShock(mach1, defl, self.gamma); 
            normalMach1 = mach1 * sind(shockAngle); 

            [~, ~, ~, ~, ~, stagPresRatio2] = flownormalshock(self.gamma, normalMach1); 
            stagPres2 = stagPres1 * stagPresRatio2;
            
            %% Terminating Normal Shock Procedure 
            [~, ~, ~, ~, ~, stagPresRatio3] = flownormalshock(self.gamma, mach2); 
            stagPres3 = stagPres2 * stagPresRatio3;
            finalStagPres = stagPres3;
            
            % Calculate Stagnation Pressure Loss
            stagPresLoss = stagPres1 - stagPres3;
            
        end
        
        function velocity = Velocity(self, mach1, altitude)
             if exist('T','var')==0
                load GRAM_Model.mat
                T = GRAM;
            end
            
            % Declare Freestream Properties
            staticTemp1 = interp1(T.Hgtkm, T.Tmean, altitude/1000);  % <K>
            velocity = mach1*sqrt(self.gamma * self.R * staticTemp1);
        end
    end
end

