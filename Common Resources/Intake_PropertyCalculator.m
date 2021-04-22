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
        
        function [intake] = Properties(self, mach1, mach4, altitude, defl)
            if exist('T','var')==0
                load GRAM_Model.mat
                T = GRAM;
            end
            
            % Station 1: Freestream
            intake.mach(1) = mach1;
            intake.staticDens(1) = interp1(T.Hgtkm, T.DensMean, altitude/1000);  % <kg/m3>
            intake.staticPres(1) = interp1(T.Hgtkm, T.PresMean, altitude/1000);   % <Pa>
            intake.staticTemp(1) = interp1(T.Hgtkm, T.Tmean, altitude/1000);       % <K>

            [~, tempRatio1, presRatio1, densRatio1, ~] = flowisentropic(self.gamma, intake.mach(1)); 
            intake.stagDens(1) = intake.staticDens(1) / densRatio1;
            intake.stagPres(1) = intake.staticPres(1) / presRatio1;
            intake.stagTemp(1) = intake.staticTemp(1) / tempRatio1;
            
            % Station 2: Oblique Shock
            [intake.mach(2), shockAngle] = obliqueShock(intake.mach(1), defl, self.gamma); 
            normalMach1 = intake.mach(1) * sind(shockAngle); 
            [~, tempRatio2, presRatio2, norm_densRatio2, ~, stagPresRatio2] = flownormalshock(self.gamma, normalMach1);
            [~, ~, ~, densRatio2, ~] = flowisentropic(self.gamma, intake.mach(2));
                        
            intake.staticDens(2)    = intake.staticDens(1) * norm_densRatio2; 
            intake.stagDens(2)      = intake.staticDens(2) / densRatio2;

            intake.staticTemp(2)    = intake.staticTemp(1) * tempRatio2;
            intake.stagTemp(2)      = intake.stagTemp(1);
            
            intake.staticPres(2)    = intake.staticPres(1) * presRatio2;
            intake.stagPres(2)      = intake.stagPres(1) * stagPresRatio2;
            
            % Station 3: Normal Shock
            [~, tempRatio3, presRatio3, norm_densRatio3, intake.mach(3), stagPresRatio3] = flownormalshock(self.gamma, intake.mach(2));
            [~, ~, ~, densRatio3, ~] = flowisentropic(self.gamma, intake.mach(3));
            
            intake.staticDens(3)    = intake.staticDens(2) * norm_densRatio3;
            intake.stagDens(3)      = intake.staticDens(3) / densRatio3;
            
            intake.staticTemp(3)    = intake.staticTemp(2) * tempRatio3;
            intake.stagTemp(3)      = intake.stagTemp(2);
            
            intake.staticPres(3)    = intake.staticPres(2) * presRatio3;
            intake.stagPres(3)      = intake.stagPres(2) * stagPresRatio3;
            
            
            % Station 4: Diffuser (mach 4 set by user)
            intake.mach(4) = mach4;
            [~, tempRatio4, presRatio4, densRatio4, ~] = flowisentropic(self.gamma, intake.mach(4));
            
            intake.stagDens(4)      = intake.stagDens(3);
            intake.staticDens(4)    = intake.stagDens(4) * densRatio4;
            
            intake.stagTemp(4)      = intake.stagTemp(3);
            intake.staticTemp(4)    = intake.stagTemp(4) * tempRatio4;
            
            intake.stagPres(4)      = intake.stagPres(3);            
            intake.staticPres(4)    = intake.stagPres(4) * presRatio4;
            
            % Solve for Pressure Recovery
            intake.PressureRecovery = intake.stagPres(4) / intake.stagPres(1);
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

