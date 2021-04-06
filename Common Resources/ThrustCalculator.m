classdef ThrustCalculator   
    properties     
    end
    
    methods
        function init = ThrustCalculator()
        end
        
        function thrust = ideal_Thrust(self,airMassFlow, f, exitVelocity, flightVelocity)
            thrust = airMassFlow*(((1+f)*exitVelocity) - flightVelocity);
        end
        
        function thrust = nonideal_Thrust(self,airMassFlow, f, exitVelocity, flightVelocity, Pe, Pa, Ae)
            thrust = airMassFlow*(((1+f)*exitVelocity) - flightVelocity) + ((Pe-Pa)*Ae);
        end
        
    end
end

