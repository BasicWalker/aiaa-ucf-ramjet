function [mach2, shockAngle] = obliqueShock(mach1, deflectionAngle, gamma)
%obliqueShock Outputs the weak shock shock angle and downstream mach number
%   inputs:
%       mach1 (unitless) mach number of flow before shock
%       deflectionAngle (degrees) half angle from ramp to horizontal
%       gamma (unitless) specific heat ratio of fluid (optional input;
%           default va1ue of 1.4)
%   method:
%       Collars method converts equation 6.18 to a cubic equation of the form  
%       x^3 + C*x^2 - A*x + (B - A*C) == 0, where x = cotd(theta).
%       iterate until convergence of theta to find the largest root this
%       theta (weak shock) is then used to find upstream mach with eq. 6.17
%       
%   references from "Gas Dynamics" third edition James E. John

% check if gamma is input
 if ~exist('gamma','var')
    % third parameter does not exist, so default it to something
    gamma = 1.4;
 end

 
% Collar's method coefficients
A = mach1^2 - 1;
B = ((gamma + 1)/2)*mach1^4*tand(deflectionAngle);
C = (1 + ((gamma + 1)/2)*mach1^2)*tand(deflectionAngle);
 
 itr = 25;  % max iteration
 x_1 = sqrt(A);  % initial guess
 for i = 1:itr
     x = sqrt(A - B/(x_1+C));  % re-arrange cubic function to solve for x
     if abs(x - x_1) < eps()
         break
     end
     x_1 = x;
 end
 
theta = acotd(x);  % solve for theta based on iteration result

shockAngle = theta;

% equation 6.17 downstream mach number
mach2 = sqrt( (1+(gamma-1)/2*mach1^2)/((gamma*mach1^2*sind(theta)^2)- (gamma - 1)/2)...
    + ( mach1^2*cosd(theta)^2/(1 + (gamma - 1)/2*mach1^2*sind(theta)^2) )) ;
if ~isreal(mach2)
    error('aero:obliqueDetach',...
        'Error. \nSolution invalid; oblique shock detached. Try lowering deflection angle');
end
end

