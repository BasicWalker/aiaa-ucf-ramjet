function [mach2] = normalShock(mach1, gamma)
%obliqueShock Outputs the weak shock angle and downstream mach number
%   inputs:
%       mach1 (unitless) mach number of flow before shock
%       gamma (unitless) specific heat ratio of fluid (optional input;
%           default va1ue of 1.4)
%   method:
%       Collars method converts equation 6.18 to a cubic equation of the form  
%       x^3 + C*x^2 - A*x + (B - A*C) == 0, where x = cotd(theta).
%       iterate until convergence of theta to find the largest root. this
%       theta (weak shock) is then used to find downstream mach with eq. 6.17
%       
%   references from "Gas Dynamics" third edition James E. John

% check if gamma is input
 if ~exist('gamma','var')
    % third parameter does not exist, so default it to something
    gamma = 1.4;
 end

% equation 4.9 flow across normal shock
mach2 = sqrt( ( mach1^2 + 2/(gamma - 1) )/...
    ( (2*gamma/(gamma - 1)*mach1^2 - 1) ) );
end

