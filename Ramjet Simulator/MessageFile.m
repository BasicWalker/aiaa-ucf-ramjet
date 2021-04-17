% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: MessageFile.m 
%
% File Description: 
% Outputs simulation messages. 
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  04/04/21  ---  Initial creation 
% ---------------------------------------------------------------------- %

% Estimate Simulation Run Time
MaxSimSteps = (fuel.DiaOuter/2 - fuel.DiaInner(1)/2)/fuel.Regression + 1;
Status = 100 - (n/MaxSimSteps)*100;
if Status > 100
    Status = 100;
end

if Burnout == true
    fprintf('Running... Fuel Depleted...\n')
else
    fprintf('Running... Fuel Remaining: %.1f%%\n',Status)               % Running Simulator indicator
end
    
    