%% --------- Senior Design - Ramjet Powered Vehicle --------- %
% Program Name: Normal Shock Placement 

% Program Description: 
%
% 
% File Name: Normalshcok_Placement.m
% 
% File Description: 
% 
% Name            Date      Description
% --------------  --------  ------------------------------
% Samer           03/1/21  Initial Creation 
% --------------------------------------------------------------------- %  
gamma = 1.4;
pressureCombustor = 500e3;  % combustion chamber pressure <kPa>
intake_stagPressure = 400e4;  % intake.stagPres(1,n
r_d = pressureCombustor/intake_stagPressure;  % stag pressure ratio needed across the normal shock

stepLim = 1000;  % convergence iteration step limit
tolerance = 0.0001;  % tolerance for residual
residual = zeros(stepLim,1);
mach_r_d = zeros(stepLim,1);
convergeFlag = 0;
ind = 0;
mach_r_d(1) = 2.5;  %intake.mach(2,n);  % first guess
mach_r_d(2) = mach_r_d(1) + 0.001;  % second guess

% find mach number needed to create proper strength of normal shock
for j=1:stepLim
    if j < 3
        % ignore first 2 steps (guesses)
    else
        change = mach_r_d(j-1) - mach_r_d(j-2);
        if abs(change) < eps()
            convergeFlag = 1;  % flag for iterative convergence
        end
        
        Delta = (residual(j-1)*change) / (residual(j-1)-residual(j-2));
        mach_r_d(j) = mach_r_d(j-1) - (Delta);  % divided by 10 to slow down convergence
    end
    [~, ~, ~, ~, ~, r_d_j(j)] = flownormalshock(gamma, mach_r_d(j), 'mach');
    residual(j) = r_d - r_d_j(j);
    if convergeFlag == 1
        ind = j;
        break
    end
    if abs(residual(j)) < tolerance
        ind = j;
        break
    end
    ind = j;
    j = j+1;
end

