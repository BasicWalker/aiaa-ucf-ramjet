%% Variable Declaration
clc; clear all; close all;

gamma = 1.4; 
mach1 = 2; 
n = 11;          % n_max is 11 for mach 2 and 1 deflection 
%n = 6;         % n_max is 6 for mach 2 and 2 deflections
%n = 4;         % n_max is 4 for mach 2 and 3 deflections           
k = 1;           % set the number of oblique shocks

defl = zeros(n*2, 1); 
mach2 = zeros(n*2, 1);
mach3 = zeros(n*2, 1);
mach4 = zeros(n*2, 1);

theta1 = zeros(n*2, 1);
theta2 = zeros(n*2, 1);
theta3 = zeros(n*2, 1);


normalMach1 = zeros(n*2, 1);
normalMach2 = zeros(n*2, 1);
normalMach3 = zeros(n*2, 1); 
normalMach4 = zeros(n*2, 1); 
normalMach5 = zeros(n*2, 1); 


staticTempRatio_21 = zeros(n*2, 1);
staticTempRatio_32 = zeros(n*2, 1);
staticTempRatio_43 = zeros(n*2, 1);
staticTempRatio_54 = zeros(n*2, 1);


staticPresRatio_21 = zeros(n*2, 1);
staticPresRatio_32 = zeros(n*2, 1);
staticPresRatio_43 = zeros(n*2, 1);
staticPresRatio_54 = zeros(n*2, 1);

stagPresRatio_21 = zeros(n*2, 1);
stagPresRatio_32 = zeros(n*2, 1);
stagPresRatio_43 = zeros(n*2, 1);
stagPresRatio_54 = zeros(n*2, 1);

presRecoveryFactor = zeros(n*2,1);

%% Solve for isentropic pressure ratio
[~, ~, presRatio_1, ~, ~] = flowisentropic(gamma, mach1, 'mach');

%% Begin calculations 

switch k
    case 1: 1
        i = 1;
        % oblique 1
        for deflectionAngle = 1 : 0.5 : 22.5
            defl(i,1) = deflectionAngle;
            [mach2(i), theta1(i)] = obliqueShockCalculator(mach1, deflectionAngle, gamma);
            if mach2(i) ~= real(mach2(i))
                mach2(i) = NaN;
                theta1(i) = NaN;
            end

            % normal 1
            normalMach1(i) = mach1 * sind(theta1(i));
            [staticTempRatio_21(i), staticPresRatio_21(i), normalMach2(i), stagPresRatio_21(i)]  ...
                = normalShockCalculator(gamma, normalMach1(i));

            % normal 2
            [staticTempRatio_32, staticPresRatio_32(i), normalMach3(i), stagPresRatio_32(i)] ...
            = normalShockCalculator(gamma, mach2(i));

            presRecoveryFactor(i) = stagPresRatio_32(i) * stagPresRatio_21(i);

            i = i + 1;
        end
        varNames = {'Deflection Angle', 'Shock Angle', 'Recovery Factor'}; 
        T = table(defl, theta1, presRecoveryFactor);
        T.Properties.VariableNames = varNames;
        
        
    case 2: 2
        i = 1;
        for deflectionAngle = 1 : 0.5 : 12.5
            defl(i,1) = deflectionAngle;
            
            % Oblique 1
            [mach2(i), theta1(i)] = obliqueShockCalculator(mach1, deflectionAngle, gamma);
            if mach2(i) ~= real(mach2(i))
                mach2(i) = NaN;
                theta1(i) = NaN;
            end

            % Normal 1
            normalMach1(i) = mach1 * sind(theta1(i));
            [staticTempRatio_21(i), staticPresRatio_21(i), normalMach2(i), stagPresRatio_21(i)]  ...
                = normalShockCalculator(gamma, normalMach1(i));
            
            % Oblique 2
            [mach3(i), theta2(i)] = obliqueShockCalculator(mach2(i), deflectionAngle, gamma);
            if mach3(i) ~= real(mach3(i))
                mach3(i) = NaN;
                theta2(i) = NaN;
            end
            
            % Normal 2
            normalMach2(i) = mach2(i) * sind(theta2(i)); 
            [staticTempRatio_32(i), staticPresRatio_32(i), normalMach3(i), stagPresRatio_32(i)]  ...
                = normalShockCalculator(gamma, normalMach2(i));
            
            % Normal 3
            [staticTempRatio_43, staticPresRatio_43(i), normalMach4(i), stagPresRatio_43(i)] ...
            = normalShockCalculator(gamma, mach3(i));

            presRecoveryFactor(i) = stagPresRatio_43(i) * stagPresRatio_32(i) * stagPresRatio_21(i);

            i = i + 1;
        end
        varNames = {'Deflection Angle1', 'Shock Angle1', 'DeflectionAngle2', ...
            'ShockAngle2', 'Recovery Factor'}; 
        T = table(defl, theta1, defl, theta2, presRecoveryFactor);
        T.Properties.VariableNames = varNames;
        
    case 3: 3
        i = 1;
        for deflectionAngle = 1 : 0.5 : 8.5    
            defl(i,1) = deflectionAngle;
            
            % Oblique 1
            [mach2(i), theta1(i)] = obliqueShockCalculator(mach1, deflectionAngle, gamma);
            if mach2(i) ~= real(mach2(i))
                mach2(i) = NaN;
                theta1(i) = NaN;
            end

            % Normal 1
            normalMach1(i) = mach1 * sind(theta1(i));
            [staticTempRatio_21(i), staticPresRatio_21(i), normalMach2(i), stagPresRatio_21(i)]  ...
                = normalShockCalculator(gamma, normalMach1(i));
            
            % Oblique 2
            [mach3(i), theta2(i)] = obliqueShockCalculator(mach2(i), deflectionAngle, gamma);
            if mach3(i) ~= real(mach3(i))
                mach3(i) = NaN;
                theta2(i) = NaN;
            end
            
            % Normal 2
            normalMach2(i) = mach2(i) * sind(theta2(i)); 
            [staticTempRatio_32(i), staticPresRatio_32(i), normalMach3(i), stagPresRatio_32(i)]  ...
                = normalShockCalculator(gamma, normalMach2(i));
            
            % Oblique 3
            [mach4(i), theta3(i)] = obliqueShockCalculator(mach3(i), deflectionAngle, gamma);
            if mach4(i) ~= real(mach4(i))
                mach3(i) = NaN;
                theta2(i) = NaN;
            end
            
            % Normal 3
            normalMach3(i) = mach3(i) * sind(theta3(i)); 
            [staticTempRatio_43(i), staticPresRatio_43(i), normalMach4(i), stagPresRatio_43(i)]  ...
                = normalShockCalculator(gamma, normalMach3(i));
            
            % Normal 4
            [staticTempRatio_54, staticPresRatio_54(i), normalMach5(i), stagPresRatio_54(i)] ...
            = normalShockCalculator(gamma, mach4(i));

            presRecoveryFactor(i) = stagPresRatio_54(i) * stagPresRatio_43(i) * stagPresRatio_32(i) * stagPresRatio_21(i);

            i = i + 1;
        end
        varNames = {'Deflection Angle1', 'Shock Angle1', 'DeflectionAngle2', ...
            'ShockAngle2', 'DeflectionAngle3', 'ShockAngle3', 'Recovery Factor'}; 
        T = table(defl, theta1, defl, theta2, defl, theta3, presRecoveryFactor);
        T.Properties.VariableNames = varNames;
end 
     
% IF WE WANT TO MAKE TABLES BUT ID RECOMMEND EXPORTING TO EXCEL
% varNames = {'Deflection Angle', 'Shock Angle', 'Recovery Factor'}; 
% T = table(defl, theta, presRecoveryFactor);
% T.Properties.VariableNames = varNames;








