% ---------- SFRJ Numerical Simulation / UCF CAPSTONE PROJECT ---------- %
% Program Name:  SFRJ Internal Ballistic Simulator
%
% Program Description:
% This program models and simulates critical performance parameters of a
% Solid-Fuel Ramjet.
%
% File Description:
% Main executive file.
%
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial creation
% Ethan Sherlock  02/14/21  005  1DOF Trajectory Initialization
% Samer & karam   --------  ---  Sim Revamp
% Ethan Sherlock  04/14/21  ---  2DOF Trajectory Update
% ---------------------------------------------------------------------- %
clear; clc; close all;
addpath(genpath(pwd))
addpath('..\Common Resources')
addpath(genpath('..\Ramjet Simulator'))
% Import data
load GRAM_Model.mat  % GRAM atmospheric model
load Constants.mat  % load in constants and conversions

% Initialize Simulation Variables
SFRJDt              = 1/100;                                % Simulation rate (Hz)
n                   = 1;                                    % Initialize counter
time                = 0.0;                                  % Initialize time (s)
StopBurn            = false;                                % Burn status flag (boolean)
Burnout             = false;                                % Burn out status flag (boolean)

cnt = 10;  % iteration count

% nozzle -----
Design.DiaThroatINCHiter = linspace(1.5,10,cnt);  % Throat Diameter <in> 1.8  % iteration parameter ****
Design.AreaRatioExititer = linspace(1.5,5,cnt); % iteration parameter ****

% combustion -----
% combustion.InletDiaINCHiter = linspace(2,5,cnt);  % <in>  % iteration parameter ****



% intake -----
Design.AreaRatioEnteriter = linspace(1.5,10,cnt);  % iteration parameter ****

% Design.EnterDiaiter = linspace(0.5,5,cnt);  % iteration parameter ****
intake.DeflAngle = 15;   % Deflection angle (deg)

% fuel -----

fuel.Length = 15.00 /constants.In2Mtr;  % Grain Length (m)
fuel.Density = 1020;  % Grain Density (kg/m^3)

% Ramjet Dimensions
vehicle.DragCoeff = 0.23;  % Drag coefficient (0.35)
vehicle.DryMass = 6.80389;  % Mass of ramjet without fuelgrain (kg)
burntimedesign = 5;  %<s>
fuel.BurnRate     = 0.001;
fuel.Regression = fuel.BurnRate*SFRJDt; % Regression Per Step (m)




% User Defined Parameters
vehicle.Mach(1)         = 2.5;  % Ramjet Initial Mach number
trajectory.Z_pos(1)     = 1000.0;  % Initial altitude for ramjet start (m)
alpha                   = 0.0;  % Launch Angle (deg) - in reference to horizon
trajectory.LiftOnOff    = 1;  % 0.0 = off; % 1.0 = on

% -------------------- Commence Ramjet trade study  --------------------- %

for i=1:numel(Design.DiaThroatINCHiter)
    for j=1:numel(Design.AreaRatioExititer)
        for k=1:numel(Design.AreaRatioEnteriter) % Design.EnterDiaiter
            n=1;
            
            nozzle.DiaThroatINCH = Design.DiaThroatINCHiter(i);
            nozzle.AreaRatioExit = Design.AreaRatioExititer(j);
%             intake.EnterDiaINCH = Design.EnterDiaiter(k);  % <in>
            

            
            nozzle.DiaThroat = nozzle.DiaThroatINCH/constants.In2Mtr;  % <m>
            nozzle.Area_throat = pi*(nozzle.DiaThroat)^2/4;  % Throat area (m^2)
            nozzle.Area_exit = nozzle.Area_throat*nozzle.AreaRatioExit;  % nozzle exit area (m^2)
            nozzle.DiaExit = sqrt(nozzle.Area_exit*4/pi);  % nozzle exit area (m)
            nozzle.DiaExitINCH = nozzle.DiaExit*constants.In2Mtr;  % nozzle exit area (in)
            
            intake.Area_enter = nozzle.Area_throat/Design.AreaRatioEnteriter(k);  % <m^2>
            
%             intake.EnterDia = intake.EnterDiaINCH/constants.In2Mtr;  % <m>
%             intake.Area_enter = pi*intake.EnterDia^2/4;  % Area of throat (m^2) - Drives mass flow rate through intake
            
            
            trajectory.Rho_a(1)     = interp1(GRAM.Hgtkm, GRAM.DensMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Density (kg/m^3)
            trajectory.pressure_a(1)= interp1(GRAM.Hgtkm, GRAM.PresMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Pressure (Pa)
            trajectory.Temp_a(1)    = interp1(GRAM.Hgtkm, GRAM.Tmean, (trajectory.Z_pos(1))/1e3);                   % Atmospheric Temperature (K)
%             
%             n=1;
%             % find area needed to slow com bustion chamber to mach 0.15
%             designmach(1) = vehicle.Mach(n);
%             designstaticPres(1) = trajectory.pressure_a(n);  % <Pa>
%             designstaticDens(1) = trajectory.Rho_a(n);  % <kg/m^3>
%             designstaticTemp(1) = trajectory.Temp_a(n);  % <K>
%             [~, designtempRatio(1), designpresRatio(1), designdensRatio(1), ~] = flowisentropic(constants.gamma, designmach(1), 'mach');  % ratios are static over stagnation
%             designstagTemp(1) = designstaticTemp(1)/designtempRatio(1);  % <K>
%             designstagPres(1) = designstaticPres(1)/designpresRatio(1);  % <Pa>
%             designstagDens(1) = designstaticDens(1)/designdensRatio(1);  % <kg/m^3>
%             designvelocity(1) = designmach(1)*sqrt(constants.gamma*constants.R*designstaticTemp(1));  % <m/s>
%             % Station 2 properties (oblique shock)-----
%             [designmach(2), designshockAngle(n)] = obliqueShock(designmach(1), intake.DeflAngle, constants.gamma);  % uniform mach number of external ramp region
%             designmach1Norm(n) = designmach(1)*sind(designshockAngle(n));  % flow component crossing normal to oblique shock
%             [~, designtempRatio(2), designpresRatio(2), designdensRatio(2), ~, designstagPresRatio(2)] = ...
%                 flownormalshock(constants.gamma, designmach1Norm(n), 'mach');  % ratios are downstream over upstream
%             designstagTemp(2) = designstagTemp(1);  % does not change over normal shock <K>
%             designstagPres(2) = designstagPres(1) * designstagPresRatio(2);  % <Pa>
%             designstaticDens(2) = designstaticDens(1) * designdensRatio(2);  % <kg/m3>
%             designstaticTemp(2) = designstaticTemp(1) * designtempRatio(2);  % <K>
%             designvelocity(2) = designmach(2)*sqrt(constants.gamma*constants.R*designstaticTemp(2));  % <m/s>
%             [~, ~, ~, ~, designA_Astar(2)] = flowisentropic(constants.gamma, designmach(2), 'mach');
%             designAstar(2) = intake.Area_enter/designA_Astar(2);  % reference area before normal shock
%             designmassFlow(2) = designstaticDens(2)*intake.Area_enter*designvelocity(2);  % air mass flow at intake opening <kg/s>
%             designchokeStagPres(n) = pressureToChoke(designmassFlow(2), nozzle.Area_throat, designstagTemp(2));  % no stag temp change from combustion
%             designr_n = designchokeStagPres(n)/designstagPres(2);  % total stagnation pressure loss ratio from normal shock
%             [~, ~, ~, ~, ~, designr_n_chokeLimit(n)] = flownormalshock(constants.gamma, designmach(2), 'mach');
%             designAstar(4) = designAstar(2)/designr_n_chokeLimit(n);
%             ideal_combustor_area =  3.91034275*designAstar(4);

%             Initialize  % set initial conditions
% Initialize BurnTime variable (s)
                        BurnTime(1)         = 0.0;
            StopBurn = 0;
            while StopBurn == 0
                if n > 3  % stop simulation at 3rd iteration
                    % add variables to track
                    
                    Design.Thrust(i,j,k) = vehicle.thrust(3);
                    Design.intakemassFlow(i,j,k) = intake.massFlow(2,3);
                    Design.MaxAvailableStag(i,j,k) = intake.MaxAvailableStag(3);
                    Design.chokeStagPres(i,j,k) = intake.chokeStagPres(3);
                    Design.combustionStagPres(i,j,k) = combustion.stagPres(2,3);
                    
                    Design.Thrustgood(i,j,k) = 1;
                    Design.intakemassFlowgood(i,j,k) = 1;
                    Design.MaxAvailableStaggood(i,j,k) = 1;
                    Design.chokeStagPresgood(i,j,k) = 1;
                    Design.combustionStagPresgood(i,j,k) = 1;
                    break
                end
                BurnTime(n) = time;                                     % Simulation Time
                if Burnout == 0
                    try
%                         RegressionRate                                      % Call Regression Rate Model
                        IntakeDesign                                              % Call Intake Model
                        design.EnterHeight(i,j,k) = intake.EnterDiaINCH;  % <in>
                        design.CowlHeight(i,j,k) = intake.CowlDiaINCH;  % <in>
                        design.IntakeArea(i,j,k) = intake.Area_enter*constants.In2Mtr^2;  %  <in^2>
                        design.ThroatArea(i,j,k) = nozzle.Area_throat*constants.In2Mtr^2;  %  <in^2>



                        % Initialize Fuel Grain Model
%                         fuel.PortArea(1)         = pi*(fuel.DiaInner(1)^2)*(1/4);  % Fuel Port Area (m^2)
                        fuel.CsxArea(1)           = pi*(fuel.DiaOuter^2)*(1/4) - fuel.PortArea(1);   % Fuel Grain Crossectional Area (m^2)
                        fuel.Volume(1)          = fuel.CsxArea(1) * fuel.Length;  % Fuel Grain Volume (m^3)
                        fuel.SurfArea(1)           = fuel.DiaInner(1)* pi * fuel.Length;  % Fuel Grain Surface Area (m^2)
                        fuel.Mass(1)         = fuel.Density*fuel.Volume(1);  % Grain fuel mass, instantaneous (kg)
                        vehicle.TotalMass(1) = vehicle.DryMass + fuel.Mass(1);  % Mass of Vehicle (Kg)
                        
                        % Initialize Chemistry Model
                        chem = Chemistry();
                        
                        % Initialize Trajectory Model
                        vehicle.thrust(1)               = 0.0;                                                                  % Thrust (N)
                        trajectory.Lift(1)      = vehicle.TotalMass(1)*constants.gravity*trajectory.LiftOnOff;                                 % Lift (N)
                        trajectory.Rho_a(1)     = interp1(GRAM.Hgtkm, GRAM.DensMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Density (kg/m^3)
                        trajectory.pressure_a(1)= interp1(GRAM.Hgtkm, GRAM.PresMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Pressure (Pa)
                        trajectory.Temp_a(1)    = interp1(GRAM.Hgtkm, GRAM.Tmean, (trajectory.Z_pos(1))/1e3);                   % Atmospheric Temperature (K)
                        trajectory.Vel(1)       = vehicle.Mach(1)*sqrt(constants.gamma*constants.R*trajectory.Temp_a(1));                             % Velocity (m/s)
                        trajectory.Vel_x(1)     = trajectory.Vel(1)*cosd(alpha);                                        % Velocity X (m/s)
                        trajectory.Vel_z(1)     = trajectory.Vel(1)*sind(alpha);                                        % Velocity Z (m/s)
                        trajectory.F_d(1)       = 0.5*trajectory.Rho_a(1)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel(1)^2;                    % Drag (N)
                        trajectory.F_dx(1)      = 0.5*trajectory.Rho_a(1)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel(1)*trajectory.Vel_x(1);  % Drag X (N)
                        trajectory.F_dz(1)      = 0.5*trajectory.Rho_a(1)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel(1)*trajectory.Vel_z(1);  % Drag Z (N)
                        trajectory.F_tx(1)      = vehicle.thrust(1)*cosd(alpha);                                            % Thrust X (N)
                        trajectory.F_tz(1)      = vehicle.thrust(1)*sind(alpha);                                            % Thrust Z (N)
                        trajectory.F_x(1)       = trajectory.F_tx(1)-trajectory.F_dx(1);                                    % Force X (N)
                        trajectory.F_z(1)       = trajectory.F_tz(1)-trajectory.F_dz(1)-vehicle.TotalMass(1)*constants.gravity+trajectory.Lift(1); % Force Z (N)
                        trajectory.F_net(1)     = sqrt(trajectory.F_x(1)^2+trajectory.F_z(1)^2);                        % Force (N)
                        trajectory.Acc(1)       = trajectory.F_net(1)/vehicle.TotalMass(1);                                          % Acceleration (m/s/s)
                        trajectory.Acc_x(1)     = trajectory.F_x(1)/vehicle.TotalMass(1);                                            % Acceleration X (m/s/s)
                        trajectory.Acc_z(1)     = trajectory.F_z(1)/vehicle.TotalMass(1);                                            % Acceleration Z (m/s/s)
                        trajectory.X_pos(1)     = 0.0;                                                                  % Initial X Position (m)
                        trajectory.Vel_t(1)     = sqrt((2*vehicle.TotalMass(1)*constants.gravity)/(vehicle.DragCoeff*trajectory.Rho_a(1)*vehicle.FrontSurfArea));                % Terminal Velocity (m/s)
                        
                        
                        
                        
                        
                        GrainGeometry                                       % Call Instantaneous Grain Geometry Model
                        Gas                                                 % Call Gas Model (And Chemistry Model)
                        CombustionChamber                                   % Call Combustion Chamber Model
                        Nozzle                                              % Call Nozzle Model
                        Thrust                                              % Call Thrust Model
                        if StopBurn == 1                                    % do not update trajectory since fuel is depleted
                            break
                        end
                        Trajectory                                          % Call Trajectory Model
                    catch ME
                        fprintf('something went wrong:( %s\n', ME.message);
                        try
                            Design.Thrust(i,j,k) = vehicle.thrust(n);
                        catch
                            Design.Thrust(i,j,k) = NaN;
                        end
                        try
                            Design.intakemassFlow(i,j,k) = intake.massFlow(2,n);
                        catch
                            Design.intakemassFlow(i,j,k) = NaN;
                        end
                        try
                            Design.MaxAvailableStag(i,j,k) = intake.MaxAvailableStag(n);
                        catch
                            Design.MaxAvailableStag(i,j,k) = NaN;
                        end
                        try
                            Design.chokeStagPres(i,j,k) = intake.chokeStagPres(n);
                        catch
                            Design.chokeStagPres(i,j,k) = NaN;
                        end
                        try
                            Design.combustionStagPres(i,j,k) = combustion.stagPres(2,n);
                        catch
                            Design.combustionStagPres(i,j,k) = NaN;
                        end
                        Design.Thrustgood(i,j,k) = 0;
                        Design.intakemassFlowgood(i,j,k) = 0;
                        Design.MaxAvailableStaggood(i,j,k) = 0;
                        Design.chokeStagPresgood(i,j,k) = 0;
                        Design.combustionStagPresgood(i,j,k) = 0;
                        break
                    end
                end
                time = time + SFRJDt;                                   % Step through simulation time
                n = n + 1;                                              % Increase Index
            end
            
        end
    end
end
clear gamma R
%% plot trade study results

% for i=1:numel(Design.DiaThroatINCHiter)
%     for j=1:numel(Design.AreaRatioExititer)
%         for k=1:numel(Design.EnterDiaiter)


figure()
plot(Design.AreaRatioEnteriter,squeeze(Design.Thrust(1,1,:)));
xlabel('Entrance Diameter <in>')
ylabel('Thrust <n>')
title('entrance dia thrust')

figure()
plot(Design.AreaRatioEnteriter,squeeze(Design.MaxAvailableStag(1,1,:))/1e3,Design.AreaRatioEnteriter,squeeze(Design.chokeStagPres(1,1,:))/1e3);
xlabel('Entrance Diameter <in>')
ylabel('stagnation pressure <kPa>')
title('entrance dia pressures')
legend('available stagnation', 'choke stagnation required')

figure()
plot(Design.DiaThroatINCHiter,squeeze(Design.Thrust(:,1,1)));
xlabel('Entrance Diameter <in>')
ylabel('Thrust <n>')
title('throat diameter thrust')

figure()
plot(Design.DiaThroatINCHiter,squeeze(Design.MaxAvailableStag(:,1,1))/1e3,Design.DiaThroatINCHiter,squeeze(Design.chokeStagPres(:,1,1))/1e3);
xlabel('Diameter <in>')
ylabel('stagnation pressure <kPa>')
title('throat diameter pressures')
legend('available stagnation', 'choke stagnation required')






% entrance Dia vs nozzle throat

[X,Y] = meshgrid(Design.AreaRatioEnteriter,Design.DiaThroatINCHiter);
xname = 'Entrance Diameter <in>';
yname = 'throat Diameter <in>';
group = append(xname,' vs ',yname);
figure()
surf(X,Y,squeeze(Design.Thrust(:,1,:)));
xlabel(xname)
ylabel(yname)
title('thrust')
figure()
surf(X,Y,squeeze(Design.intakemassFlow(:,1,:)));
xlabel(xname)
ylabel(yname)
title('massflow')
figure()
surf(X,Y,squeeze(Design.MaxAvailableStag(:,1,:)),'FaceAlpha',0.2,'FaceColor','r');
hold on
surf(X,Y,squeeze(Design.chokeStagPres(:,1,:)),'FaceColor','b');
xlabel(xname)
ylabel(yname)
title('Pressures')
legend('available stagnation', 'choke stagnation required')
hold off
figure()
surf(X,Y,squeeze(Design.combustionStagPres(:,1,:)));
xlabel(xname)
ylabel(yname)
title('combustion stag')
% 
% Design.EnterDiaiter
% % nozzle throat vs nozzle Area ratio --------
% [X,Y] = meshgrid(Design.DiaThroatINCHiter,Design.AreaRatioExititer);
% xname = 'throat Diameter <in>';
% yname = 'nozzle area ratio';
% group = append(xname,' vs ',yname);
% figure()
% surf(X,Y,Design.Thrust(:,:,1));
% xlabel(xname)
% ylabel(yname)
% title('thrust')
% figure()
% surf(X,Y,Design.intakemassFlow(:,:,1));
% xlabel(xname)
% ylabel(yname)
% title('massflow')
% figure()
% surf(X,Y,Design.MaxAvailableStag(:,:,1));
% xlabel(xname)
% ylabel(yname)
% title(append('available Stag')
% figure()
% surf(X,Y,Design.chokeStagPres(:,:,1));
% xlabel(xname)
% ylabel(yname)
% title('choke pressure')
% figure()
% surf(X,Y,Design.combustionStagPres(:,:,1));
% xlabel(xname)
% ylabel(yname)
% title('combustion stag')
%

