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


preserve_data = 0;  % boolean to try to recover data if iteration failed, else data stored as NaN
cnt = 10;  % iteration count

% nozzle -----
Design.DiaThroatINCHiter = linspace(1.5,3,cnt);  % Throat Diameter <in> 1.8  % iteration parameter ****
% Design.DiaThroatINCHiter = linspace(0.5,5,cnt);  % Throat Diameter <in> 1.8  % iteration parameter ****
Design.AreaRatioExititer = linspace(1,5,cnt); % iteration parameter ****
% Design.AreaRatioExititer = linspace(1.5,5,cnt); % iteration parameter ****

% intake -----
Design.AreaRatioEnteriter = linspace(1,5,cnt);  % iteration parameter ****
% Design.AreaRatioEnteriter = linspace(1.5,5,cnt);  % iteration parameter ****
% Design.AreaEnteriter = linspace(0.01,5,cnt);  % iteration parameter ****

% Design.EnterDiaiter = linspace(0.5,5,cnt);  % iteration parameter ****
intake.DeflAngle = 15;   % Deflection angle (deg)

% fuel -----
fuel.Length = 15.00 /constants.In2Mtr;  % Grain Length (m)
fuel.Density = 1020;  % Grain Density (kg/m^3)
DesignStepHeight = 0.4;  % ratio

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
        for k=1:numel(Design.AreaRatioEnteriter) % Design.EnterDiaiter Design.AreaEnteriter Design.AreaRatioEnteriter
            n=1;
            vehicle.Mach(1)         = 2.5;  % Ramjet Initial Mach number
            trajectory.Z_pos(1)     = 1000.0;  % Initial altitude for ramjet start (m)
            
            nozzle.DiaThroatINCH = Design.DiaThroatINCHiter(i);
            nozzle.AreaRatioExit = Design.AreaRatioExititer(j);
%             intake.Area_enter = Design.AreaEnteriter(k);
%             intake.EnterDiaINCH = Design.EnterDiaiter(k);  % <in>

            

            
            nozzle.DiaThroat = nozzle.DiaThroatINCH/constants.In2Mtr;  % <m>
            nozzle.Area_throat = pi*(nozzle.DiaThroat)^2/4;  % Throat area (m^2)
            nozzle.Area_exit = nozzle.Area_throat*nozzle.AreaRatioExit;  % nozzle exit area (m^2)
            nozzle.DiaExit = sqrt(nozzle.Area_exit*4/pi);  % nozzle exit area (m)
            nozzle.DiaExitINCH = nozzle.DiaExit*constants.In2Mtr;  % nozzle exit area (in)
            
            intake.Area_enter = nozzle.Area_throat/Design.AreaRatioEnteriter(k);  % <m^2>
            
            
            trajectory.Rho_a(1)     = interp1(GRAM.Hgtkm, GRAM.DensMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Density (kg/m^3)
            trajectory.pressure_a(1)= interp1(GRAM.Hgtkm, GRAM.PresMean, (trajectory.Z_pos(1))/1e3);                % Atmospheric Pressure (Pa)
            trajectory.Temp_a(1)    = interp1(GRAM.Hgtkm, GRAM.Tmean, (trajectory.Z_pos(1))/1e3);                   % Atmospheric Temperature (K)


%             Initialize  % set initial conditions
            % Initialize BurnTime variable (s)
            BurnTime(1)         = 0.0;
            StopBurn = 0;
            while StopBurn == 0
                if n > 3  % stop simulation at 3rd iteration
                    % add variables to track

                    % store data if 3 iterations ran succsesfully
                    Design.intakemassFlow(i,j,k) = intake.massFlow(2,3);
                    Design.FreestreamStag(i,j,k) = intake.stagPres(1,3);
                    Design.MaxAvailableStag(i,j,k) = intake.MaxAvailableStag(3);
                    Design.chokeStagPres(i,j,k) = intake.chokeStagPres(3);
                    Design.combustionStagPres(i,j,k) = combustion.stagPres(2,3);
                    Design.Thrust(i,j,k) = vehicle.thrust(3);
                    Design.TSFC(i,j,k) = vehicle.TSFC(3);
                    Design.Specific_Thrust(i,j,k) = vehicle.Specific_Thrust(3);
                    Design.Phi(i,j,k) = combustion.Phi(3);
                    Design.AFT(i,j,k) = vehicle.AFT(3);
                    Design.FuelMassFlow(i,j,k) = fuel.MassFlow(3);
                    Design.AmbientPressure(i,j,k) = trajectory.pressure_a(1);
                    Design.NozzleExitPressure(i,j,k) = nozzle.staticPres(2,3);
                    Design.NozzleExitMach(i,j,k) = nozzle.mach(2,3);
                    Design.NozzleExitVelocity(i,j,k) = nozzle.velocity(2,3);
                    Design.StaticTempCombustion(i,j,k) = intake.staticTemp(5,3);
                    
                    
                    
                    % flag data good
                    Design.intakemassFlowgood(i,j,k) = 1;
                    Design.FreestreamStaggood(i,j,k) = 1;
                    Design.MaxAvailableStaggood(i,j,k) = 1;
                    Design.chokeStagPresgood(i,j,k) = 1;
                    Design.combustionStagPresgood(i,j,k) = 1;
                    Design.Thrustgood(i,j,k) = 1;
                    Design.TSFCgood(i,j,k) = 1;
                    Design.Specific_Thrustgood(i,j,k) = 1;
                    Design.Phigood(i,j,k) = 1;
                    Design.AFTgood(i,j,k) = 1;
                    Design.FuelMassFlowgood(i,j,k) = 1;
                    Design.AmbientPressuregood(i,j,k) = 1;
                    Design.NozzleExitPressuregood(i,j,k) = 1;
                    Design.NozzleExitMachgood(i,j,k) = 1;
                    Design.NozzleExitVelocitygood(i,j,k) = 1;
                    break
                end
                BurnTime(n) = time;                                     % Simulation Time
                if Burnout == 0
                    try
%                         RegressionRate                                      % Call Regression Rate Model
                        IntakeDesign                                              % Call Intake Model
                        % heights in <in> for plotting
                        Design.EnterHeight(i,j,k) = intake.EnterDiaINCH;  % <in>
                        Design.CowlHeight(i,j,k) = intake.CowlDiaINCH;  % <in>
                        Design.ThroatHeight(i,j,k) = nozzle.DiaThroatINCH;  % <in>
                        Design.ExitHeight(i,j,k) = nozzle.DiaExitINCH;  % <in>

                        % areas in <in^2> for plotting
                        Design.IntakeArea(i,j,k) = intake.Area_enter*constants.In2Mtr^2;  %  <in^2>
                        Design.ThroatArea(i,j,k) = nozzle.Area_throat*constants.In2Mtr^2;  %  <in^2>
                        Design.ExitArea(i,j,k) = nozzle.Area_exit*constants.In2Mtr^2;  %  <in^2>



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
                        
                        
                    % try to recover data from failed run for analysis
                    catch ME
%                         fprintf('something went wrong:( %s\n', ME.message);
                        if preserve_data == 1
                            try
                                Design.intakemassFlow(i,j,k) = intake.massFlow(2,n);
                            catch
                                Design.intakemassFlow(i,j,k) = NaN;
                            end
                            try
                                Design.FreestreamStag(i,j,k) = intake.stagPres(1,n);
                            catch
                                Design.FreestreamStag(i,j,k) = NaN;
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
                            try
                                Design.Thrust(i,j,k) = vehicle.thrust(n);
                            catch
                                Design.Thrust(i,j,k) = NaN;
                            end
                            try
                                Design.TSFC(i,j,k) = vehicle.TSFC(n);
                            catch
                                Design.TSFC(i,j,k) = NaN;
                            end
                            try
                                Design.Specific_Thrust(i,j,k) = vehicle.Specific_Thrust(n);
                            catch
                                Design.Specific_Thrust(i,j,k) = NaN;
                            end
                            %new
                            try
                                Design.FreestreamStag(i,j,k) = intake.stagPres(1,n);
                            catch
                                Design.FreestreamStag(i,j,k) = NaN;
                            end
                            try
                                Design.Phi(i,j,k) = combustion.Phi(n);
                            catch
                                Design.Phi(i,j,k) = NaN;
                            end
                            try
                                Design.AFT(i,j,k) = vehicle.AFT(n);
                            catch
                                Design.AFT(i,j,k) = NaN;
                            end
                            try
                                Design.FuelMassFlow(i,j,k) = fuel.MassFlow(n);
                            catch
                                Design.FuelMassFlow(i,j,k) = NaN;
                            end
                            try
                                Design.AmbientPressure(i,j,k) = trajectory.pressure_a(1);
                            catch
                                Design.AmbientPressure(i,j,k) = NaN;
                            end
                            try
                                Design.NozzleExitPressure(i,j,k) = nozzle.staticPres(2,n);
                            catch
                                Design.NozzleExitPressure(i,j,k) = NaN;
                            end
                            try
                                Design.NozzleExitMach(i,j,k) = nozzle.mach(2,n);
                            catch
                                Design.NozzleExitMach(i,j,k) = NaN;
                            end                            
                            try
                                Design.NozzleExitVelocity(i,j,k) = nozzle.velocity(2,n);
                            catch
                                Design.NozzleExitVelocity(i,j,k) = NaN;
                            end                            
                                                
                    
                            

                        else
                            Design.intakemassFlow(i,j,k) = NaN;
                            Design.FreestreamStag(i,j,k) = NaN;
                            Design.MaxAvailableStag(i,j,k) = NaN;
                            Design.chokeStagPres(i,j,k) = NaN;
                            Design.combustionStagPres(i,j,k) = NaN;
                            Design.Thrust(i,j,k) = NaN;
                            Design.TSFC(i,j,k) = NaN;
                            Design.Specific_Thrust(i,j,k) = NaN;
                            Design.Phi(i,j,k) = NaN;
                            Design.AFT(i,j,k) = NaN;
                            Design.FuelMassFlow(i,j,k) = NaN;
                            Design.AmbientPressure(i,j,k) = trajectory.pressure_a(1);
                            Design.NozzleExitPressure(i,j,k) = NaN;
                            Design.NozzleExitMach(i,j,k) = NaN;
                            Design.NozzleExitVelocity(i,j,k) = NaN;
                        end
                        % flag data bad
                        Design.intakemassFlowgood(i,j,k) = 0;
                        Design.FreestreamStaggood(i,j,k) = 0;
                        Design.MaxAvailableStaggood(i,j,k) = 0;
                        Design.chokeStagPresgood(i,j,k) = 0;
                        Design.combustionStagPresgood(i,j,k) = 0;
                        Design.Thrustgood(i,j,k) = 0;
                        Design.TSFCgood(i,j,k) = 0;
                        Design.Specific_Thrustgood(i,j,k) = 0;
                        Design.Phigood(i,j,k) = 0;
                        Design.AFTgood(i,j,k) = 0;
                        Design.FuelMassFlowgood(i,j,k) = 0;
                        Design.AmbientPressuregood(i,j,k) = trajectory.pressure_a(1);
                        Design.NozzleExitPressuregood(i,j,k) = 0;
                        Design.NozzleExitMachgood(i,j,k) = 0;
                        Design.NozzleExitVelocitygood(i,j,k) = 0;
                        break
                    end
                end
                time = time + SFRJDt;  % Step through simulation time
                n = n + 1; % Increase Index
            end
            
        end
    end
end
clear gamma R

% plot trade study results
%% individual studies-------

% Entrance
group = 'Entrance Opening vs. ';
xname = 'Height <in>';
indx = floor(cnt/2);
X = squeeze(Design.CowlHeight(indx,indx,:)) - squeeze(Design.EnterHeight(indx,indx,:));
figure()
plot(X,squeeze(Design.Thrust(indx,indx,:)));
xlabel(xname)
ylabel('Thrust <N>')
title(append(group, 'thrust'))
figure()
plot(X,squeeze(Design.TSFC(indx,indx,:)));
xlabel(xname)
ylabel('Thrust specific fuel consumption <Kg/(N*s)>')
title(append(group, 'TSFC'))
figure()
plot(X,squeeze(Design.intakemassFlow(indx,indx,:)));
xlabel(xname)
ylabel('Massflow <kg/s>')
title(append(group, 'Massflow'))
figure()
hold on 
plot(X,squeeze(Design.chokeStagPres(indx,indx,:)),'*');
plot(X,squeeze(Design.combustionStagPres(indx,indx,:)));
plot(X,squeeze(Design.MaxAvailableStag(indx,indx,:)));
plot(X,squeeze(Design.FreestreamStag(indx,indx,:)));
hold off
legend('pressure required to choke','chamber pressure','available pressure (after shocks)','freestream pressure')
xlabel(xname)
ylabel('Pressure (stag) <Pa>')
title(append(group, 'Pressures'))
figure()
plot(X,squeeze(Design.CowlHeight(indx,indx,:)));
xlabel(xname)
ylabel('Outer Diameter  <in>')
title(append(group, 'Ramjet Outer Diameter'))


% Throat
group = 'Nozzle Throat Diameter vs. ';
xname = 'Diameter <in>';
X = squeeze(Design.ThroatHeight(:,indx,indx));
figure()
plot(X,squeeze(Design.Thrust(:,indx,indx)));
xlabel(xname)
ylabel('Thrust <N>')
title(append(group, 'thrust'))
figure()
plot(X,squeeze(Design.TSFC(:,indx,indx)));
xlabel(xname)
ylabel('Thrust specific fuel consumption <Kg/(N*s)>')
title(append(group, 'TSFC'))
figure()
plot(X,squeeze(Design.intakemassFlow(:,indx,indx)));
xlabel(xname)
ylabel('Massflow <kg/s>')
title(append(group, 'Massflow'))
figure()
hold on 
plot(X,squeeze(Design.chokeStagPres(:,indx,indx)),'*');
plot(X,squeeze(Design.combustionStagPres(:,indx,indx)));
plot(X,squeeze(Design.MaxAvailableStag(:,indx,indx)));
plot(X,squeeze(Design.FreestreamStag(:,indx,indx)));
hold off
legend('pressure required to choke','chamber pressure','available pressure (after shocks)','freestream pressure')
xlabel(xname)
ylabel('Pressure (stag) <Pa>')
title(append(group, 'Pressures'))
figure()
plot(X,squeeze(Design.CowlHeight(:,indx,indx)));
xlabel(xname)
ylabel('Outer Diameter  <in>')
title(append(group, 'Ramjet Outer Diameter'))
figure()
plot(X,(squeeze(Design.CowlHeight(:,indx,indx)) - squeeze(Design.EnterHeight(:,indx,indx))));
xlabel(xname)
ylabel('Opening height  <in>')
title(append(group, 'Entrance Opening'))
figure()
plot(X,squeeze(Design.IntakeArea(:,indx,indx)));
xlabel(xname)
ylabel('entrance area  <in^2>')
title(append(group, 'Entrance Area'))
% freestyle
figure()
plot(X,squeeze(Design.IntakeArea(:,indx,indx)),X,squeeze(Design.ThroatArea(:,indx,indx)),X,squeeze(Design.ExitArea(:,indx,indx)));
xlabel(xname)
ylabel('area <in^2>')
legend('intake area','throat area','exit area')
title(append(group, 'Exit Velocity'))


% Exit
group = 'Exit Diameter vs. ';
xname = 'Diameter <in>';
X = squeeze(Design.ExitHeight(indx,:,indx));
figure()
plot(X,squeeze(Design.Thrust(indx,:,indx)));
xlabel(xname)
ylabel('Thrust <N>')
title(append(group, 'thrust'))
figure()
plot(X,squeeze(Design.TSFC(indx,:,indx)));
xlabel(xname)
ylabel('Thrust specific fuel consumption <Kg/(N*s)>')
title(append(group, 'TSFC'))
figure()
plot(X,squeeze(Design.NozzleExitMach(indx,:,indx)));
xlabel(xname)
ylabel('mach number')
title(append(group, 'Exit Mach'))
figure()
plot(X,squeeze(Design.NozzleExitVelocity(indx,:,indx)));
xlabel(xname)
ylabel('velocity <m/s>')
title(append(group, 'Exit Velocity'))
figure()
hold on 
plot(X,squeeze(Design.AmbientPressure(indx,:,indx)),'*');
plot(X,squeeze(Design.NozzleExitPressure(indx,:,indx)));
hold off
legend('ambient pressure', 'exit pressure')
xlabel(xname)
ylabel('Pressure (static) <Pa>')
title(append(group, 'Pressures'))

%% 2 dimension trade studies-------

% -- entrance vs nozzle throat
indx = floor(cnt/2);
entrance_opening = squeeze(Design.CowlHeight(indx,indx,:)) - squeeze(Design.EnterHeight(indx,indx,:));
[X,Y] = ndgrid(squeeze(Design.ThroatHeight(:,indx,indx)),entrance_opening);  % Design.IntakeArea
xname = 'throat Diameter';
yname = 'Entrance Opening';
xunit = ' <in>';
yunit = ' <in>';
group = append(xname,' vs ',yname);
% thrust
figure()
surf(X,Y,squeeze(Design.Thrust(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('thrust <N>')
title(append(group,' vs Thrust'))
% TSFC
figure()
surf(X,Y,squeeze(Design.TSFC(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('TSFC <Kg/(N*s)>')
title(append(group,' vs TSFC'))
% specific thrust
figure()
surf(X,Y,squeeze(Design.Specific_Thrust(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('specific thrust <(N*s)/Kg>')
title(append(group,' vs Specific Thrust'))
% air massflow
figure()
surf(X,Y,squeeze(Design.intakemassFlow(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('massflow <Kg/s>')
title(append(group,' vs Air Massflow'))
% fuel massflow
figure()
surf(X,Y,squeeze(Design.FuelMassFlow(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('massflow <Kg/s>')
title(append(group,' vs Fuel Massflow'))
% pressure
figure()
surf(X,Y,squeeze(Design.MaxAvailableStag(:,indx,:)),'FaceAlpha',0.2,'FaceColor','r');
hold on
surf(X,Y,squeeze(Design.chokeStagPres(:,indx,:)),'FaceColor','b');
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('pressure (stag) <Pa>')
title(append(group,' vs Pressures'))
legend('available stagnation (after shock)', 'choke stagnation required')
hold off
% Equivalence Ratio
figure()
surf(X,Y,squeeze(Design.Phi(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('\phi')
title(append(group,' vs Equivalence Ratio'))
% AFT
figure()
surf(X,Y,squeeze(Design.AFT(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('Temperature <K>')
title(append(group,' vs Adiabatic Flame Temperature'))
% chamber size
figure()
surf(X,Y,squeeze(Design.CowlHeight(:,indx,:)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('Outer Diameter <in>')
title(append(group,' vs Ramjet Outer Diameter'))



% -- nozzle throat vs Exit Diameter
[X,Y] = ndgrid(squeeze(Design.ThroatHeight(:,indx,indx)),squeeze(Design.ExitHeight(indx,:,indx)));
xname = 'throat Diameter';
yname = 'exit Diameter';
xunit = ' <in>';
yunit = ' <in>';
group = append(xname,' vs ',yname);
% thrust
figure()
surf(X,Y,squeeze(Design.Thrust(:,:,indx)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('thrust <N>')
title(append(group,' vs Thrust'))
% TSFC
figure()
surf(X,Y,squeeze(Design.TSFC(:,:,indx)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('TSFC <Kg/(N*s)>')
title(append(group,' vs TSFC'))
% specific thrust
figure()
surf(X,Y,squeeze(Design.Specific_Thrust(:,:,indx)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('specific thrust <(N*s)/Kg>')
title(append(group,' vs Specific Thrust'))
% pressure
figure()
surf(X,Y,squeeze(Design.AmbientPressure(:,:,indx)),'FaceAlpha',0.2,'FaceColor','r');
hold on
surf(X,Y,squeeze(Design.NozzleExitPressure(:,:,indx)),'FaceColor','b');
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('pressure (static) <Pa>')
title(append(group,' vs Pressures'))
legend('ambient pressure', 'exit pressure')
hold off
% mach
figure()
surf(X,Y,squeeze(Design.NozzleExitMach(:,:,indx)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('mach number')
title(append(group,' vs Exit Mach'))
% velocity
figure()
surf(X,Y,squeeze(Design.NozzleExitVelocity(:,:,indx)));
xlabel(append(xname, xunit))
ylabel(append(yname, yunit))
zlabel('velocity <m/s>')
title(append(group,' vs Exit Velocity'))



                    

