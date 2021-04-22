% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Trajectory.m 
% 
% File Description: 
% Trajectory model for the SFRJ
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  005  1DOF trajectory update
% Ethan Sherlock  04/03/21  ---  Update to 2DOF
% ---------------------------------------------------------------------- %
if Burnout == 1
    fuel.Mass(n) = 0.0;
    F_t(n) = 0.0;
end


% resolve kinematic components
trajectory.F_tx(n)    = vehicle.thrust(n)*cosd(alpha);                                                  % Thrust X (N)
trajectory.F_tz(n)    = vehicle.thrust(n)*sind(alpha);                                                  % Thrust Z (N)
trajectory.F_x(n) = trajectory.F_tx(n)-trajectory.F_dx(n);                                        % Sum of Force X (N)
trajectory.F_z(n) = trajectory.F_tz(n)-trajectory.F_dz(n)-vehicle.TotalMass(n)*constants.gravity+trajectory.Lift(n);     % Sum of Force Y (N)
trajectory.Vel_mag(n)   = sqrt(abs(trajectory.Vel_x(n)^2 + trajectory.Vel_z(n)^2));                   % Velocity magnitude (m/s)
trajectory.F_mag(n) = sqrt(abs(trajectory.F_x(n)^2 + trajectory.F_z(n)^2));                       % Force magnitude (N)
trajectory.Acc_mag(n)    = trajectory.F_mag(n)/vehicle.TotalMass(n);                                         % Acceleration magnitude (m/s/s)

% update trajectory
trajectory.Vel_t(n+1) = sqrt((2*vehicle.TotalMass(n)*constants.gravity)/(vehicle.DragCoeff*trajectory.Rho_a(n)*vehicle.FrontSurfArea));                    % Terminal Velocity (m/s)
trajectory.Vel_x(n+1) = trajectory.Acc_x(n)*SFRJDt + trajectory.Vel_x(n);                           % Velocity X (m/s)
trajectory.Vel_z(n+1) = trajectory.Acc_z(n)*SFRJDt + trajectory.Vel_z(n);                           % Velocity Z (m/s)
trajectory.Acc_x(n+1) = trajectory.F_x(n)/vehicle.TotalMass(n);                                                % Acceleration X (m/s/s)
trajectory.Acc_z(n+1) = trajectory.F_z(n)/vehicle.TotalMass(n);                                                % Acceleration Z (m/s/s)
trajectory.X_pos(n+1) = 0.5*trajectory.Acc_x(n)*SFRJDt^2 + trajectory.Vel_x(n)*SFRJDt + trajectory.X_pos(n);     % Downrange Position of Vehicle (m)
trajectory.Z_pos(n+1) = 0.5*trajectory.Acc_z(n)*SFRJDt^2 + trajectory.Vel_z(n)*SFRJDt + trajectory.Z_pos(n);     % Altitude Postion of Vehicle (m) 
trajectory.Rho_a(n+1) = interp1(GRAM.Hgtkm, GRAM.DensMean, (trajectory.Z_pos(n))/1e3);                % Interpolated Lookup Table, atm density 
trajectory.Temp_a(n+1)= interp1(GRAM.Hgtkm, GRAM.Tmean, (trajectory.Z_pos(n))/1e3);                   % Interpolated Lookup Table, atm temp
trajectory.pressure_a(n+1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (trajectory.Z_pos(n))/1e3);           % Interpolated Lookup Table, atm pressure 
trajectory.F_d(n+1)   = 0.5*trajectory.Rho_a(n)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel_mag(n)^2;                            % Drag (N)
trajectory.F_dx(n+1)  = abs(0.5*trajectory.Rho_a(n)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel_x(n)^2);     % Drag X (N)
trajectory.F_dz(n+1)  = abs(0.5*trajectory.Rho_a(n)*vehicle.DragCoeff*vehicle.FrontSurfArea*trajectory.Vel_z(n)^2);     % Drag Z (N)
vehicle.Mach(n+1) = trajectory.Vel_mag(n)/sqrt(constants.gamma*constants.R*trajectory.Temp_a(n));                                       % Vehicle Mach
vehicle.TotalMass(n+1)   = vehicle.DryMass + fuel.Mass(n);                                                                 % Vehicle Mass
vehicle.Weight(n+1) = constants.gravity*vehicle.TotalMass(n);                                                                        % Vehicle weight
trajectory.Lift(n+1)   = vehicle.TotalMass(n)*constants.gravity*trajectory.LiftOnOff;                                        % Lift (N) - equal to weight of vehicle

% Stop simulation *crash*
if trajectory.Z_pos(n) <= 0.0
    StopBurn = true; 
    warning('trajectory:crash',...
        'Ground impact; altitude <= 0\nsimulation stopped prematurely \n')
end


