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
    FuelMass(n) = 0.0;
    F_t(n) = 0.0;
end

if (n > 1) 
    thrust.F_tx(n-1)    = thrust.F_t(n-1)*cosd(alpha);                                                  % Thrust X (N)
    thrust.F_tz(n-1)    = thrust.F_t(n-1)*sind(alpha);                                                  % Thrust Z (N)
    trajectory.F_x(n-1) = thrust.F_tx(n-1)-trajectory.F_dx(n-1);                                        % Sum of Force X (N)
    trajectory.F_z(n-1) = thrust.F_tz(n-1)-trajectory.F_dz(n-1)-Mass(n-1)*gravity+trajectory.Lift(n-1); % Sum of Force Y (N)
    trajectory.Acc_x(n) = trajectory.F_x(n-1)/Mass(n-1);                                                % Acceleration X (m/s/s)
    trajectory.Acc_z(n) = trajectory.F_z(n-1)/Mass(n-1);                                                % Acceleration Z (m/s/s)
    trajectory.Vel_x(n) = trajectory.Acc_x(n)*SFRJDt + trajectory.Vel_x(n-1);                           % Velocity X (m/s)
    trajectory.Vel_z(n) = trajectory.Acc_z(n)*SFRJDt + trajectory.Vel_z(n-1);                           % Velocity Z (m/s)
    trajectory.Vel(n)   = sqrt(trajectory.Vel_x(n)^2 + trajectory.Vel_z(n)^2);                          % Total Velocity (m/s)
    trajectory.Vel_t(n) = sqrt((2*Mass(n-1)*gravity)/(c_d*trajectory.Rho_a(n-1)*S));                    % Terminal Velocity (m/s)
    trajectory.X_pos(n) = 0.5*trajectory.Acc_x(n)*SFRJDt^2 + trajectory.Vel_x(n-1)*SFRJDt + trajectory.X_pos(n-1);     % Downrange Position of Vehicle (m)
    trajectory.Z_pos(n) = 0.5*trajectory.Acc_z(n)*SFRJDt^2 + trajectory.Vel_z(n-1)*SFRJDt + trajectory.Z_pos(n-1);     % Altitude Postion of Vehicle (m) 
    trajectory.Rho_a(n) = interp1(GRAM.Hgtkm, GRAM.DensMean, (trajectory.Z_pos(n))/1e3);                % Interpolated Lookup Table, atm density 
    trajectory.Temp_a(n)= interp1(GRAM.Hgtkm, GRAM.Tmean, (trajectory.Z_pos(n))/1e3);                   % Interpolated Lookup Table, atm temp
    trajectory.F_d(n)   = 0.5*trajectory.Rho_a(n)*c_d*S*trajectory.Vel(n)^2;                            % Drag (N)
    trajectory.F_dx(n)  = abs(0.5*trajectory.Rho_a(n)*c_d*S*trajectory.Vel(n)*trajectory.Vel_x(n));     % Drag X (N)
    trajectory.F_dz(n)  = abs(0.5*trajectory.Rho_a(n)*c_d*S*trajectory.Vel(n)*trajectory.Vel_z(n));     % Drag Z (N)
    Mach_f(n) = trajectory.Vel(n)/sqrt(k*R*trajectory.Temp_a(n));                                       % Vehicle Mach
    Mass(n)   = dry_mass + FuelMass(n);                                                                 % Vehicle Mass
    Weight(n) = gravity*Mass(n);                                                                        % Vehicle weight
    trajectory.Lift(n)   = Mass(n)*gravity*trajectory.LiftOnOff;                                        % Lift (N) - equal to weight of vehicle
    trajectory.pressure_a(n) = interp1(GRAM.Hgtkm, GRAM.PresMean, (trajectory.Z_pos(n))/1e3);           % Interpolated Lookup Table, atm pressure
    trajectory.pressure_a(n) = trajectory.pressure_a(n)/Pa2kPa;                                         % Convert to Kpa
    
    if trajectory.F_x(n-1) < 0.0 || trajectory.F_z(n-1) < 0.0
        trajectory.F_net(n-1)= -1*sqrt(trajectory.F_x(n-1)^2 + trajectory.F_z(n-1)^2);                  % Net Force (N)
        trajectory.Acc(n)    = -1* trajectory.F_net(n-1)/Mass(n-1);                                     % Acceleration (m/s/s)
    else
        trajectory.F_net(n-1)= sqrt(trajectory.F_x(n-1)^2 + trajectory.F_z(n-1)^2);                     % Net Force (N)
        trajectory.Acc(n)    = trajectory.F_net(n-1)/Mass(n-1);                                         % Acceleration (m/s/s)
    end
    
    if trajectory.Z_pos(n) <= 0.0
        StopBurn = true;                                                                                % Stop simulation
    end
end

