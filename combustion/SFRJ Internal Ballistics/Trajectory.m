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
    %fprintf('Running... Motor Burnout Detected\n')               % Running Simulator indicator
end

if (n > 1) 
    F_tx(n-1) = F_t(n-1)*cosd(alpha);                                       % Thrust X (N)
    F_tz(n-1) = F_t(n-1)*sind(alpha);                                       % Thrust Z (N)
    F_x(n-1)  = F_tx(n-1) - F_dx(n-1);                                      % Sum of Force X (N)
    F_z(n-1)  = F_tz(n-1) - F_dz(n-1) - Mass(n-1)*gravity;                  % Sum of Force Y (N)
    Acc_x(n)  = F_x(n-1)/Mass(n-1);                                         % Acceleration X (m/s/s)
    Acc_z(n)  = F_z(n-1)/Mass(n-1);                                         % Acceleration Z (m/s/s)
    Vel_x(n)  = Acc_x(n)*SFRJDt + Vel_x(n-1);                               % Velocity X (m/s)
    Vel_z(n)  = Acc_z(n)*SFRJDt + Vel_z(n-1);                               % Velocity Z (m/s)
    Vel(n)    = sqrt(Vel_x(n)^2 + Vel_z(n)^2);                              % Total Velocity (m/s)
    Vel_t(n)  = sqrt((2*Mass(n-1)*gravity)/(c_d*Rho_a(n-1)*S));             % Terminal Velocity (m/s)
    X_pos(n)  = 0.5*Acc_x(n)*SFRJDt^2 + Vel_x(n-1)*SFRJDt + X_pos(n-1);     % Downrange Position of Vehicle (m)
    Z_pos(n)  = 0.5*Acc_z(n)*SFRJDt^2 + Vel_z(n-1)*SFRJDt + Z_pos(n-1);     % Altitude Postion of Vehicle (m) 
    Rho_a(n)  = interp1(GRAM.Hgtkm, GRAM.DensMean, (Z_pos(n))/1e3);         % Interpolated Lookup Table, atm density 
    Temp_a(n) = interp1(GRAM.Hgtkm, GRAM.Tmean, (Z_pos(n))/1e3);            % Interpolated Lookup Table, atm temp
    F_d(n)    = 0.5*Rho_a(n)*c_d*S*Vel(n)^2;                                % Drag (N)
    F_dx(n)   = abs(0.5*Rho_a(n)*c_d*S*Vel(n)*Vel_x(n));                    % Drag X (N)
    F_dz(n)   = abs(0.5*Rho_a(n)*c_d*S*Vel(n)*Vel_z(n));                    % Drag Z (N)
    Mach_f(n) = Vel(n)/sqrt(k*R*Temp_a(n));                                 % Vehicle Mach
    Mass(n)   = dry_mass + FuelMass(n);                                     % Vehicle Mass
    Weight(n) = gravity*Mass(n);                                            % Vehicle weight
    pressure_atm(n) = interp1(GRAM.Hgtkm, GRAM.PresMean, (Z_pos(n))/1e3);   % Interpolated Lookup Table, atm pressure
    pressure_atm(n) = pressure_atm(n)/Pa2kPa;                               % Convert to Kpa
    
%     if Burnout == 1
%        if Vel_z(n) < -1*Vel_t(n)
%           Vel_z(n) = -1*Vel_t(n);
%        end
%     end
    
    if F_x(n-1) < 0.0 || F_z(n-1) < 0.0
        F_net(n-1)= -1*sqrt(F_x(n-1)^2 + F_z(n-1)^2);                       % Net Force (N)
        Acc(n)    = -1* F_net(n-1)/Mass(n-1);                               % Acceleration (m/s/s)
    else
        F_net(n-1)= sqrt(F_x(n-1)^2 + F_z(n-1)^2);                          % Net Force (N)
        Acc(n)    = F_net(n-1)/Mass(n-1);                                   % Acceleration (m/s/s)
    end
    
    if Z_pos(n) <= 0.0
        StopBurn = true;
    end
    
end

