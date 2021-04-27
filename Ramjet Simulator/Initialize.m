% Initialize BurnTime variable (s)
BurnTime(1)         = 0.0; 

% Initialize Fuel Grain Model
fuel.PortArea(1)         = pi*(fuel.DiaInner(1)^2)*(1/4);  % Fuel Port Area (m^2)
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