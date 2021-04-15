% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: PlotData.m 
%
% File Description: 
% Plots all results
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% Ethan Sherlock  02/14/21  001  Updated plots for SCR 001
% Ethan Sherlock  02/14/21  005  1DOF model update
% Ethan Sherlock  03/12/21  002  Chem model update
% Ethan Sherlock  04/14/21  ---  2DOF model update
% ---------------------------------------------------------------------- %
thrust.F_t(n-1)         = trajectory.F_t(n-2);
PC(n-1)                 = 0.0;
OFRatio(n-1)            = 0.0;
InltPres(n-1)           = 0.0;
PCreq(n-1)              = 0.0;
PC_TAFT(n-1)            = 0.0;
trajectory.F_d(n-1)     = trajectory.F_d(n-2);
trajectory.Vel(n-1)     = 0.0;
trajectory.Acc(n-1)     = 0.0;
trajectory.F_net(n-1)   = trajectory.F_net(n-2);
MassFlow(n-1)           = 0.0;
MdotAir(n-1)            = 0.0;
MdotFuel(n-1)           = 0.0;
InltPres_stag(n-1)      = 0.0;
PCreq(n-1)              = 0.0;
PC_TAFT(n-1)            = 0.0;
phi_eqv(n-1)            = 0.0;
T_stag(n-1)             = 0.0;
gamma_t(n-1)            = 0.0;
TotallImp(n-1)          = 0.0;

fprintf('Simulation Complete \n')

figure('Name','O/F Ratio')
plot(BurnTime(1:index-1),OFRatio(1:index-1))
title('O/F Ratio')
xlabel('Time (s)')
ylabel('')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Force Profile')
plot(BurnTime,thrust.F_t,BurnTime,trajectory.F_d, BurnTime, Mass*gravity, BurnTime, trajectory.F_net)
title('Thrust & Drag vs Time')
xlabel('Time (s)')
ylabel('Force (N)')
legend('Thrust Curve','Drag','Weight','Net Force')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Mass Flow Rate')
plot(BurnTime(1:index-1),MassFlow(1:index-1), BurnTime(1:index-1), MdotAir(1:index-1), BurnTime(1:index-1), MdotFuel(1:index-1)) 
title('Mass Flow Rate')
xlabel('Time (s)')
ylabel('Mass (kg/s)')
legend('Total Mdot','Mdot Air','Mdot Fuel')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Pressure Plots')
plot(BurnTime(1:index-1), InltPres_stag(1:index-1), BurnTime(1:index-1), PCreq(1:index-1), BurnTime(1:index-1), pressure_a(1:index-1), BurnTime(1:index-1), PC_TAFT(1:index-1))
title('Pressure Plots')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
legend('Inlet Pressure (stag)','Required Pressure - Choked Flow (stag)','Back Pressure (stag)', 'Chamber Pressure - Adiabatic Flame Temp (stag)')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Flight Mach')
plot(BurnTime,Mach_f)
title('Vehicle Mach vs Time')
xlabel('Time (s)')
ylabel('Mach No.')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Equivalence Ratio')
plot(BurnTime(1:index-1),phi_eqv(1:index-1))
title('Equivalence Ratio vs Time')
xlabel('Time (s)')
ylabel('Equivalence Ratio (phi)')
grid on
% ylim([0 6])
set(gcf,'position',[550,200,800,700])

figure('Name','Adiabatic Flame Temperature')
yyaxis left
plot(BurnTime(1:index-1),T_stag(1:index-1))
title('Adiabatic Flame Temperature & Specific Heat Ratio')
xlabel('Time (s)')
ylabel('Temperature <k>')
yyaxis right
plot(BurnTime(1:index-1),gamma_t(1:index-1))
ylabel('Specific Heat Ratio')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Trajectory Profile')
axis equal
subplot(3,1,1)
plot(trajectory.X_pos,trajectory.Z_pos)
title('Trajectory')
xlabel('Range (m)')
ylabel('Altitude (m)')
grid on
subplot(3,1,2)
plot(BurnTime,trajectory.Z_pos)
title('Altitude vs Time')
xlabel('Time(s)')
ylabel('Altitude(m)')
grid on
subplot(3,1,3)
plot(BurnTime,trajectory.X_pos)
title('Range vs Time')
xlabel('Time(s)')
ylabel('Range(m)')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Velocity/Acceleration Profile')
axis equal
subplot(2,1,1)
grid on
title('Vertical Profile (Z)')
yyaxis left
plot(BurnTime, trajectory.Vel_z)
ylabel('Vz (m/s)')
yyaxis right
plot(BurnTime, trajectory.Acc_z)
ylabel('Az (m/s^2)')
subplot(2,1,2)
grid on;
title('Horizontal Profile (X)')
yyaxis left
plot(BurnTime, trajectory.Vel_x)
xlabel('Time (s)')
ylabel('Vx (m/s)')
yyaxis right
plot(BurnTime, trajectory.Acc_x)
ylabel('Ax (m/s^2)')
set(gcf,'position',[550,200,800,700])

if (mean(InltMach) > 0.2)
    fprintf(2,'\nWARNING: Inlet Mach number is too high.\n')
end
if (mean(phi_eqv) > 3.0)
    fprintf(2,'WARNING: Equivalence ratio exceeds flamability limit.\n')
end
if (Mach_f(n-1) < 2)
    fprintf(2,'WARNING: Vehicle does not reach target Mach Number.\n')
end

fprintf('------------ Simulation Results ------------\n')
fprintf('Burn Time:                  %.2f    (s)\n', BurnTime(n-1))
fprintf('Average Thrust:             %.2f   (N)\n', mean(trajectory.F_t))
fprintf('Average Drag Force:         %.2f   (N)\n',mean(trajectory.F_d))
fprintf('Total Impulse:              %.2f  (Ns)\n', max(TotallImp))
fprintf('Air Mass Flow Rate:         %.3f    (kg/s)\n', mean(MdotAir))
fprintf('Chamber Pressure From TAFT: %.2f   (kPa)\n',mean(PC_TAFT))
fprintf('Initial Step Height:        %.2f     (in) \n', StepHeight(1)*In2Mtr)
fprintf('Average Inlet Velocity:     %.2f   (m/s)\n', mean(InltVel))
fprintf('Max Altitude(During Boost): %.2f (m)\n',max(trajectory.Z_pos))
fprintf('Average Inlet Mach:         %.2f  \n', mean(InltMach))
fprintf('Average O/F Ratio:          %.2f  \n',mean(OFRatio))
fprintf('Average Equivalence Ratio:  %.2f  \n',mean(phi_eqv))
fprintf('Average Gamma Ratio:        %.2f  \n',mean(gamma_t))
fprintf('--------------------------------------------\n')
toc