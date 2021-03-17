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
% Ethan Sherlock  02/14/21  005  1DOF trajectory update
% ---------------------------------------------------------------------- %
Thrustdlvd(n-1)     = 0.0;
PC(n-1)             = 0.0;
OFRatio(n-1)        = 0.0;
InltPres(n-1)       = 0.0;
PCreq(n-1)          = 0.0;
PC_TAFT(n-1)        = 0.0;
drag(n-1)           = 0.0;
velocity(n-1)       = 0.0;
acceleration(n-1)   = 0.0;

figure('Name','O/F Ratio')
plot(BurnTime,OFRatio)
title('O/F Ratio')
xlabel('Time (s)')
ylabel('')
grid on
% ylim([0 12])

% figure('Name','A/F Ratio')
% hold on
% plot(BurnTime,AFRatio, BurnTime,AFRst)
% title('Air to Fuel Ratios')
% xlabel('Time (s)')
% ylabel('')
% legend('A/F Results','A/F Target')
% grid on

figure('Name','Force Profile')
plot(BurnTime,Thrustdlvd,BurnTime,drag)
title('Thrust & Drag vs Time')
xlabel('Time (s)')
ylabel('Force (N)')
legend('Thrust Curve','Drag')
grid on

figure('Name','Mass Flow Rate')
plot(BurnTime,MassFlow, BurnTime, MdotAir, BurnTime, MdotFuel) %
title('Mass Flow Rate')
xlabel('Time (s)')
ylabel('Mass (kg/s)')
legend('Total Mdot','Mdot Air','Mdot Fuel')
grid on

% figure('Name','Step Height')
% plot(BurnTime,StepHeight)
% title('Step Height')
% xlabel('Time (s)')
% ylabel('Height (m)')
% grid on
% ylim([0 .05])

figure('Name','Pressure Plots')
plot(BurnTime, InltPres_stag, BurnTime, PCreq, BurnTime, pressure_atm, BurnTime, PC_TAFT)
title('Pressure Plots')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
legend('Inlet Pressure','Required Pressure - Choked Flow','Back Pressure', 'Chamber Pressure - Adiabatic Flame Temp')
grid on

figure('Name','Velocity & Acceleration Profiles')
yyaxis left
plot(BurnTime, velocity)
title('Velocity and Acceleration Profiles')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
yyaxis right
plot(BurnTime, acceleration)
ylabel('Acceleration (m/s^2)')
legend('Velocity','Acceleration')
grid on

figure('Name','Trajectory Profile')
plot(BurnTime,altitude)
title('Altitude vs Time')
xlabel('Time (s)')
ylabel('Altitude (m)')
grid on
ylim([0 max(altitude)+500])

figure('Name','Flight Mach')
plot(BurnTime,flight_mach)
title('Vehicle Mach vs Time')
xlabel('Time (s)')
ylabel('Mach No.')
grid on

figure('Name','Equivalence Ratio')
plot(BurnTime,phi_eqv)
title('Equivalence Ratio vs Time')
xlabel('Time (s)')
ylabel('Equivalence Ratio (phi)')
grid on
ylim([0 11])

figure('Name','Adiabatic Flame Temperature')
plot(BurnTime,T_stag)
title('Adiabatic Flame Temperature vs Time')
xlabel('Time (s)')
ylabel('Temperature <k>')
grid on
% ylim([0 8])

fprintf('\n------------ Simulation Results ------------\n')
fprintf('Burn Time:                 %.2f   (s)\n', BurnTime(n-1))
fprintf('Average Thrust:            %.2f  (N)\n', mean(Thrustdlvd))
fprintf('Average Drag Force:        %.2f  (N)\n',mean(drag))
fprintf('Total Impulse:             %.2f (Ns)\n', TotallImp(n-1))
fprintf('Inlet Mass Flow Rate:      %.3f   (kg/s)\n', mean(InltMassFlw))
fprintf('PC From TAFT:              %.2f  (kPa)\n',mean(PC_TAFT))
fprintf('Initial Step Height:       %.2f    (in) \n', StepHeight(1)*In2Mtr)
fprintf('Average Inlet Velocity:    %.2f  (m/s)\n', mean(InltVel))
fprintf('Average Inlet Mach:        %.2f  \n', mean(InltMach))
fprintf('Average O/F Ratio:         %.2f      \n',mean(OFRatio))
fprintf('Average Equivalence Ratio: %.2f  \n',mean(phi))
fprintf('--------------------------------------------\n')


