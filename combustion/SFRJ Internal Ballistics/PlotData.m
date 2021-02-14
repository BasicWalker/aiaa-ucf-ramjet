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
% ---------------------------------------------------------------------- %
Thrustdlvd(n-1) = 0.0;
PC(n-1)         = 0.0;
OFRatio(n-1)    = 0.0;
InltPres(n-1)   = 0.0;
PCreq(n-1)      = 0.0;
PC_TAFT(n-1)    = 0.0;


% figure
% hold on
% plot(BurnTime,InltPres)
% plot(BurnTime,PC)
% plot(BurnTime,PCres)
% plot(BurnTime,BackPres)
% title('Pressure Plots')
% xlabel('Time (s)')
% ylabel('Pressure (kPa)')
% grid on
% legend('Inlet Pressure','Chamber Pressure','Required Pressure','Back Pressure')
% ylim([0 900])

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

figure('Name','Thrust Curve')
plot(BurnTime,Thrustdlvd)
title('Thrust Curve')
xlabel('Time (s)')
ylabel('Thrust (N)')
grid on
% ylim([0 4000])

% figure('Name','Instantaneous Mass Properties')
% plot(BurnTime,MassGen, BurnTime, MairGen, BurnTime, MFuelGen)
% title('Instantaneous Mass Properties')
% xlabel('Time (s)')
% ylabel('Mass (kg)')
% legend('Total Mass','Air Mass','Fuel Mass')
% grid on

figure('Name','Mass Flow Rate')
plot(BurnTime,MassFlow, BurnTime, MdotAir, BurnTime, MdotFuel)
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

% figure('Name','Pressure Plots')
% plot(BurnTime, InltPres, BurnTime, PC, BurnTime, PCreq, BurnTime, BackPres, BurnTime, PC_TAFT)
% title('Pressure Plots')
% xlabel('Time (s)')
% ylabel('Pressure (kPa)')
% legend('Inlet Pressure','Chamber Pressure - Legacy','Required Pressure - Choked Flow','Back Pressure', 'Chamber Pressure - Adiabatic Flame Temp')
% grid on

figure('Name','Pressure Plots')
plot(BurnTime, InltPres, BurnTime, PCreq, BurnTime, BackPres, BurnTime, PC_TAFT)
title('Pressure Plots')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
legend('Inlet Pressure','Required Pressure - Choked Flow','Back Pressure', 'Chamber Pressure - Adiabatic Flame Temp')
grid on


fprintf('\n------------ Simulation Results ------------\n')
fprintf('Burn Time:             %.2f   (s)\n', BurnTime(n-1))
fprintf('Average Thrust:        %.2f  (N)\n', mean(Thrustdlvd))
fprintf('Total Impulse:         %.2f  (Ns)\n', TotallImp(n-1))
fprintf('Inlet Velocity:        %.2f   (m/s)\n', InltVel(n-1))
fprintf('Inlet Mass Flow Rate:  %.4f  (kg/s)\n', InltMassFlw)
fprintf('PC From TAFT:          %.1f  (kPa)\n',mean(PC_TAFT))