% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: PlotData.m 
%
% File Description: 
% Plots all results
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %
Thrustdlvd(n-1) = 0.0;
PC(n-1)         = 0.0;
OFRatio(n-1)    = 0.0;
InltPres(n-1)   = 0.0;

figure
plot(BurnTime,FuelSA)
title('Fuel Grain Surface Area vs Time')
xlabel('Time (s)')
ylabel('Surface Area (m^2)')
grid on
% ylim([0 0.15])

figure
hold on
plot(BurnTime,PC)
plot(BurnTime,InltPres)
title('Combustion Chamber Pressure vs Inlet Pressure')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
grid on
legend('Chamber Pressure','Inlet Pressure')
ylim([0 900])

figure
plot(BurnTime,MFuelGen)
title('Fuel Mass Generated')
xlabel('Time (s)')
ylabel('Mass (kg)')
grid on

figure
plot(BurnTime,OFRatio)
title('O/F Ratio')
xlabel('Time (s)')
ylabel('')
grid on
% ylim([0 12])

figure
hold on
plot(BurnTime,AFRatio) 
plot(BurnTime,AFRst)
title('Air to Fuel Ratio')
xlabel('Time (s)')
ylabel('')
legend('A/F Results','A/F Target')
grid on

figure
plot(BurnTime,TotallImp)
title('Total Impulse')
xlabel('Time (s)')
ylabel('N*s')
grid on


figure
plot(BurnTime,Thrustdlvd)
title('Thrust Curve')
xlabel('Time (s)')
ylabel('Thrust (N)')
grid on
% ylim([0 4000])

figure
plot(BurnTime,MassGen)
title('Total Mass Generated')
xlabel('Time (s)')
ylabel('Mass Flow (kg)')
grid on
% ylim([0.01 0.03])

figure
plot(BurnTime,FuelMass)
title('Fuel Grain Mass Left')
xlabel('Time (s)')
ylabel('Mass (kg)')
grid on
% ylim([0 1.4])

figure
plot(BurnTime,MassFlow)
title('Total Mass Flow')
xlabel('Time (s)')
ylabel('Mass (kg)')
grid on

figure
plot(BurnTime,StepHeight)
title('Step Height')
xlabel('Time (s)')
ylabel('Height (m)')
grid on
ylim([0 .05])

fprintf('\n------------ Simulation Results ------------\n')
fprintf('Burn Time:             %.2f   (s)\n', BurnTime(n-1))
fprintf('Average Thrust:        %.2f  (N)\n', mean(Thrustdlvd))
fprintf('Total Impulse:         %.2f  (Ns)\n', TotallImp(n-1))
fprintf('Inlet Velocity:        %.2f   (m/s)\n', InltVel(n-1))
fprintf('Inlet Mass Flow Rate:  %.4f  (kg/s)\n', InltMassFlw)

