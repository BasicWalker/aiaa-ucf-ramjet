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

fprintf('Simulation Complete \n')

figure('Name','O/F Ratio')
plot(BurnTime(1:index-1),combustion.OFRatio(1:index-1))
title('O/F Ratio')
xlabel('Time (s)')
ylabel('')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Force Profile')
plot(BurnTime,vehicle.thrust,BurnTime,trajectory.F_d, BurnTime, vehicle.TotalMass*constants.gravity, BurnTime, trajectory.F_net)
title('Thrust & Drag vs Time')
xlabel('Time (s)')
ylabel('Force (N)')
legend('Thrust Curve','Drag','Weight','Net Force')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Mass Flow Rate')
plot(BurnTime(1:index-1),vehicle.TotalMassFlow(1:index-1), BurnTime(1:index-1), intake.massFlow(end,1:index-1), BurnTime(1:index-1), fuel.MassFlow(1:index-1)) 
title('Mass Flow Rate')
xlabel('Time (s)')
ylabel('Mass (kg/s)')
legend('Total Mdot','Mdot Air','Mdot Fuel')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Pressure Plots')
plot(BurnTime(1:index-1), intake.stagPres(1,1:index-1)/1e3, BurnTime(1:index-1), intake.chokeStagPres(1:index-1)/1e3,'*', BurnTime(1:index-1), combustion.stagPres(1,1:index-1)/1e3, BurnTime(1:index-1), trajectory.pressure_a(1:index-1)/1e3,BurnTime(1:index-1), intake.stagPres(2,1:index-1)/1e3,BurnTime(1:index-1), intake.stagPres(4,1:index-1)/1e3,BurnTime(1:index-1),intake.MaxAvailableStag(1:index-1)/1e3)%,BurnTime(1:index-1),intake.chokeStagPresTEST(1:index-1)/1e3)
title('Stagnation Pressure Plots')
xlabel('Time (s)')
ylabel('Stagnation Pressure (kPa)')
legend('Freestream Pressure','Choked Flow Pressure Required', 'Chamber Pressure (before stag loss)','Back Pressure','after oblique','after normal','max stag av')%,'test choke')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Flight Mach')
plot(BurnTime,vehicle.Mach)
title('Vehicle Mach vs Time')
xlabel('Time (s)')
ylabel('Mach No.')
grid on
set(gcf,'position',[550,200,800,700])

figure('Name','Equivalence Ratio')
plot(BurnTime(1:index-1),combustion.Phi(1:index-1))
title('Equivalence Ratio vs Time')
xlabel('Time (s)')
ylabel('Equivalence Ratio (phi)')
grid on
% ylim([0 6])
set(gcf,'position',[550,200,800,700])

figure('Name','Adiabatic Flame Temperature')
yyaxis left
plot(BurnTime(1:index-1),vehicle.AFT(1:index-1))
title('Adiabatic Flame Temperature & Specific Heat Ratio')
xlabel('Time (s)')
ylabel('Temperature <k>')
yyaxis right
plot(BurnTime(1:index-1),nozzle.gamma(1:index-1))
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

%%    geometry plots
% solve for ramp length
spikeLength = intake.EnterDiaINCH/2/sind(intake.DeflAngle);
subDiffuserLength = spikeLength/2;
totalLength = spikeLength + subDiffuserLength;
subDiffuserSlope = -(intake.EnterDiaINCH/2)/subDiffuserLength;
topslopeCowlx = totalLength + ((intake.CowlDiaINCH/2)-(combustion.InletDiaINCH/2))/subDiffuserSlope;  % used to match slope of sub diffuser line to intersect with cowl
buffer = totalLength/10;
% solve for shocks
obliqueShockHgt = sind(intake.shockAngle(1))*spikeLength;

f1 = figure('Name','ramp shock structure');
hold on
% plot surfaces
plot([0,totalLength+buffer],[0,0],':k', 'LineWidth',2);  % centerline
plot([spikeLength,totalLength+buffer],[(intake.CowlDiaINCH/2),(intake.CowlDiaINCH/2)],'k', 'LineWidth',2);  % cowl line
plot([0,spikeLength],[0,(intake.EnterDiaINCH/2)], 'k', 'LineWidth',2);  % spike ramp line
plot([spikeLength,totalLength],[(intake.EnterDiaINCH/2),0], 'k', 'LineWidth',2);  % subsonic diffuser ramp line
plot([topslopeCowlx,totalLength],[(intake.CowlDiaINCH/2),(combustion.InletDiaINCH/2)], 'k', 'LineWidth',2);  % subsonic top slope ramp line
plot([totalLength,totalLength],[(combustion.InletDiaINCH/2),(intake.CowlDiaINCH/2)], 'k', 'LineWidth',2);  % combustor inlet

% plot shocks
plot([0,spikeLength],[0,obliqueShockHgt], '--c');  % oblique shock
plot([spikeLength,spikeLength],[(intake.EnterDiaINCH/2),(intake.CowlDiaINCH/2)], '--c');  % Normal Shcok
title('Ramjet Intake Geometry')
xlabel('<in>')
ylabel('<in>')
hold off



%%



if (mean(combustion.mach(1)) > 0.2)
    fprintf(2,'\nWARNING: Inlet Mach number is too high.\n')
end
if (mean(combustion.Phi) > 3.0)
    fprintf(2,'WARNING: Equivalence ratio exceeds flamability limit.\n')
end
if (vehicle.Mach(n-1) < 2)
    fprintf(2,'WARNING: Vehicle does not reach target Mach Number.\n')
end

fprintf('------------ Simulation Results ------------\n')
fprintf('Burn Time:                  %.2f    (s)\n', BurnTime(n-1))
fprintf('Average Thrust:             %.2f   (N)\n', mean(vehicle.thrust))
fprintf('Average Drag Force:         %.2f   (N)\n',mean(trajectory.F_d))
fprintf('Total Impulse:              %.2f  (Ns)\n', max(vehicle.Total_Impulse))
fprintf('Air Mass Flow Rate:         %.3f    (kg/s)\n', mean(intake.massFlow(end)))
fprintf('Chamber Pressure From TAFT: %.2f   (kPa)\n',mean(combustion.stagPres(end))/1e3)
fprintf('Initial Step Height:        %.2f     (in) \n', fuel.StepHeight(1)*constants.In2Mtr)
fprintf('Average Inlet Velocity:     %.2f   (m/s)\n', mean(intake.velocity(end)))
fprintf('Max Altitude(During Boost): %.2f (m)\n',max(trajectory.Z_pos))
fprintf('Average Inlet Mach:         %.2f  \n', mean(intake.mach(end)))
fprintf('Average O/F Ratio:          %.2f  \n',mean(combustion.OFRatio(2:end)))
fprintf('Average Equivalence Ratio:  %.2f  \n',mean(combustion.Phi))
fprintf('Average Gamma Ratio:        %.2f  \n',mean(nozzle.gamma(2:end)))
fprintf('--------------------------------------------\n')
toc