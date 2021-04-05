function [throatArea,throatPres,throatTemp,Vel] = nozzleThroatArea(chamberPres,chamberVel,chamberTemp,chamberArea)
%nozzleThroatArea Minimum Nozzle Throat Area Finder Function
%   Detailed explanation goes here

gamma = 1.4;
R = 287.1;

a = sqrt(gamma*R*chamberTemp);
mach = chamberVel/a;

[m_i, t_i, p_i, d_i, a_i] = flowisentropic(gamma, mach, 'mach');

throatArea = chamberArea / a_i;
throatPres = chamberPres / p_i;
throatTemp = chamberTemp / t_i;
Vel = chamberVel;



end