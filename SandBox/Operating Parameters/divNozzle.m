function [exitMach, exitPres, exitVel, exitTemp] = divNozzle(expansionRatio,throatPres, throatTemp)
%Divergent Nozzle Calculations 
%   Detailed explanation goes here

gamma = 1.4;
R = 287.1;


[m_e, t_e, p_e, d_e, a_e] = flowisentropic(gamma, expansionRatio, 'sup');

exitTemp = throatTemp * t_e;

exitA = sqrt(gamma * R * exitTemp);

exitMach = m_e;

exitPres = p_e * throatPres;

exitVel = exitMach * exitA;

exitTemp = throatTemp * t_e;

end

