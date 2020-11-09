%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            AIAA UCF Ramjet Intake Geometry Design Script                %
%                                                                        %
%                              Authored by                                %
%        Samer Armaly, Karam Paul, Jared Durlak, Matthew Aubertin         %
%                           October 28, 2020                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

%%
defl = [1 1 1];  % deflection angle array
n = size(defl,2);  % number of oblique shocks


% pre-allocate arrays
mach = zeros(n+2,1);
shockAngle = zeros(n,1);

%input known terms 
mach(1) = 1.5;  % freestream mach
%ctrbdy_dia = 1;  % max diameter for centerbody (inches)
cowl_inner_dia = 2.5;  % inner diameter for cowl (inches)






%% oblique shocks
%solving for oblique shocks
for i = 1:n
[mach(i+1), shockAngle(i)] = obliqueShock(mach(i), defl(i));
end

%% normal shock
%solve for normal shock





