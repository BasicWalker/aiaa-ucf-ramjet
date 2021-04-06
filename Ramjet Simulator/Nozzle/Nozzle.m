% ------ SFRJ Internal Ballistic Simulator / UCF CAPSTONE PROJECT ------ %
% File Name: Nozzle.m 
% 
% File Description: 
% Numerical solver to design the optimal nozzle geometry
% 
% Name            Date      SCR  Description
% --------------  --------  ---  ------------------------------
% Ethan Sherlock  01/22/21  000  Initial Creation 
% ---------------------------------------------------------------------- %


nozzle.Area_ratio(2) = nozzle.Area_exit/nozzle.Area_throat;

% nozzle.mach(1,n) = combustion.mach(2,n);
nozzle.stagTemp(1,n) = combustion.stagTemp(2,n);
nozzle.stagDens(1,n) = combustion.stagDens(2,n);
nozzle.stagPres(1,n) = combustion.stagPres(2,n);

[nozzle.mach(1,n), nozzle.tempRatio(1,n), nozzle.presRatio(1,n), nozzle.densRatio(1,n), nozzle.areaRatio(1)]...
    = flowisentropic(gamma, 1, 'mach');  % ratios are static over stagnation

nozzle.staticTemp(1,n) = nozzle.tempRatio(1,n)*nozzle.stagTemp(1,n);
nozzle.staticDens(1,n) = nozzle.densRatio(1,n)*nozzle.stagDens(1,n);
nozzle.staticPres(1,n) = nozzle.presRatio(1,n)*nozzle.stagPres(1,n);
nozzle.velocity(1,n) = sqrt(gamma*R*nozzle.staticTemp(1,n));
nozzle.massFlow(1,n) = nozzle.staticDens(1,n)*nozzle.Area_throat*nozzle.velocity(1,n);

[nozzle.mach(2,n), nozzle.tempRatio(2,n), nozzle.presRatio(2,n), nozzle.densRatio(2,n), ~]...
    = flowisentropic(gamma, nozzle.Area_ratio(2), 'sup');  % ratios are static over stagnation

nozzle.stagTemp(2,n) = nozzle.stagTemp(1,n);
nozzle.stagDens(2,n) = nozzle.stagDens(1,n);
nozzle.stagPres(2,n) = nozzle.stagPres(1,n);
nozzle.staticTemp(2,n) = nozzle.tempRatio(2,n)*nozzle.stagTemp(2,n);
nozzle.staticDens(2,n) = nozzle.densRatio(2,n)*nozzle.stagDens(2,n);
nozzle.staticPres(2,n) = nozzle.presRatio(2,n)*nozzle.stagPres(2,n);
nozzle.velocity(2,n) = nozzle.mach(2,n)*sqrt(gamma*R*nozzle.staticTemp(2,n));
nozzle.massFlow(2,n) = nozzle.staticDens(2,n)*nozzle.Area_exit*nozzle.velocity(2,n);


