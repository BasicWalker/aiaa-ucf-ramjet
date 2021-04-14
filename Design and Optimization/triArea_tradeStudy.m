clc;clear;
if exist('GRAM','var')==0                      % GRAM atmospheric model
    load GRAM_Model.mat                     
end
def = 5:1:10;  % Deflection angle (deg)
intake.Area_enter = 0.0005:0.0001:0.0009;  % Area of throat (m^2) - Drives mass flow rate through intake
nozzle.Area_throat = 0.001:0.001:0.0508;  % Throat area (m^2)
nozzle.Area_exit =   0.0508:0.001:0.0698;  % nozzle exit area (m^2)

DesignThrust = zeros(numel(def),numel(intake.Area_enter),numel(nozzle.Area_throat),numel(nozzle.Area_exit));
DesignintakemassFlow = zeros(numel(def),numel(intake.Area_enter),numel(nozzle.Area_throat),numel(nozzle.Area_exit));
DesigncombustionStaticPres = zeros(numel(def),numel(intake.Area_enter),numel(nozzle.Area_throat),numel(nozzle.Area_exit));
DesigncombustionStagPres = zeros(numel(def),numel(intake.Area_enter),numel(nozzle.Area_throat),numel(nozzle.Area_exit));


% --------------------- Constants and Conversion Variables ---------------------- %
gravity         = 9.81;                     % gravitation acceleration constant (m/s^2)
In2Mtr          = 39.3701;                  % Inch to meter converstion 
Bar2kPa         = 100.0;                    % Bar to kPa conversion
Pa2kPa          = 1000.0;                   % Pa to Kpa
C2K             = 273.15;                   % Celcius to Kelvin conversion
R               = 287.05;                   % Universal Gas Constant for air
% --------------- Environmental User Defined Parameters --------------- %
flight_mach(1)  = 2;                      % Booster max mach
altitude(1)     = 5000;                     % Initial altitude for ramjet start (m)
c_d             = 0.23;                     % Drag coefficient (0.35)
S               = 0.008119;                 % Frontal surface area (m^2)
gamma_atm       = 1.4;                      % Specific heat ratio
dry_mass        = 6.80389;                  % Mass of ramjet without fuelgrain (kg)
OxPercent       = 0.2314;                   % Density percentage of oxygen in air by mass
% --------------- Fuel Grain User Defined Parameters --------------- %
GrainOD         =  2.75 /In2Mtr;                        % Grain OD (m)
GrainID(1)      =  2.00 /In2Mtr;                        % Grain ID (m)
GrainL          = 15.00 /In2Mtr;                        % Grain Length (m)
FuelRho         = 1020;                                 % Grain Density (kg/m^3)
PortArea(1)     = pi*(GrainID(1)^2)*(1/4);              % Fuel Port Area (m^2)
FuelCS(1)       = pi*(GrainOD^2)*(1/4) - PortArea(1);   % Fuel Grain Crossectional Area (m^2)
FuelVol(1)      = FuelCS(1) * GrainL;                   % Fuel Grain Volume (m^3)
FuelSA(1)       = GrainID(1)* pi * GrainL;              % Fuel Grain Surface Area (m^2)
FuelMass(1)     = FuelRho*FuelVol(1);                   % Grain fuel mass, instantaneous (kg)
MdotFuel = 0.09;
% ----------------- Combustor User Defined Parameters ----------------- %
InltD           = 1.4 / In2Mtr;            % Diameter of Combustor inlet (m)
gamma_Inlt      = 1.3845;                   % Specific heat ratio of air 
radius_combustor= InltD/2;                  % Radius of the combustor inlet (m)
Area_combustor  = pi*radius_combustor^2;    % Area of the combustor inlet (m^2)
gamma           = 1.4;                      % Specific heat ratio (atm)
% --------------- Chemistry User Defined Parameters ---------------- %
chem = Chemistry();
% ----------------- Trajectory Initial Conditions ------------------ %                  
Rho_atm(1)      = interp1(GRAM.Hgtkm, GRAM.DensMean, (altitude(1))/1e3);    % Atmospheric Density (kg/m^3)
pressure_atm(1) = interp1(GRAM.Hgtkm, GRAM.PresMean, (altitude(1))/1e3);    % Atmospheric Pressure (Pa)
pressure_atm(1) = pressure_atm(1)*(1/Pa2kPa);                               % Atmospheric Pressure (kPa)
Temp_atm(1)     = interp1(GRAM.Hgtkm, GRAM.Tmean, (altitude(1))/1e3);       % Atmospheric Temperature (K)
velocity(1)     = flight_mach(1)*sqrt(gamma_atm*R*Temp_atm(1));             % Atmospheric Velocity (m/s)
drag(1)         = c_d*0.5*Rho_atm(1)*velocity(1)^2*S;                       % Induced Drag (N)
Thrustdlvd(1)   = 0.0;                                                      % Initialize First Thrust Value 
mass(1)         = dry_mass + FuelMass(1);                                   % Mass of Vehicle
weight(1)       = gravity*mass(1);                                          % Weight of Vehicle
acceleration(1) = (Thrustdlvd(1) - drag(1) - weight(1))/ mass(1);           % Initial Acceleration
SFRJDt          = 1/100;                    % Simulation rate (Hz) 
n = 2;
for i=1:numel(def)
    for j=1:numel(intake.Area_enter)
        for k=1:numel(nozzle.Area_throat)
            for q=1:numel(nozzle.Area_exit)
%                 try
                    
                    
                    % *****intake*****
                    
                    
                    % Station 1 properties (free-stream)-----
                    intake.mach(1) = flight_mach;
                    intake.staticPres(1) = pressure_atm*1e3;  % <Pa>
                    intake.staticDens(1) = Rho_atm;  % <kg/m^3>
                    intake.staticTemp(1) = Temp_atm;  % <K>
                    [~, intake.tempRatio(1), intake.presRatio(1), intake.densRatio(1), ~] = flowisentropic(gamma, intake.mach(1), 'mach');  % ratios are static over stagnation
                    intake.stagTemp(1) = intake.staticTemp(1)/intake.tempRatio(1);  % <K>
                    intake.stagPres(1) = intake.staticPres(1)/intake.presRatio(1);  % <Pa>
                    intake.stagDens(1) = intake.staticDens(1)/intake.densRatio(1);  % <kg/m^3>
                    intake.velocity(1) = intake.mach(1)*sqrt(gamma*R*intake.staticTemp(1));  % <m/s>
                    
                    % Station 2 properties (oblique shock)-----
                    [intake.mach(2), intake.shockAngle] = obliqueShock(intake.mach(1), def(i), gamma);  % uniform mach number of external ramp region
                    intake.mach1Norm = intake.mach(1)*sind(intake.shockAngle);  % flow component crossing normal to oblique shock
                    [~, intake.tempRatio(2), intake.presRatio(2), intake.densRatio(2), ~, intake.stagPresRatio(2)] = ...
                        flownormalshock(gamma, intake.mach1Norm, 'mach');  % ratios are downstream over upstream
                    intake.stagPres(2) = intake.stagPres(1) * intake.stagPresRatio(2);  % <Pa>
                    intake.stagDens(2) = intake.stagDens(1) * intake.stagPresRatio(2);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
                    intake.stagTemp(2) = intake.stagTemp(1);  % does not change over normal shock <K>
                    intake.staticPres(2) = intake.staticPres(1) * intake.presRatio(2);  % <Pa>
                    intake.staticDens(2) = intake.staticDens(1) * intake.densRatio(2);  % <kg/m3>
                    intake.staticTemp(2) = intake.staticTemp(1) * intake.tempRatio(2);  % <K>
                    intake.velocity(2) = intake.mach(2)*sqrt(gamma*R*intake.staticTemp(2));  % <m/s>
                    [~, ~, ~, ~, intake.A_Astar(2)] = flowisentropic(gamma, intake.mach(2), 'mach');
                    intake.Astar(2) = intake.Area_enter(j)/intake.A_Astar(2);  % reference area before normal shock

                    DesignintakemassFlow(i,j,k,q) = intake.staticDens(2)*intake.Area_enter(j)*intake.velocity(2);  % air mass flow at intake opening <kg/s>

                    intake.chokeStagPres = pressureToChoke(DesignintakemassFlow(i,j,k,q), nozzle.Area_throat(k), intake.stagTemp(2));  % no stag temp change from combustion

                    % check for choked conditions
                    if intake.stagPres(2) < intake.chokeStagPres
                        error('intake:obliqueStagLoss',...
                            'Error. \nOblique shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry decreasing altitude, deflection angle or increasing speed.');
                    end
                    [~, ~, ~, ~, ~, intake.r_d_chokeLimit] = flownormalshock(gamma, intake.mach(2), 'mach');  % ratios are downstream over upstream
                    if (intake.stagPres(2)*intake.r_d_chokeLimit) < intake.chokeStagPres
                        error('intake:normalStagLoss',...
                            'Error. \nNormal shock stagnation pressure loss too great\ncannot provide enough stagnation pressure; flow is no longer choked. \nTry reducing the strength of the normal shock.');
                    end
                    % Station 3 properties (normal shock)-----
                    % solve normal shock placement
                    intake.r_d = intake.chokeStagPres/intake.stagPres(2);
                    [intake.machShock, intake.tempRatio(3), intake.presRatio(3), intake.densRatio(3), intake.mach(3), intake.stagPresRatio(3), ~]...
                        = flownormalshock(gamma, intake.r_d, 'totalp');
                    if intake.machShock < intake.mach(2)
                        error('intake:shockEjected',...
                            'Error. \nNormal shock ejected; losing mass flow. \ndiffuser ratio (r_d) too high.');
                    end
                    [~, ~, ~, ~, intake.Ashock_A2star] = flowisentropic(gamma, intake.machShock, 'mach');
                    intake.AShock = intake.Ashock_A2star*intake.Astar(2);   % area of where the normal shock occurs
                    
                    intake.stagPres(3) = intake.stagPres(2) * intake.stagPresRatio(3);  % <Pa>
                    intake.stagDens(3) = intake.stagDens(2) * intake.stagPresRatio(3);  % changes with the same rate as pressure (ideal gas law) <kg/m3>
                    intake.stagTemp(3) = intake.stagTemp(2);  % does not change over normal shock <K>
                    intake.staticPres(3) = intake.staticPres(2) * intake.presRatio(3);  % <Pa>
                    intake.staticDens(3) = intake.staticDens(2) * intake.densRatio(3);  % <kg/m3>
                    intake.staticTemp(3) = intake.staticTemp(2) * intake.tempRatio(3);  % <K>
                    intake.velocity(3) = intake.mach(3)*sqrt(gamma*R*intake.staticTemp(3));  % <m/s>
                    [~, ~, ~, ~, intake.A_Astar(3)] = flowisentropic(gamma, intake.mach(3), 'mach');
                    intake.Astar(3) = intake.Astar(2)/intake.stagPresRatio(3);  % reference area after normal shock
                    % intake.AShock = intake.A_Astar(3)*intake.Astar(3);
                    intake.massFlow(3) = intake.staticDens(3)*intake.AShock*intake.velocity(3);
                    
                    % Station 4 properties (subsonic iscentropic expansion)-----
                    intake.A_Astar(4) = Area_combustor/intake.Astar(3);  % Area ratio of the combustor inlet
                    [intake.mach(4), intake.tempRatio(4), intake.presRatio(4), intake.densRatio(4), ~] = ...
                        flowisentropic(gamma, intake.A_Astar(4), 'sub');  % ratios are static over stagnation
                    intake.stagPres(4) = intake.stagPres(3);  % does not change; iscentropic <Pa>
                    intake.stagDens(4) = intake.stagDens(3);  % does not change; iscentropic <kg/m3>
                    intake.stagTemp(4) = intake.stagTemp(3);  % does not change; iscentropi <K>
                    intake.staticPres(4) = intake.stagPres(4) * intake.presRatio(4);  % <Pa>
                    intake.staticDens(4) = intake.stagDens(4) * intake.densRatio(4);  % <kg/m3>
                    intake.staticTemp(4) = intake.stagTemp(4) * intake.tempRatio(4);  % <K>
                    intake.velocity(4) = intake.mach(4)*sqrt(gamma*R*intake.staticTemp(4));  % <m/s>
                    intake.massFlow(4) = intake.staticDens(4)*Area_combustor*intake.velocity(4);
                    
                    % Pass variables to ballistic simulator
                    InltPres_stag = intake.stagPres(4);  % Inlet stagnation pressure
                    InltTemp_stag = intake.stagTemp(4);  % Inlet stagnation temperature
                    InltTemp_dens = intake.stagDens(4);  % Inlet stagnation density
                    InltPres = intake.staticPres(4);  % Inlet static pressure
                    InltRho = intake.staticDens(4);  % Inlet static density
                    InltTemp = intake.staticTemp(4);  % Inlet static temp
                    InltVel = intake.velocity(4);  % Inlet velocity
                    Inltmach = intake.mach(4);  % Inlet mach number
                    
                    
                    % *****RegressionRate*****
                    
                    
                    BurnRt     = 0.001;         % Burn Rate (m/s)
                    RgrsPerStp = BurnRt*SFRJDt; % Regression Per Step (m)
                    
                    
                    % *****GrainGeometry*****
                    
                    

                    GrainID = GrainID(n-1) + 2 * RgrsPerStp;     % Fuel Grain Inner Diameter (m)
                    PortArea = pi*(GrainID^2)*(1/4);          % Fuel Port Area (m^2)
                    FuelCS = pi*(GrainOD^2)*(1/4) - PortArea; % Fuel Grain Crossectional Area (m^2)
                    FuelVol = FuelCS * GrainL;                % Fuel Grain Volume (m^3)
                    FuelSA = GrainID* pi * GrainL;            % Fuel Grain Surface Area (m^2)

                    
                    StepHeight = (GrainID - InltD)/2;             % Rearward Step Height
                    
                    % Fuel Mass Properties
                    MFuelGen = RgrsPerStp*FuelRho*FuelSA;         % Fuel mass generated every time step (kg)
                    MdotFuel = MFuelGen/SFRJDt;                   % Fuel mass flow rate (kg/s)
                    FuelMass = FuelRho*FuelVol;                   % Grain fuel mass, instantaneous (kg)
                    
                    
                    % *****Gas*****
                    
                    
                    % Calculate air mass properties
                    MairGen = DesignintakemassFlow(i,j,k,q)*SFRJDt;                                     % Mass of air generated (kg)
                    MOxdzrGen = MairGen * OxPercent;                              % Mass of oxygen generated by weight % (kg)
                    MdotTotal = DesignintakemassFlow(i,j,k,q) + MdotFuel;                            % Total mass flow (kg/s)
                    f_yield = (MFuelGen/MairGen);                              % fuel to air ratio
                    [phi, T_stag, gamma_T, R_t] = chem.phiSolver(f_yield,intake.staticTemp(4));     % Call Chemistry Model, need to add T_air before combustion chamber 475
                    gamma_t = gamma_T;                                               % Gamma at the nozzle throat (old:gamma_nzl = 1.3845)
                    phi_eqv = phi;                                                   % Grab phi value
                    % Fuel Mass Prop & O/F Ratio calculation
                    OFRatio = MOxdzrGen/MFuelGen;                              % O/F Ratio
                    MassGen = MairGen + MFuelGen;                              % Total mass generated (kg)
                    MassFlow = MassGen/SFRJDt;                                    % Total mass flow (kg/s)
                    
                    
                    % *****combustor*****
                    
                    
                    [~, ~, ~, ~, combustion.A_Astar(1)] = flowisentropic(gamma, Inltmach, 'mach');
                    combustion.Astar(1) = Area_combustor / combustion.A_Astar(1);
                    combustion.A_Astar(2) = PortArea / combustion.Astar(1);
                    %
                    %
                    [combustion.mach(1), combustion.TempRatio(1), combustion.PresRatio(1), combustion.DensRatio(1), ~] = ...
                        flowisentropic(gamma, combustion.A_Astar(2), 'sub');  % ratios are static over stagnation
                    combustion.stagPres(1) = InltPres_stag;  % does not change; iscentropic <Pa>
                    combustion.stagDens(1) = InltTemp_dens;  % does not change; iscentropic <kg/m3>
                    combustion.stagTemp(1) = InltTemp_stag;  % does not change; iscentropi <K>
                    combustion.staticPres(1) = combustion.stagPres(1) * combustion.PresRatio(1);  % <Pa>
                    combustion.staticDens(1) = combustion.stagDens(1) * combustion.DensRatio(1);  % <kg/m3>
                    combustion.staticTemp(1) = combustion.stagTemp(1) * combustion.TempRatio(1);  % <K>
                    combustion.velocity(1) = combustion.mach(1)*sqrt(gamma*R*combustion.staticTemp(1));  % <m/s>
                    combustion.massFlow(1) = combustion.staticDens(1)*combustion.velocity(1)*PortArea;
                    
                    % use rayleigh flow to solve for properties at combustion chamber exit
                    [combustion.mach(2), combustion.staticTemp(2), DesigncombustionStaticPres(i,j,k,q), combustion.stagTemp(2), DesigncombustionStagPres(i,j,k,q), combustion.stagPresLoss] = ...
                        RayleighFlow(gamma, combustion.mach(1), InltTemp_stag, T_stag, InltPres_stag);
                    
                    combustion.staticDens(2) = DesigncombustionStaticPres(i,j,k,q)/(combustion.staticTemp(2)*R);
                    combustion.stagDens(2) = DesigncombustionStagPres(i,j,k,q)/(combustion.stagTemp(2)*R);
                    combustion.velocity(2) = combustion.mach(2)*sqrt(1.2*600*combustion.staticTemp(2));
                    combustion.massFlow(2) = combustion.staticDens(2)*combustion.velocity(2)*PortArea + MdotFuel;  % add in the fuel mass flow to rayleigh heated air mass flow
                    
                    
                    % *****Nozzle*****
                    
                    
                    nozzle.Area_ratio(2) = nozzle.Area_exit(q)/nozzle.Area_throat(k);
                    
                    % nozzle.mach(1) = combustion.mach(2);
                    nozzle.stagTemp(1) = combustion.stagTemp(2);
                    nozzle.stagDens(1) = combustion.stagDens(2);
                    nozzle.stagPres(1) = DesigncombustionStagPres(i,j,k,q);
                    
                    [nozzle.mach(1), nozzle.tempRatio(1), nozzle.presRatio(1), nozzle.densRatio(1), nozzle.areaRatio(1)]...
                        = flowisentropic(gamma, 1, 'mach');  % ratios are static over stagnation
                    
                    nozzle.staticTemp(1) = nozzle.tempRatio(1)*nozzle.stagTemp(1);
                    nozzle.staticDens(1) = nozzle.densRatio(1)*nozzle.stagDens(1);
                    nozzle.staticPres(1) = nozzle.presRatio(1)*nozzle.stagPres(1);
                    nozzle.velocity(1) = sqrt(gamma*R*nozzle.staticTemp(1));
                    nozzle.massFlow(1) = nozzle.staticDens(1)*nozzle.Area_throat(k)*nozzle.velocity(1);
                    
                    [nozzle.mach(2), nozzle.tempRatio(2), nozzle.presRatio(2), nozzle.densRatio(2), ~]...
                        = flowisentropic(gamma, nozzle.Area_ratio(2), 'sup');  % ratios are static over stagnation
                    
                    nozzle.stagTemp(2) = nozzle.stagTemp(1);
                    nozzle.stagDens(2) = nozzle.stagDens(1);
                    nozzle.stagPres(2) = nozzle.stagPres(1);
                    nozzle.staticTemp(2) = nozzle.tempRatio(2)*nozzle.stagTemp(2);
                    nozzle.staticDens(2) = nozzle.densRatio(2)*nozzle.stagDens(2);
                    nozzle.staticPres(2) = nozzle.presRatio(2)*nozzle.stagPres(2);
                    nozzle.velocity(2) = nozzle.mach(2)*sqrt(gamma*R*nozzle.staticTemp(2));
                    nozzle.massFlow(2) = nozzle.staticDens(2)*nozzle.Area_exit(q)*nozzle.velocity(2);
                    
                    
                    % *****Thrust*****
                    DesignThrust(i,j,k,q) = (DesignintakemassFlow(i,j,k,q)*intake.velocity(2)) - nozzle.massFlow(2)*nozzle.velocity(2);


%                 catch ME
%                     fprintf('something went wrong:( %s\n', ME.message);
%                     DesignThrust(i,j,k,q) = NaN;
%                     DesignintakemassFlow(i,j,k,q) = NaN;
%                     DesigncombustionStaticPres(i,j,k,q) = NaN;
%                     DesigncombustionStagPres(i,j,k,q) = NaN;
%                 end
            end
        end
    end
end





