%Ideal gas variable calculator

%----------ABS------------%

    molf = 1; %Moles of fuel (stoichiometric)
    molair = 24.107;
    molh2o = 2.425;
    molco2 = 3.85;
    moln2 = 19.26;
    
    molr = molf + molair;
    molp = molh2o + molco2 + moln2;
    
    %Molecular weights needed (kg/kmol)
    %Reactants
    mwabs = 57.07;
    mwair = 28.96;
    
    %Products
    mwh2o = 18.02;
    mwco2 = 44.01;
    mwn2 = 28.013;

    cpair = 29.67; 
    
for (temp = 2500:100:3500)
    
    fprintf("Temp = %d", temp);
     cpabs = 1.423; %Average value given, won't change with temperature

    %Molar Cp values
    mcpabs = cpabs*mwabs
    mcpair = cpair

    %Molar CP must be calculated for high temperatures
    x = temp/100;
    mcph2o = 143.05 - (183.54 * (x^0.25)) + (82.751 * (x^ 0.5)) - (3.6989 * x)
    mcpco2 = -3.7357 + (30.529 * (x^0.5)) - (4.1034 * x) + (.024198 * (x^2))
    mcpn2 = 39.06 - (512.79 * (x^(-1.5))) + (1072.7 * (x^(-2))) - (820.4 * (x^(-3)))
    
    %Molar basis to mass basis
    cph2o = mcph2o/mwh2o
    cpco2 = mcpco2/mwco2
    cpn2 = mcpn2/mwn2
    
    %Mass Fractions
    massh2o = molh2o * mwh2o;
    massco2 = molco2 * mwco2;
    massn2 = moln2 * mwn2;
    
    massp = massh2o + massco2 + massn2;
    
    mafh2o = massh2o/massp;
    mafco2 = massco2/massp;
    mafn2 = massn2/massp;
    
    %Average Molar Cp values
    mcpr = (mfabs * mcpabs)  + (mfair * mcpair); %Reactants
    mcpp = (mfh2o * mcph2o) + (mfco2 * mcpco2) + (mfn2 * mcpn2);
    
    
    %Average Mass Cp value
    macpp = (mafh2o * cph2o) + (mafco2 * cpco2) + (mafn2 * cpn2)

    %---IDEAL GAS VALUES---%

    Ru = 8.314; %kJ/kmol * K
    R_abs = Ru/mwabs;
    R_air = Ru/mwair;
    R_h2o = Ru/mwh2o
    R_co2 = Ru/mwco2
    R_n2 = Ru/mwn2

    Cv_abs = mcpabs - R_abs;
    Cv_air = mcpair - R_air;
    Cv_h2o = cph2o - R_h2o;
    Cv_co2 = cpco2 - R_co2;
    Cv_n2 = cpn2 - R_n2;

    macvp = (mafh2o * Cv_h2o) + (mafco2 * Cv_co2) + (mafn2 * Cv_n2);

    k = macpp/macvp;

    figure(5)
    plot(temp, macpp, '*'); hold on
    grid on;
    xlabel('Product Temperature (K)');
    ylabel('Average Cp of Products (kJ/kg * K)');

    figure(6)
    plot(temp, macvp, '*'); hold on
    grid on;
    xlabel('Product Temperature (K)');
    ylabel('Average Cv of Products (kJ/kg * K)');

    figure(7)
    plot(temp, k, '*'); hold on
    grid on;
    xlabel('Temperature (K)');
    ylabel('Specific Heat Ratio of Products');
    
end
 