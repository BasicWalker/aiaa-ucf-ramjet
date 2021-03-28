%% Adiabatic flame temperature calculator %%
%% Using hydrocarbon fuels (C, H, N O) %%
%---------Notes---------%
%TO AUTOMATE:
%Make an array containing products and an array containing reactants
%Inputs required will include mole numbers for each molecule (stoichiometric combustion), molar mass of molecules, heats of formation of molecules, and molar Cp of products

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

    
    for(reacttemp = 400:100:1000)
        
        fprintf("------------REACTANT TEMPERATURE = %d K---------------", reacttemp);
        %---HEAT OF COMBUSTION CALCULATION---%

        %NEEDED:
        
        %Average molar Cp of products and reactants
        %Number of moles (total, reactants, products)
        %Temperature of reactants
        %Heats of formation of all molecules
        %Total moles of products
        %Mole fractions for each molecule
        
        cpair = 29.67 
        
        %HEAT OF FORMATION FOR EACH MOLECULE%

        %hf(molecule) = heat of formation for each molecule
        %Values found online

        hfn2 = 0;
        hfo2 = 0;
        hfco2 = -393520; %kJ/kmol
        hfh2o = -241820; %kJ/kmol
        hfabs = 62630; %kJ/kmol
        
        
        tstp = 273.15; %K
        
        %Molar Cp of molecules needed (CP VALUES CHANGE WITH TEMPERATURE,WILL BE INCORPORATED AT A LATER TIME)
        
        %Mass Cp values (kJ/kg*K)
        %Other values taken at 300K
        cpabs = 1.423; %Average value given, won't change with temperature
        
        %Molar Cp values
        mcpabs = cpabs*mwabs
        mcpair = cpair
        
        %Molar CP must be calculated for high temperatures
        x = reacttemp/100;
        mcph2o = 143.05 - (183.54 * (x^0.25)) + (82.751 * (x^ 0.5)) - (3.6989 * x)
        mcpco2 = -3.7357 + (30.529 * (x^0.5)) - (4.1034 * x) + (.024198 * (x^2))
        mcpn2 = 39.06 - (512.79 * (x^(-1.5))) + (1072.7 * (x^(-2))) - (820.4 * (x^(-3)))
        
        %Mole Fractions
        mfabs = molf/molr
        mfair = molair/molr
        mfh2o = molh2o/molp
        mfco2 = molco2/molp
        mfn2 = moln2/molp
        
        %Average Molar Cp values
        mcpr = (mfabs * mcpabs)  + (mfair * mcpair); %Reactants
        mcpp = (mfh2o * mcph2o) + (mfco2 * mcpco2) + (mfn2 * mcpn2)
        
        %HEAT OF COMBUSTION*
        
        %Long equation will be simplified into terms to eliminate chances of parentheses errors.
        
        %Followiing terms are number of moles multiplied by heat of formation
        %HEATS OF FORMATION FOR N2 AND O2 ARE ZERO
        
        %Products
        h2oterm = molh2o * hfh2o;
        co2term = molco2 * hfco2;
        
        
        %Reactants
        absterm = molf * hfabs;
        
        %Heat of reaction determined, then multiplied by -1 to get heat of combustion.
        
        molhrxn = h2oterm + co2term - absterm + (molp * mcpp * (reacttemp - tstp)) - (molr * mcpr * (reacttemp - tstp))
        
        hrxn = molhrxn/mwabs;
        
        hoc = hrxn * -1;
        
        fprintf("Heat of Combustion: %f kJ/kg\n", hoc)
        
       
        %---ADIABATIC FLAME TEMPERATURE---%
       
        
        %AFT Calculation
        
        aft = reacttemp + hoc/((molh2o * mcph2o) + (molco2 * mcpco2) + (moln2 * mcpn2));
        
        fprintf("ADIABATIC FLAME TEMP = %f K\n\n\n\n", aft);
       
        
        %---PLOTS---%
        
        figure(1)
        plot(reacttemp, aft, '*'); hold on
        grid on;
        xlabel('Reactant Temperature (K)');
        ylabel('Adiabatic Flame Temperature(K)');
        
%         figure(2)
%         plot(reacttemp, mcph2o, '*'); hold on
%         grid on;
%         xlabel('Temperature (K)');
%         ylabel('Molar Cp of H2O (kJ/kmol * K)');
%         
%         figure(3)
%         plot(reacttemp, mcpco2, '*'); hold on
%         grid on;
%         xlabel('Temperature (K)');
%         ylabel('Molar Cp of CO2 (kJ/kmol * K)');
%         
%         figure(4)
%         plot(reacttemp, mcpn2, '*'); hold on
%         grid on;
%         xlabel('Temperature (K)');
%         ylabel('Molar Cp of N2 (kJ/kmol * K)');
               
    end
    
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
 