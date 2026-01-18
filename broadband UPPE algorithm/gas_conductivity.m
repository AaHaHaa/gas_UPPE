function conductivity = gas_conductivity(pressure_ratio,material,temperature)
%GAS_CONDUCTIVITY Summary of this function goes here
%   Detailed explanation goes here

num_gas = length(material);
conductivity = 0; % initialization
for gas_i = 1:num_gas
    switch material{gas_i}
        case 'H2'
            % Colozza and Jakupca, "Thermal System Sizing Comparison of a
            % PEM and Solid Oxide Fuel Cell Systems on Mars," (2019)
            conductivity_i = -5.2500e-3 + 7.8560e-4*(temperature-273.15) - 4.9000e-7*(temperature-273.15)^2; % W/m/K
        case 'D2'
            % Saxena and Saxena, "Thermal conductivity data for hydrogen
            % and deuterium in the range 100-1100 degrees C," J. Phys. A:
            % Gen. Phys. 3(3), 309-320 (1970)
            %
            % This has H2, D2, and N2.
            conductivity_i = (16.8 + 5.68e-2*temperature + 2.354e-6*temperature^2)/1e5; % cal/cm/s/K
            conductivity_i = conductivity_i*4.184/0.01; % W/m/K
        case 'N2'
            % Saxena and Chen, "Thermal conductivity of nitrogen in the
            % temperature range 350--2500 K," Mol. Phys. 29(5), 1507-1519
            % (1975)
            conductivity_i = (12.18 + 0.05224*temperature - 0.6482e-6*temperature^2 - 0.2765e-9*temperature^3)*1e-3; % W/m/K
        case 'O2'
            % Colozza and Jakupca, "Thermal System Sizing Comparison of a
            % PEM and Solid Oxide Fuel Cell Systems on Mars," (2019)
            conductivity_i = 4.8000e-4 + 9.2414e-5*(temperature-273.15) - 2.2143e-8*(temperature-273.15)^2; % W/m/K
        case 'air'
            % The Engineering ToolBox (2009). Air Properties - Thermal
            % Conductivity vs. Temperature and Pressure Charts and
            % Calculator.
            % [online] Available at:
            % https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
            T = [-190,-150,-100,-75,-50,-25,-15,-10,-5,0,5,10,15,20,25,30,40,50,60,80,100,125,150,175,200,225,300,412,500,600,700,800,900,1000,1100]' + 273.15; % K
            kappa = [7.82,11.69,16.20,18.34,20.41,22.41,23.20,23.59,23.97,24.36,24.74,25.12,25.50,25.87,26.24,26.62,27.35,28.08,28.80,30.23,31.62,33.33,35.00,36.64,38.25,39.83,44.41,50.92,55.79,61.14,66.32,71.35,76.26,81.08,85.83]'*1e-3; % W/m/K
            % Fit to a cubic polynomial
            p = polyfit(T,kappa,3);
            conductivity_i = polyval(p,temperature);
        case 'Xe'
            % Saxena and Saxena, "Thermal Conductivity of Krypton and Xenon
            % in the Temperature Range 350--1500°K," J. Chem. Phys. 51(8), 
            % 3361-3368 (1969)
            conductivity_i = (0.2772 + 0.4147e-2*temperature - 0.5748e-6*temperature^2)/1e5; % cal/cm/s/K
            conductivity_i = conductivity_i*4.184/0.01; % W/m/K
        case 'Ar'
            % Saxena and Chen, "Thermal conductivity of argon in the
            % temperature range 350--2500 K," Mol. Phys. 29(2), 455-466
            % (1975)
            conductivity_i = 5.465 + 0.04729*temperature - 0.1111e-4*temperature^2 + 0.1599e-8*temperature^3; % W/m/K
        case 'Ne'
            % Springer and Wingeier, "Thermal conductivity of neon, argon,
            % and xenon at high temperatures," J. Chem. Phys. 59(5),
            % 2747-2750 (1973)
            conductivity_i = 9.683e-4*temperature^0.685; % W/m/K
        case 'He'
            % Blais and Mann, "Thermal Conductivity of Helium and Hydrogen
            % at High Temperatures," J. Chem. Phys. 32(5), 1459-1465 (1960)
            conductivity_i = (991+0.678*(temperature-1200))*1e-6; % cal/cm/s/K
            conductivity_i = conductivity_i*4.184/0.01; % W/m/K
        case 'Kr'
            % Saxena and Saxena, "Thermal Conductivity of Krypton and Xenon
            % in the Temperature Range 350--1500°K," J. Chem. Phys. 51(8), 
            % 3361-3368 (1969)
            conductivity_i = (0.476 + 0.634e-2*temperature - 0.888e-6*temperature^2)/1e5; % cal/cm/s/K
            conductivity_i = conductivity_i*4.184/0.01; % W/m/K
        case 'CH4'
            % Afshar, Cogley and Saxena, "Thermal Conductivity of Methane 
            % at Atmospheric Pressure in the Temperature Range of 360--1275
            % K," J. Heat Transf. 102(1), 163-167 (1980)
            conductivity_i = (-23.35 + 0.1698*temperature + 1.893e-5*temperature^2)*1e-3; % W/m/K
        case 'N2O'
            % Saxena and Gupta, "Thermal conductivity of nitrous oxide in 
            % the temperature range 50° to 900° C," Chem. Phys. Lett. 4(5), 
            % 291-294 (1969)
            conductivity_i = (-2.674 + 2.37e-2*temperature - 5.454e-6*temperature^2)/1e5; % cal/cm/s/K
            conductivity_i = conductivity_i*4.184/0.01; % W/m/K
        case 'CO2'
            % The Engineering ToolBox (2018). Carbon dioxide - Thermal 
            % Conductivity vs. Temperature and Pressure.
            % [online] Available at:
            % https://www.engineeringtoolbox.com/carbon-dioxide-thermal-conductivity-temperature-pressure-d_2019.html
            T = [220,240,250,260,280,300,320,340,360,400,450,500,600,650,850,1050]'; % K
            kappa = [10.90,12.24,12.95,13.68,15.20,16.79,18.42,20.09,21.77,25.14,29.35,33.49,41.55,45.47,60.30,73.84]'*1e-3; % W/m/K
            % Fit to a cubic polynomial
            p = polyfit(T,kappa,3);
            conductivity_i = polyval(p,temperature);
    end
    conductivity = conductivity + conductivity_i*pressure_ratio(gas_i);
end

end

