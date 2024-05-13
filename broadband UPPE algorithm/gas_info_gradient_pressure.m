function [fiber,sim,gas] = gas_info_gradient_pressure(fiber,sim,gas,wavelength)
%GAS_INFO Summary of this function goes here
%   wavelength: m
%
% The gradient pressure p(x) = sqrt( p_in^2 + z/L*(p_out^2 - p_in^2) )

%% Preset
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));
addpath([upper_folder 'HCF/helper functions'],[upper_folder 'gas absorption spectra']);

permittivity0 = 8.8541878176e-12; % F/m
c = 299792458; % m/s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius

%% Raman parameters

% Number density of the gas
% The final gas density is gas.Ng*eta.
% eta is included during the propagation.
gas.Ng = pressure0/k/temperature0; % m^(-3)

if ismember(gas.gas_material,{'H2','N2','O2','air'})
    % eta is calculated with the unit, amagat
    % Ideal gas law is used here.
    % Because lower pressure creates a longer dephasing time, I take the
    % minimum pressure below to find the corresponding dephasing time, T2.
    % This is important to determine which propagation model to use.
    % In my UPPE, I have two Raman models.
    % One is to use aperiodic convolution with a small time window, which
    % is much smaller than the dephasing time of the Raman response; the
    % other is to use a normal periodic convolution when the time window is
    % much larger than the dephasing time.
    eta = min(gas.pressure_in,gas.pressure_out)/pressure0*temperature0/gas.temperature;
    
    gas = Raman_model(gas,eta); % Obtain the Raman parameters according to the gas
else % no Raman
    sim.Raman_model = 0;
end

%% Refractive index

% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)

[a,b] = Sellmeier_coefficients(gas.gas_material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
switch gas.gas_material
    case 'H2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        gas.permittivity_r_120 = n_from_Sellmeier(0.120).^2;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,wavelength,1);
        gas.Raman_absorption = Raman_absorption./(2*pi./wavelength);
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        gas.permittivity_r_400 = n_from_Sellmeier(0.4).^2;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,wavelength,1);
        gas.Raman_absorption = Raman_absorption./(2*pi./wavelength);
    case {'air','N2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,wavelength,1);
        gas.Raman_absorption = Raman_absorption./(2*pi./wavelength);
    case {'He','Ne','Ar'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
    case {'Kr','Xe'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        gas.permittivity_r_113 = n_from_Sellmeier(0.113).^2;
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
end

%% Propatation constant (updated with Raman parameters)
% eta is calculated with the unit, amagat
% Ideal gas law is used here.
%pressure = sqrt((gas.pressure_in^2 + gas.pressure_out^2)/2); % average pressure is calculated by taking sqrt(integral(p^2,x)/L)
%eta = pressure/pressure0*temperature0/gas.temperature;
eta = 0;

switch gas.fiber_type
    case 'Ag_coating' % Typical capillary has an Ag coating inside the capillary to prevent the loss.
                      % Silica alone generates extremely lossy guided modes.
                      % Ag coating provides excellent guiding at visible/NIR.
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Ag_coating_beta_func(wavelength,eta,sim,gas);
    case 'no_coating' % No coating: silica alone
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_no_coating_beta_func(wavelength,eta,sim,gas);
    case 'MWLW_coating' % MWLW coating
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_MWLW_coating_beta_func(wavelength,eta,sim,gas);
    case 'AR_HC_PCF' % anti-resonant hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_AR_HC_PCF_beta_func(wavelength,eta,sim,gas);
    case 'Kagome' % Kagome hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Kagome_beta_func(wavelength,eta,sim,gas);
    otherwise
        error('fiber_type is wrong.');
end
n_eff = real(fiber.betas./(2*pi./wavelength));

gas.beta_no_gas = fiber.betas;
gas.k0 = 2*pi./wavelength; % 1/m
gas.coeff_beta_with_eta = gas.k0.^2./gas.beta_no_gas;

%% Nonlinear coefficient
switch gas.gas_material
    case 'H2' % m^2/(W*atm)
              % This value is taken from Wahlstrand, et al., 
              % "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2" (2015)
              % Its value, after converted into X3, is close to the paper by Belli et al.,
              % "Vacuum-ultraviolet to infrared supercontinuum in hydrogen-filled photonic crystal fiber" Optica (2015)
              % with X3 = 2.206e-26 m^2/V^2 at standard conditions
        n2 = 0.65e-23;
    case 'N2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'O2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'air' % Calculate n2 for N2 and O2
               % Add them up with 79% N2 and 21% O2
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2_N2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_N2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2_O2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_O2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2 = 0.79*n2_N2 + 0.21*n2_O2; %clear n2_N2 n2_O2
        n2(isinf(n2)) = 0; % avoid the singularity
    case 'Xe' % m^2/(W*atm)
              % From Shu-Zee Alencious Lo, et al.,
              % "Pulse propagation in hollow-core fiber at high-pressure regime: application to compression of tens of ?J pulses and determination of nonlinear refractive index of xenon at 1.03um" Applied Optics (2018)
        n2 = 50.1e-24;
    case 'Ar' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 7.96e-24;
    case 'Ne' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.85e-24;
    case 'He' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.34e-24;
    case 'Kr' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 18.9e-24;
end
sim.X3_prefactor = n2/3*4*permittivity0.*n_eff.^2*c; % m^2/V^2

%% Ionization potential
if sim.photoionization_model ~= 0
    switch gas.gas_material % from "NIST Chemistry WebBook"
        case 'H2'
            ionization_energy = 15.42593; % eV
        case 'N2'
            ionization_energy = 15.581; % eV
        case 'O2'
            ionization_energy = 12.0697; % eV
        case 'CH4'
            ionization_energy = 12.61; % eV
        case 'He'
            ionization_energy = 24.58741; % eV
        case 'Ne'
            ionization_energy = 21.56454; % eV
        case 'Ar'
            ionization_energy = 15.759; % eV
        case 'Kr'
            ionization_energy = 13.99961; % eV
        case 'Xe'
            ionization_energy = 12.12987; % eV
        otherwise
            error('This code doesn''t support the ionization computation of other materials yet');
    end
    e = 1.60217663e-19; % Coulomb
    gas.ionization_energy = ionization_energy*e; % J
end

%% uint32 is required in the cuda computation later
gas.wavelength_order = uint32(gas.wavelength_order);

end