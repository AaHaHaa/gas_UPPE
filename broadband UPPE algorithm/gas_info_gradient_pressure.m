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
addpath([upper_folder 'HCF/helper functions'],[upper_folder 'Gas absorption spectra']);

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

if ismember(gas.material,{'H2','N2','O2','air','N2O','CO2'})
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
    sim.include_Raman = false;
end

%% Refractive index

% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)

[a,b] = Sellmeier_coefficients(gas.material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
switch gas.material
    case 'H2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        gas.permittivity_r_120 = n_from_Sellmeier(0.120).^2;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,1);
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        gas.permittivity_r_400 = n_from_Sellmeier(0.4).^2;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,1);
    case {'air','N2','N2O'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,1);
    case {'Ar','CO2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;

        Raman_absorption = 0;
    case {'He','Ne'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

        Raman_absorption = 0;
    case {'Kr','Xe'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        gas.permittivity_r_113 = n_from_Sellmeier(0.113).^2;

        Raman_absorption = 0;
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

        Raman_absorption = 0;
end
gas.permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
gas.Raman_absorption = Raman_absorption./(2*pi./wavelength);

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
n2 = gas_n2(gas.material,wavelength);
fiber.X3_prefactor = n2/3*4*permittivity0.*n_eff.^2*c; % m^2/V^2

%% Ionization potential
if sim.photoionization_model ~= 0
    [ionization_energy,l,Z] = gas_photoionization_parameters(gas.material);
    e = 1.60217663e-19; % Coulomb
    gas.ionization = struct('energy',ionization_energy*e,... % J
                            'l',l,... % quantum number l
                            'Z',Z); % effective charge
end

%% uint32 is required in the cuda computation later
gas.wavelength_order = uint32(gas.wavelength_order);

end