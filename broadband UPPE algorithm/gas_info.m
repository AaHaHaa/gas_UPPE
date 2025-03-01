function [fiber,sim,gas] = gas_info(fiber,sim,gas,wavelength)
%GAS_INFO It loads parameters related to gases
%   wavelength: m

%% Preset
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));
addpath([upper_folder 'HCF/helper functions']);

permittivity0 = 8.8541878176e-12; % F/m
c = 299792458; % m/s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius
% eta is calculated with the unit, amagat
% Ideal gas law is used here.
eta = gas.pressure/pressure0*temperature0/gas.temperature;
if isfield(gas,'mixed_material') && ~isempty(gas.mixed_material) % mixed gases
    mixed_eta = zeros(1,length(gas.mixed_material));
    for mat_i = 1:length(gas.mixed_material)
        mixed_eta(mat_i) = gas.mixed_pressure/pressure0*temperature0/gas.temperature;
    end
else
    mixed_eta = [];
end

%% Raman parameters

% Number density of the gas
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

if ismember(gas.material,{'H2','N2','O2','air','CH4','N2O','CO2'})
    gas = Raman_model(gas,eta); % Obtain the Raman parameters according to the gas
else % no Raman
    sim.include_Raman = false;
end

%% Propatation constant (updated with Raman parameters)
switch gas.fiber_type
    case 'Ag_coating' % Typical capillary has an Ag coating inside the capillary to prevent the loss.
                      % Silica alone generates extremely lossy guided modes.
                      % Ag coating provides excellent guiding at visible/NIR.
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Ag_coating_beta_func(wavelength,eta,sim,gas,mixed_eta);
    case 'no_coating' % No coating: silica alone
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_no_coating_beta_func(wavelength,eta,sim,gas,mixed_eta);
    case 'MWLW_coating' % MWLW coating
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_MWLW_coating_beta_func(wavelength,eta,sim,gas,mixed_eta);
    case 'AR_HC_PCF' % anti-resonant hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_AR_HC_PCF_beta_func(wavelength,eta,sim,gas,mixed_eta);
    case 'Kagome' % Kagome hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Kagome_beta_func(wavelength,eta,sim,gas,mixed_eta);
    otherwise
        error('fiber_type is wrong.');
end
n_eff = real(fiber.betas./(2*pi./wavelength));

%% Nonlinear coefficient
n2 = gas_n2(gas.material,wavelength);
fiber.X3 = (n2/3*4*permittivity0.*n_eff.^2*c)*eta; % m^2/V^2

if isfield(gas,'mixed_material') && ~isempty(gas.mixed_material) % mixed gases
    for mat_i = 1:length(gas.mixed_material)
        mixed_n2 = gas_n2(gas.mixed_material{mat_i},wavelength);
        fiber.X3 = fiber.X3 + (mixed_n2/3*4*permittivity0.*n_eff.^2*c)*mixed_eta(mat_i); % m^2/V^2
    end
end

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