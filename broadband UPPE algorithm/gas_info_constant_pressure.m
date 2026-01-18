function [fiber,sim,gas] = gas_info_constant_pressure(fiber,sim,gas,wavelength)
%GAS_INFO_CONSTANT_PRESSURE It loads parameters related to constant-pressure gases
%   wavelength: m

%% Error check
if length(gas.pressure) ~= length(gas.material)
    error('gas_info_constant_pressure:Pressure_MaterialError',...
          'The length of the gas-pressure (numeric) array should be the same as the length of the gas-material (cell) array');
end
num_gas = length(gas.material);

if any(ismember(gas.material,'air'))
    if any(ismember(gas.material,{'N2','O2'}))
        error('gas_info_constant_pressure:GasMaterialError',...
              ['If gas.material contains ''air'', other materials cannot be ''N2'' or ''O2''.\n',...
               'If necessary, just specify N2 and O2 independently, rather than using ''air''.']);
    end
end

if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('gas_info_constant_pressure:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

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
eta = zeros(1,num_gas);
for gas_i = 1:num_gas
    eta(gas_i) = gas.pressure(gas_i)/pressure0*temperature0/gas.temperature;
end

%% Raman parameters

% Number density of the gases
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

if any(ismember(gas.material,{'H2','D2','N2','O2','air','CH4','N2O','CO2'}))
    gas = Raman_model(gas,eta); % obtain the Raman parameters according to the gas
else % no Raman
    sim.include_Raman = false;
end

%% Propatation constant (updated with Raman parameters)
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
        error('gas_info_constant_pressure:FiberTypeError',...
              'fiber_type is wrong.');
end
n_eff = real(fiber.betas./(2*pi./wavelength));

%% Nonlinear coefficient
n2 = gas_n2(gas.material,wavelength,eta);
fiber.X3 = n2/3*4*permittivity0.*n_eff.^2*c; % m^2/V^2

%% Ionization potential
if sim.include_photoionization
    gas = gas_photoionization_parameters(gas);
end

%% uint32 is required in the cuda computation later
gas.wavelength_order = uint32(gas.wavelength_order);

end