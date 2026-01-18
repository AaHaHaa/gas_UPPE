function [fiber,sim,gas] = gas_info_gradient_pressure(fiber,sim,gas,wavelength)
%GAS_INFO_GRADIENT_PRESSURE It loads parameters related to gradient-pressure gases
%   wavelength: m
%
% The gradient pressure p(x) = sqrt( p_in^2 + z/L*(p_out^2 - p_in^2) )

%% Error check
if length(gas.pressure_in) ~= length(gas.material) || length(gas.pressure_out) ~= length(gas.material)
    error('gas_info_gradient_pressure:Pressure_MaterialError',...
          'The length of the gas-pressure (numeric) array should be the same as the length of the gas-material (cell) array');
end

num_gas = length(gas.material);

if any(ismember(gas.material,'air'))
    if any(ismember(gas.material,{'N2','O2'}))
        error('gas_info_gradient_pressure:GasMaterialError',...
              ['If gas.material contains ''air'', other materials cannot be ''N2'' or ''O2''.\n',...
               'If necessary, just specify N2 and O2 independently, rather than using ''air''.']);
    end
end

if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('gas_info_gradient_pressure:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

if sum(gas.pressure_in) == 0
    if sum(gas.pressure_out) == 0 % no gas and thus no gradient
        gas.pressure_ratio = ones(1,length(gas.pressure_in))/length(gas.pressure_in); % just make it all equal ratio
    else
        gas.pressure_ratio = gas.pressure_out/sum(gas.pressure_out);
    end
else
    gas.pressure_ratio = gas.pressure_in/sum(gas.pressure_in);

    if sum(gas.pressure_out) ~= 0
        pressure_ratio_out = gas.pressure_out/sum(gas.pressure_out);
        if any((gas.pressure_ratio - pressure_ratio_out) > eps(1))
            error('gas_info_gradient_pressure:PressureError',...
                  'The input and output pressure content should have the same ratio of gases.');
        end
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
addpath([upper_folder 'HCF/helper functions'],[upper_folder 'Gas absorption spectra']);

permittivity0 = 8.8541878176e-12; % F/m
c = 299792458; % m/s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius

%% Raman parameters

% Number density of the gas
% The final gas density is gas.Ng.*eta.
% eta is included during the propagation.
gas.Ng = gas.pressure_ratio*pressure0/k/temperature0; % m^(-3)

% eta is calculated with the unit, amagat
% Ideal gas law is used here.
% Because lower pressure creates a longer dephasing time, I take the
% close-to-minimum pressure below to find the corresponding dephasing time,
% T2. The reason why I don't pick the minimum pressure is to avoid zero
% pressure where gas is neglected.
% This is important to determine which propagation model to use.
% In my UPPE, I have two Raman models.
% One is to use aperiodic convolution with a small time window, which
% is much smaller than the dephasing time of the Raman response; the
% other is to use a normal periodic convolution when the time window is
% much larger than the dephasing time.
eta_min = gas.pressure_ratio*min(sum(gas.pressure_in),sum(gas.pressure_out))/pressure0*temperature0/gas.temperature;
eta_max = gas.pressure_ratio*max(sum(gas.pressure_in),sum(gas.pressure_out))/pressure0*temperature0/gas.temperature;
if any(ismember(gas.material,{'H2','D2','N2','O2','air','CH4','N2O','CO2'}))
    gas = Raman_model(gas,eta_min + 0.1*(eta_max-eta_min)); % obtain the Raman parameters according to the gas; 0.1 is to take the close-to-minimum pressure
else % no Raman
    sim.include_Raman = false;
end

%% Refractive index at 1 amagat
Dpermittivity_r = 0; % initialization
Raman_absorption = 0; % initialization
for gas_i = 1:num_gas
    [n_gas_i,Raman_absorption_i] = find_n_gas(gas.material{gas_i},wavelength,gas.pressure_ratio(gas_i));
    Dpermittivity_r_i = real(n_gas_i).^2-1; % difference from 1 of each gas

    Dpermittivity_r = Dpermittivity_r + Dpermittivity_r_i;
    Raman_absorption = Raman_absorption + Raman_absorption_i;
end

gas.permittivity_r = Dpermittivity_r + 1;
gas.Raman_absorption = Raman_absorption./(2*pi./wavelength);

%% Propatation constant (updated with Raman parameters)
% Because of the gradient-pressure distribution, the system pressure has a
% mean that lies toward the maximum end.
eta_mean = eta_min + 0.6*(eta_max-eta_min);
switch gas.fiber_type
    case 'Ag_coating' % Typical capillary has an Ag coating inside the capillary to prevent the loss.
                      % Silica alone generates extremely lossy guided modes.
                      % Ag coating provides excellent guiding at visible/NIR.
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Ag_coating_beta_func(wavelength,eta_mean,sim,gas);
    case 'no_coating' % No coating: silica alone
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_no_coating_beta_func(wavelength,eta_mean,sim,gas);
    case 'MWLW_coating' % MWLW coating
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_MWLW_coating_beta_func(wavelength,eta_mean,sim,gas);
    case 'AR_HC_PCF' % anti-resonant hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_AR_HC_PCF_beta_func(wavelength,eta_mean,sim,gas);
    case 'Kagome' % Kagome hollow-core fiber
        [fiber.betas,fiber.SR,sim.mode_profiles] = solve_for_EH_Kagome_beta_func(wavelength,eta_mean,sim,gas);
    otherwise
        error('gas_info:FiberTypeError',...
              'fiber_type is wrong.');
end
n_eff = real(fiber.betas./(2*pi./wavelength));

gas.beta_mean_gas = fiber.betas;
gas.n_gas_mean = sqrt((gas.permittivity_r - 1)*sum(eta_mean) + 1); % refractive index of the gas at eta_mean
gas.k0 = 2*pi./wavelength; % 1/m

%% Nonlinear coefficient
n2 = gas_n2(gas.material,wavelength,gas.pressure_ratio);
fiber.X3_prefactor = n2/3*4*permittivity0.*n_eff.^2*c; % m^2/V^2

%% Ionization potential
if sim.include_photoionization
    gas = gas_photoionization_parameters(gas);
end

%% uint32 is required in the cuda computation later
gas.wavelength_order = uint32(gas.wavelength_order);

end