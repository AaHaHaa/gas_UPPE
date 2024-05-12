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

%% Raman parameters

% Number density of the gas
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

if ismember(gas.gas_material,{'H2','N2','O2','air','CH4'})
    gas = Raman_model(gas,eta); % Obtain the Raman parameters according to the gas
else % no Raman
    sim.Raman_model = 0;
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
        error('fiber_type is wrong.');
end
n_eff = real(fiber.betas./(2*pi./wavelength));

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
    case 'CH4'
        n2 = 3.118e-23; % m^2/(W*atm)
end
sim.X3 = (n2/3*4*permittivity0.*n_eff.^2*c)*eta; % m^2/V^2

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