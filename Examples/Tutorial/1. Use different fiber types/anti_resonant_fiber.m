% This code simulates soliton self-frequency shift (SSFS) in a 30-um-core
% anti-resonant hollow-core fiber filled with 15-bar H2.
% The input pulse is 50 fs at 1030 nm.
% This is a scalar simulation and doesn't include polarization modes. The
% input pulse is linearly polarized.

close all;  clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^13;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.ellipticity = 0; % linearly polarized for 0 and circularly polarized for 1

num_save = 100;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
% Load default parameters like 
%
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_deltaZ.threshold = 1e-8; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.Raman_model = 1; Use the Raman model (exist only in Raman-active gases)
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.num_noise_photon_per_bin = 1; Include photon shot noise
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'gas_UPPE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

% Configure gas parameters for the gas_info().
% These parameters vary based on different fiber type to use.
% Please check each example of each fiber for details regarding what
% parameters are required.
gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 15*1.01325e5; % Pa
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1030e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF'; % anti-resonant fiber
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.gas_material).(Raman_type).(Raman_parameters)
% 
% e.g.
%    gas.H2.R.T1 - H2's upper population lifetime for the excited rotational level (unused)
%    gas.H2.R.T2 - H2's coherence decay time between the ground state and the excited rotational level
%    gas.H2.R.gamma - H2's polarizability anisotropy
%    gas.H2.R.polarization_calibration - calibration factor to the paper's polarizability anisotropy value (gamma)
%    gas.H2.R.B0 - H2's B0 value for rotational Raman scattering
%    gas.H2.R.D0 - H2's D0 value for rotational Raman scattering
%    gas.H2.R.alpha_e - H2's alpha_e value for rotational Raman scattering
%    gas.H2.R.omega - H2's rotational transition angular frequencies
%    gas.H2.R.preR - H2's prefactors, representing Raman strength, of rotational Raman scattering
%
%    gas.H2.V.T1 - H2's upper population lifetime for the excited vibrational level (unused)
%    gas.H2.V.T2 - H2's coherence decay time between the ground state and the excited vibrational level
%    gas.H2.V.Dalpha - H2's polarizability derivative
%    gas.H2.V.Dgamma - H2's polarizability-anisotropy derivative
%    gas.H2.V.polarization_calibration - calibration factor to the paper's polarizability-derivative value (Dalpha,Dgamma)
%    gas.H2.V.Omega - H2's vibrational transition frequency for the ground rotational level Q(0), v = 0 --> 1, J = 0 --> 0 (unchanged J)
%    gas.H2.V.omega - H2's vibrational transition angular frequencies (calibrated with rotational transitions from Omega)
%    gas.H2.V.preR - H2's prefactors, representing Raman strength, of vibrational Raman scattering
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition
tfwhm = 0.050; % ps
total_energy = 1e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Plot
freq_range = [166,400];
wave_range = 3e5./fliplr(freq_range);

% Spectrum 
spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength = squeeze(spectrum_wavelength(:,1,:)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength); log_spectrum_wavelength1 = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
%spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
%spectrum = squeeze(spectrum(:,1,:)).';
%log_spectrum = 10*log10(spectrum); log_spectrum1 = log_spectrum - max(log_spectrum(:));

% Spectral evolution
figure;
pcolor(lambda,prop_output.z,log_spectrum_wavelength1); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim(wave_range);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (m)');
title('Spectral evolution');

% Output spectra
figure;
h = plot(lambda/1e3,spectrum_wavelength(end,:).'.*time_window^2/1e6);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('Spectrum (\muJ/\mum)');
xlim(wave_range/1e3);
title('Spectra');