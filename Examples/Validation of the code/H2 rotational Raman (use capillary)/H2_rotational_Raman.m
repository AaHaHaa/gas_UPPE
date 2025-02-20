% This code duplicates the simulations of H2 under a 1540-nm pulse in 
% Fig. 14(c) in "Pure rotational stimulated Raman scattering in H2-filled 
% hollow-core photonic crystal fibers (2020)"

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [1.38,1.75]*1e-6; % m
Nt = 2^22;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);
sim.f0 = f0;
sim.progress_bar_name = 'H2 S(1)';

num_save = 1;
fiber.L0 = 10; % m; propagation length
sim.save_period = 0;%fiber.L0/num_save;

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
% sim.adaptive_dz.threshold = 1e-8; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = true; Consider the Raman (exist only in Raman-active gases)
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
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

gas.core_radius = 4.5e-6/0.64; % m; their fiber has 9-um MFD, so I need to back-calculate the actual core radius.
gas.temperature = 300; % K
gas.pressure = 16*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1540e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = 'H2';
gas.fiber_type = 'no_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
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

% Set up the loss to be consistent with the paper
fiber.betas(lambda<1600) = real(fiber.betas(lambda<1600)) + 1i*log(10^(0.04/10/2));
fiber.betas(lambda>1600) = real(fiber.betas(lambda>1600)) + 1i*log(10^(0.11/10/2));

% Consider only the S(0) and S(1)
gas.H2.R.preR = [gas.H2.R.preR([1,2]),zeros(1,length(gas.H2.R.preR)-2)];
gas.H2.V.preR = zeros(size(gas.H2.V.preR)); % no vib

%% Initial condition
tfwhm = 12e3; % ps
total_energy = 1000; % nJ
pump_wavelength = 1540e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagation
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Plot
filter_lambda = 1600;
probe = struct('dt',prop_output.dt,'fields',prop_output.fields(:,:,end));
residual = edgepass_filter('lowpass', probe, sim.f0, filter_lambda,0.3);
Raman = edgepass_filter('highpass', probe, sim.f0, filter_lambda,0.3);
residual_energy = trapz(t,abs(residual.fields).^2)/1e3; % nJ
Raman_energy = trapz(t,abs(Raman.fields).^2)/1e3; % nJ

spectrum = abs(fftshift(ifft(probe.fields(:,1,:)),1)).^2*time_window^2/1e6.*3e5./(3e5./f).^2; % uJ/nm
spectrum = squeeze(spectrum(:,1,:)).';

final_spectrum = spectrum(end,:)';
final_spectrum = final_spectrum/max(final_spectrum);
figure;
plot(3e5./f,final_spectrum,'linewidth',2);
xlim([1500,1750]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');

%% Save the results
save(sprintf('S1_Raman_gain_%u.mat',total_energy),'-v7.3','t','filter_lambda','lambda','probe','pump_wavelength','Raman_energy','residual_energy','tfwhm','time_window','total_energy');