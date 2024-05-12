% This code demonstrates the dispersive-wave generation in Fig. 5 of the
% following paper:
%
% Tani et al., "Multimode ultrafast nonlinear optics in optical waveguides:
% numerical modeling and experiments in kagome photonic-crystal fiber," J.
% Opt. Soc. Am. B 31 (2), 311-320 (2014).
%
% Due to the phase-matching relation, dispersive waves appear at different
% frequencies for different spatial modes.
%
% It uses a Xe-filled Kagome fiber.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.10,4]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'DW';
sim.midx = 1:4; % include four HE01,02,03,04 modes

num_save = 30;
fiber.L0 = 0.2; % m; propagation length
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
gas.core_radius = 13.5e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 2.7*1e5; % Pa
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 800e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'Xe';
gas.delta = 300e-9; % m; wall thickness of Kagome fibers
gas.fiber_type = 'Kagome';
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
tfwhm = 0.04; % ps

total_energy = 0.7e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,4,Nt,{'ifft',freq_shift},[1,1e-2,1e-2,1e-2]);

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Plot
spectrum = sum(abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2,2);
spectrum = squeeze(spectrum).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(f/1e3,prop_output.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); colormap(jet);
caxis([-60,0]); xlim([0.1,1.5]);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
ylabel('Propagation distance (cm)');
title('Spectrum');

% Draw white dashed lines
hold on;
x1 = [0.9,0.9,1.03,1.03,0.9]';
y1 = [10,20,20,10,10]';
plot(x1,y1,'--','Color','w','Linewidth',2);
text(0.9,8,'HE_{11}','Color','w','fontsize',12);
x2 = [1.095,1.095,1.224,1.224,1.095]';
y2 = [10,20,20,10,10]';
plot(x2,y2,'--','Color','w','Linewidth',2);
text(1.095,8,'HE_{12}','Color','w','fontsize',12);
x3 = [1.3,1.3,1.374,1.374,1.3]';
y3 = [10,20,20,10,10]';
plot(x3,y3,'--','Color','w','Linewidth',2);
text(1.3,8,'HE_{13}','Color','w','fontsize',12);
hold off;