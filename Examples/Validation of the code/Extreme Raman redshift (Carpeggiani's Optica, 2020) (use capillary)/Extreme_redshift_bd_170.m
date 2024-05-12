% This code aims to duplicate Fig. 1(b) and (d) in the Raman redshifting
% paper by Carpeggiani et al.
% It's used to calibrate our code with N2.
%
% I assume that their pulse might have strong ASE. In their supplement,
% they mentioned that it doesn't contribute to nonlinear interactions.
%
% Carpeggiani et al., "Extreme Raman red shift: ultrafast multimode
% nonlinear space-time dynamics, pulse compression, and broadly tunable
% frequency conversion," Optica 71349-1354 (2020)

close all;  clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,3]*1e-6; % m
Nt = 2^12;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.midx = 1:5;

num_save = 30;
fiber.L0 = 6; % m; propagation length
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
gas.core_radius = 265e-6; % m
gas.temperature = 300; % K
gas.pressure = 4*1e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1030e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'N2';
gas.fiber_type = 'no_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.gas_material).(Raman_type).(Raman_parameters)
% 
% e.g.
%    gas.N2.R.T2 - N2's coherence decay time between the ground state and the excited rotational level
%    gas.N2.R.gamma - N2's polarizability anisotropy
%    gas.N2.R.polarization_calibration - calibration factor to the paper's polarizability anisotropy value (gamma)
%    gas.N2.R.B0 - N2's B0 value for rotational Raman scattering
%    gas.N2.R.D0 - N2's D0 value for rotational Raman scattering
%    gas.N2.R.alpha_e - N2's alpha_e value for rotational Raman scattering
%    gas.N2.R.omega - N2's rotational transition angular frequencies
%    gas.N2.R.preR - N2's prefactors, representing Raman strength, of rotational Raman scattering
%
%    gas.N2.V.T2 - N2's coherence decay time between the ground state and the excited vibrational level
%    gas.N2.V.Dalpha - N2's polarizability derivative
%    gas.N2.V.Dgamma - N2's polarizability-anisotropy derivative
%    gas.N2.V.polarization_calibration - calibration factor to the paper's polarizability-derivative value (Dalpha,Dgamma)
%    gas.N2.V.Omega - N2's vibrational transition frequency for the ground rotational level Q(0), v = 0 --> 1, J = 0 --> 0 (unchanged J)
%    gas.N2.V.omega - N2's vibrational transition angular frequencies (calibrated with rotational transitions from Omega)
%    gas.N2.V.preR - N2's prefactors, representing Raman strength, of vibrational Raman scattering
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition
tfwhm = 0.170;
total_energy = 0.7e6*1.08;
freq_shift = 299792.458/1030 - sim.f0;
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt, {'ifft',freq_shift}, sqrt([1,1e-2*ones(1,length(sim.midx)-1)]));

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Calibrated pulse
% This is used to find the pump pulse time delay without any nonlinearity
% to correctly demonstrate how the output pulses are delayed w.r.t. the
% pump pulse.
fiber_calc_pulse_speed = fiber;
fiber_calc_pulse_speed.SR = fiber.SR*1e-30;
output_calc_pulse_speed = UPPE_propagate(fiber_calc_pulse_speed,initial_condition,sim,gas);

%% Save the results
save(sprintf('redshift_%s_%.1fm_%uum_%2.1fbar.mat',gas.gas_material,fiber.L0,gas.core_radius*2*1e6,gas.pressure/1e5),'-v7.3',...
     'c','dt','f','f_range','fiber','gas','lambda','Nt','num_save','prop_output','sim','t','time_window','wavelength_range','output_calc_pulse_speed');