% This code duplicates the results of double-pulse Raman generation done by
% Konyashchenko et al.
%
% The technique of Raman conversion of sub-100-fs laser pulses based on 
% excitation of active medium by two orthogonally polarized pulses has been
% developed for Raman lasers with a glass capillary. 52 fs Stokes pulse at 
% the wavelength of 1200 nm has been generated by stimulated Raman 
% scattering of 48 fs Ti:sapphire laser pulse at the wavelength of 800 nm 
% in hydrogen. 13% energy conversion efficiency has been achieved at pulse 
% repetition rate up to 2 kHz.
%
% Konyashchenko et al., "Frequency shifting of sub-100 fs laser pulses by 
% stimulated Raman scattering in a capillary filled with pressurized gas," 
% Appl. Phys. B 93, 455-461 (2008)

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.3,10]*1e-6; % m
Nt = 2^17;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = '130-bar Raman from 800nm';
sim.pulse_centering = false;

num_save = 30;
fiber.L0 = 0.1; % m; propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

coupling_efficiency = 0.72; % coupling efficiency of the input pulse into the capillary

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
gas.core_radius = 50e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 130*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 800e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
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

%% Single-pulse pumping
tfwhm = 0.048; % ps
total_energy = 70e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1);
tfwhm2 = 1.5; % ps
func = calc_chirp;
[~,chirped_field] = func.Gaussian( tfwhm2,ifftshift(2*pi*f,1),ifft(input_field.fields),1 );
input_field.fields = chirped_field*sqrt(coupling_efficiency);

prop_output{1} = UPPE_propagate(fiber,input_field,sim,gas);

%% Two-pulse pumping
tfwhm = 0.048; % ps
total_energy = 45e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field(1) = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,-5);
tfwhm2 = 1.5; % ps
func = calc_chirp;
[~,chirped_field] = func.Gaussian( tfwhm2,ifftshift(2*pi*f,1),ifft(input_field(1).fields),1 );
input_field.fields = chirped_field;

total_energy = 30e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field(2) = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,5);
tfwhm2 = 1.5; % ps
func = calc_chirp;
[~,chirped_field] = func.Gaussian( tfwhm2,ifftshift(2*pi*f,1),ifft(input_field(2).fields),1 );
input_field(2).fields = chirped_field;

initial_condition.fields = (input_field(1).fields + input_field(2).fields)*sqrt(coupling_efficiency);
initial_condition.dt = dt;
prop_output{2} = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Save the results
save('Raman_52fs.mat');