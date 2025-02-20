% This code duplicates the simulations of Xe in Fig. 4 in Tani's "Multimode
% ultrafast nonlinear optics in optical waveguides: numerical modeling and 
% experiments in kagom� photonic-crystal fiber"

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.15,5]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.midx = [1,2]; % consider two spatial modes, EH11 and EH21

num_save = 100;
fiber.L0 = 0.3; % m; propagation length
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
gas.core_radius = 9e-6; % m
gas.temperature = 288.15; % K
gas.pressure = 25*1.013e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 800e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = 'Xe';
gas.fiber_type = 'no_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

fiber.betas = real(fiber.betas); % ignore the loss

%% Initial condition
tfwhm = 0.15; % ps
total_energy = 0.6e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,-time_window/12);

input_field.fields = [input_field.fields,input_field.fields/1e2];

%% Propagation
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Save the results
save('Xe.mat');