% This code simulates the photonionization blueshift in a 10-atm He-filled
% anti-resonant fiber.
%
% Constant pressure is applied. It undergoes more blueshift than in a
% gradient pressure.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.2,5]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.gpuDevice.Index = 1;
sim.pulse_centering = false;
sim.progress_bar_name = 'constant pressure';
sim.photoionization_model = true; % activate the photoionization effect

num_save = 1;
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
gas.core_radius = 13e-6; % m
gas.temperature = 300; % K
gas.pressure = 10*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 532e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = 'He';
gas.fiber_type = 'AR_HC_PCF';
gas.num_tubes = 7; % the number of tubes in the anti-resonant fiber
gas.r_tube = 12.5e-6; % m; the tube radius (not core!)
gas.t_tube = 215e-9; % m; the tube's wall thickness of anti-resonant fibers
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition
tfwhm = 0.02; % ps
total_energy = 6e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagation
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Plot
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum = squeeze(spectrum);

figure;
h = plot(3e5./f,spectrum(:,end));
xlim([400,1200]);
set(h,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (dB)');