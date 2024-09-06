% This code demonstrates the self-phase modulation, which results in pulse 
% compression, of a pulse in a Ar-filled HCF.
%
% It aims to duplicate the result from the reference paper below and
% calibrates the code for Ar.
%
% Suda et al., "Generation of sub-10-fs, 5-mJ-optical pulses using a hollow
% fiber with a pressure gradient," Appl. Phys. Lett. 86, 111116 (2005)
%
% For Ar, the following two papers seem to require around 7 times weaker n2
% to match their experiments. We have done Ar experiments ourselves which
% show that this reducing factor isn't necessary. We suspect that these old
% works didn't do capillary experiments correctly. It's difficult to align
% the beam into the capillary, which I've been experiencing and is a huge 
% pain. Super straight capillary is required. Super high quality beam is
% also required. Without all these, higher-order modes arise, quickly
% deteriorating the beam spatial quality and the polarization extinction
% ratio. Strong polarization coupling in HE modes in a capillary has
% significant influence to nonlinear propagation.
% 1. Sartania et al.,
%    "Generation of 0.1-TW 5-fs optical pulses at a 1-kHz repetition rate," Opt. Lett. 22(20), 1562-1564 (1997)
% 2. Suda et al.,
%    "Generation of sub-10-fs, 5-mJ-optical pulses using a hollow fiber with a pressure gradient," Appl. Phys. Lett. 86, 111116 (2005)

close all;  clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.4,2]*1e-6; % m
Nt = 2^10;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.gpu_yes = false;

num_save = 10;
fiber.L0 = 2.2; % m; propagation length
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
gas.core_radius = 250e-6; % m
gas.temperature = 288; % K
gas.pressure_in = 0; % Pa
gas.pressure_out = 40e3; % Pa
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 780e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'Ar';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.gas_material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info_gradient_pressure(fiber,sim,gas,lambda*1e-9);

sim.X3 = sim.X3/7;

%% Initial condition
tfwhm = 0.040; % ps
total_energy = 6.5e6; % nJ
pump_wavelength = 780e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('t (ps)');
ylabel('Power (W)');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([-0.2,0.2]);

% Spectrum
figure;
c = 299792.458; % nm/ps
wavelength = c./f; % nm
factor_correct_unit = time_window^2/1e6; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e6" is to make pJ into uJ
factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2*factor_correct_unit.*factor;
h = plot(lambda,spectrum_wavelength(:,:,end));
xlabel('Wavelength (nm)');
ylabel('PSD (\muJ/nm)');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

%% Save the data
%save('SPM.mat');