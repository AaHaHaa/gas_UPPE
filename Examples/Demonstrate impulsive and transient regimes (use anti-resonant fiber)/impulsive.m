% This code shows the Raman-induced index change in the impulsive regime.

close all; clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.12,20]*1e-6; % m
Nt = 2^9;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'impulsive';
sim.gpu_yes = false;

num_save = 100;
fiber.L0 = 0.05; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 12.5e-6; % m
gas.temperature = 288; % K
gas.pressure = 5*1.013e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 800e-9; % m
gas.gas_material = 'H2';
gas.delta = 250e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 1e-2; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%% Initial condition and Propagate
tfwhm = 0.010; % ps
total_energy = 1e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagate
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
save('impulsive.mat');