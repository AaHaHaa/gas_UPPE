% This code duplicates Huang's "Continuously wavelength-tunable 
% blueshifting soliton generated in gas-filled photonic crystal fibers," 
% Opt. Lett. 7, 1805-1808 (2019)

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.2,5]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'photoionization';
sim.pulse_centering = false;
sim.photoionization_model = true;
%sim.gpu_yes = false;

num_save = 100;
fiber.L0 = 0.2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 12.3e-6; % m
gas.temperature = 288; % K
gas.pressure = 10*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 800e-9; % m
gas.gas_material = 'He';
gas.fiber_type = 'AR_HC_PCF';
gas.delta = 215e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 1e-4; % loss factor
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.02; % ps
total_energy = 3.2e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagate
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
save(sprintf('photoionization_blueshift_%2.1fuJ.mat',total_energy/1e3));