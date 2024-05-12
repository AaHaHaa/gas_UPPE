% This code demonstrates the dispersive-wave generation in Fig. 5 of the
% following paper:
%
% Tani et al., "Multimode ultrafast nonlinear optics in optical waveguides:
% numerical modeling and experiments in kagome photonic-crystal fiber," J.
% Opt. Soc. Am. B 31 (2), 311-320 (2014).
%
% Due to the phase-matching relation, dispersive waves appear at different
% frequencies for different spatial modes.

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
sim.midx = 1:4;

num_save = 30;
fiber.L0 = 0.2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 13.5e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 2.7*1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 800e-9; % m
gas.gas_material = 'Xe';
gas.delta = 300e-9; % m; wall thickness of Kagome fibers
gas.fiber_type = 'Kagome';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.04; % ps

total_energy = 0.7e3; % nJ
pump_wavelength = 800e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,4,Nt,{'ifft',freq_shift},[1,1e-2,1e-2,1e-2]);

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
save('DW.mat');