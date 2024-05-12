% This code tends to simulate soliton self-frequency shift (SSFS) in a 
% 30-um-core anti-resonant hollow-core fiber filled with 15-bar H2.
% The input pulse is 50 fs at 1030 nm.
% The simulation includes polarization modes where the input pulse is
% circularly polarized.
%
% Compared to linearly polarized simulations in this folder, SSFS doesn't
% happen because strong polarization coupling distorts the SSFS process.

close all;  clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.3,3]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.scalar = false; % run with polarized fields
sim.ellipticity = 1; % linearly polarized for 0 and circularly polarized for 1

num_save = 100;
fiber.L0 = 0.2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 15*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.050; % ps
total_energy = 1e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});
initial_condition.fields = initial_condition.fields.*sqrt([1,1e-2]);

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Save
save('SSFS_50fs_pC_H2.mat');