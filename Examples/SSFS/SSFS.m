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
sim.ellipticity = 0; % linearly polarized for 0 and circularly polarized for 1

num_save = 100;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % um
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

% This is used to find the pump pulse time delay without any nonlinearity
% to correctly demonstrate how the output pulses are delayed w.r.t. the
% pump pulse.
fiber_calc_pulse_speed = fiber;
fiber_calc_pulse_speed.SR = 1e-30;
output_calc_pulse_speed = UPPE_propagate(fiber_calc_pulse_speed,initial_condition,sim,gas);

%% Save
save('SSFS_50fs_pL_H2.mat');