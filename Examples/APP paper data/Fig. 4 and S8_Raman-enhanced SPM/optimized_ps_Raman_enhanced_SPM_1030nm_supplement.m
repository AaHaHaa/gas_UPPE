% This code simulates the Raman-enhanced SPM effect with a ps
% transform-limited pulse, aiming to demonstrate the strongest pulse
% compression effect in a transient Stokes generation process.
% Please see Fig. S8(c) and (d) in the supplement for details.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,3]*1e-6; % m
Nt = 2^12;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SPM';
sim.gpuDevice.Index = 1;
sim.pulse_centering = false;

num_save = 1;
fiber.L0 = 0.5; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 150e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 20*1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 800e-9; % m
gas.material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 1;

total_energy = 2000e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);


%%
save('optimized_ps_Raman_enhanced_SPM_1030nm.mat');