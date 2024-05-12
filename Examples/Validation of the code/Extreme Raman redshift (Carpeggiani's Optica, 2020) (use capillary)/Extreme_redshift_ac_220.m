% This code aims to duplicate Fig. 1(a) and (c) in the Raman redshifting
% paper by Carpeggiani et al.
% It's used to calibrate our code with N2.
%
% I assume that their pulse might have strong ASE. In their supplement,
% they mentioned that it doesn't contribute to nonlinear interactions.
%
% Carpeggiani et al., "Extreme Raman red shift: ultrafast multimode
% nonlinear space-time dynamics, pulse compression, and broadly tunable
% frequency conversion," Optica 71349-1354 (2020)

close all;  clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,3]*1e-6; % m
Nt = 2^12;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.midx = 1:5;

num_save = 30;
fiber.L0 = 5.5; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 500e-6; % m
gas.temperature = 300; % K
gas.pressure = 0.9*1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'N2';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.220;
total_energy = 3e6*0.95;
freq_shift = 299792.458/1030 - sim.f0;
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt, {'ifft',freq_shift}, sqrt([1,1e-4*ones(1,length(sim.midx)-1)]));

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

% This is used to find the pump pulse time delay without any nonlinearity
% to correctly demonstrate how the output pulses are delayed w.r.t. the
% pump pulse.
fiber_calc_pulse_speed = fiber;
fiber_calc_pulse_speed.SR = fiber.SR*1e-30;
output_calc_pulse_speed = UPPE_propagate(fiber_calc_pulse_speed,initial_condition,sim,gas);

%%
save(sprintf('redshift_%s_%.1fm_%uum_%2.1fbar.mat',gas.gas_material,fiber.L0,gas.core_radius*2*1e6,gas.pressure/1e5),'-v7.3',...
     'c','dt','f','f_range','fiber','gas','lambda','Nt','num_save','prop_output','sim','t','time_window','wavelength_range','output_calc_pulse_speed');