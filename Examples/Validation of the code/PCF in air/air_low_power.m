% This code duplicates the simulations of air in Fig. 4 in "Nonlinear
% dynamic of picosecond pulse propagation in atmospheric air-filled hollow
% core fibers (2018)"
%
% Mousavi et al., "Nonlinear dynamic of picosecond pulse propagation in 
% atmospheric air-filled hollow core fibers," Opt. Express 26, 8866-8882
% (2018)

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.8,3]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;

num_save = 50;
fiber.L0 = 10; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 21.5e-6; % m
gas.temperature = 288.15; % K
gas.pressure = 1*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'air';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%% Initial condition and Propagate
tfwhm = 6; % ps
total_energy = 30e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,-time_window/12);

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
save('air.mat');