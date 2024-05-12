% This code simulates the SSFS process with a 10-fs pulse in a H2-filled 
% HCF.

close all;  clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^17;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.scalar = true;
sim.ellipticity = 0; % linearly polarized for 0 and circularly polarized for 1
%sim.gpu_yes = false;

num_save = 100;
fiber.L0 = 100; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % um
gas.temperature = 288; % K
gas.pressure = 10*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%% Initial condition and Propagate
pump_wavelength = 1030e-9; % m

tfwhm = 0.01;
beta2 = diff(fiber.betas,2)/(2*pi*(f(2)-f(1)))^2;
beta2_pump = beta2(find(f>3e5/(pump_wavelength*1e9),1)); % ps^2/m
fiiii = fiber.SR;
lambda0 = sim.lambda0;
freq_shift = 3e5/(pump_wavelength*1e9) - sim.f0;
n2_pump = 0.65e-23*gas.pressure/1.013e5;
initial_condition = build_MMsoliton(tfwhm, beta2_pump, fiiii, lambda0, time_window, 1, Nt, {'ifft',freq_shift},1,0,n2_pump);

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Save
save('SSFS_H2.mat');