% Gao et al., "From Raman frequency combs to supercontinuum generation in 
% nitrogen-filled hollow-core anti-resonant fiber," Laser & Photonics
% Reviews 16, 2100426 (2022)

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.2,10]*1e-6; % m
Nt = 2^19;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.gpuDevice.Index = 2;
sim.progress_bar_name = 'N2';

num_save = 1;
fiber.L0 = 10; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 13e-6; % m
gas.temperature = 300; % K
gas.pressure_in = 1*1.01325e5; % Pa
gas.pressure_out = 40*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 532e-9; % m
gas.gas_material = 'N2';
gas.fiber_type = 'AR_HC_PCF';
gas.delta = 210e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 3.7e-3; % loss factor
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info_gradient_pressure(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas) + 1i*log(10^(12.3/1e3/10/2));
fiber.betas(lambda>1170) = real(fiber.betas(lambda>1170)) + 1i*log(10^(180/1e3/10/2)); % high loss

%% Initial condition and Propagate
tfwhm = 20; % ps
total_energy = 3.8e3; % nJ
pump_wavelength = 532e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
save(sprintf('comb_%2.1fuJ.mat',total_energy/1e3));