close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.4,5]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.pulse_centering = false;

num_save = 20;
fiber.L0 = 0.6; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Initial condition and Propagate
tfwhm = 0.3; % ps
total_energy = 30e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

tfwhm2 = 0.6; % ps
func = calc_chirp;
omega = ifftshift(2*pi*f,1); % 2*pi*THz
[~,chirped_pulse] = func.Gaussian( tfwhm2,omega,ifft(input_field.fields),1 );
input_field.fields = chirped_pulse;

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 33e-6; % m
gas.temperature = 288.15; % K
gas.pressure = 23*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.delta = 310e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%%
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
save(sprintf('SM_%ubar.mat',gas.pressure/1.01325e5));