close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^15;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'Temporal modulation';
sim.gpuDevice.Index = 2;
sim.pulse_centering = false;

num_save = 1;
fiber.L0 = 1; % propagation length
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
gas.gas_material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm0 = 0.3;
tfwhm = 3;

% Pump
total_energy = 1000e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
pump_pulse = build_MMgaussian(tfwhm0,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,1);

% Chirp them
func = calc_chirp();
[~,chirped_pump] = func.General( tfwhm,ifftshift(f*2*pi),ifft(pump_pulse.fields),1 );

initial_condition.dt = pump_pulse.dt;
initial_condition.fields = chirped_pump;

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
save('Stokes_temporal_modulation_chirp_without_seed.mat');