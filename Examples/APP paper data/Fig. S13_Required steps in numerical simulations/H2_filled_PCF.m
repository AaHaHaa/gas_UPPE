% I stopped in the MATLAB's debug mode to plot separately the spectrum for
% showing some figures in Fig. S13.
% How this code is set up isn't important.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.6,2]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'any name';
sim.pulse_centering = false;
sim.gpu_yes = false;

num_save = 50;
fiber.L0 = 1; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 10*1.013e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'H2';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas); % remove the loss
gas.H2.R.preR = gas.H2.R.preR*0;

%% Initial condition and Propagate
tfwhm1 = 0.05; % ps
total_energy = 30e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm1,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,1);

func = calc_chirp();
[~,chirped_pulse] = func.General( 1,ifftshift(f*2*pi),ifft(input_field.fields),1 );
initial_condition.dt = input_field.dt;
initial_condition.fields = chirped_pulse;

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);