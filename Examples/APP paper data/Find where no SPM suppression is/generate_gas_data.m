% This code generates the mat file so that the other script can use it to 
% find where Kerr-induced suppression effect doesn't take into effect. 

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.6,5]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'David''s duplication';
sim.pulse_centering = false;
sim.num_photon_noise_per_band = 1;

num_save = 1;
fiber.L0 = 0.01; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 1*1.013e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'H2';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas); % remove the loss

%% Initial condition and Propagate
tfwhm1 = 1; % ps
total_energy = 30e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm1,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
save('H2.mat');