% This code simulates the double-pulse Raman generation with the first
% pulse at 1090 nm.
% At this first-pulse wavelength, phonon absorption is preferred in the
% second pulse.

close all; clearvars;

addpath('../../../../user_helpers','../../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,30]*1e-6; % m
Nt = 2^16;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = '1.09,2 LWIR';
sim.gpuDevice.Index = 2;
sim.pulse_centering = false;

num_save = 100;
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
gas.mode_profile_wavelength = 2000e-9; % m
gas.material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%% Initial condition and Propagate
tfwhm1_1 = 0.05; % ps
total_energy1 = 2000e3; % nJ
pump_wavelength1 = 1090e-9; % m
freq_shift = c/pump_wavelength1 - sim.f0;
time_delay1 = -20;
input_field1 = build_MMgaussian(tfwhm1_1,time_window,total_energy1,1,Nt,{'ifft',freq_shift},1,time_delay1);
tfwhm2_1 = 10; % ps
func = calc_chirp;
omega = ifftshift(2*pi*f,1); % 2*pi*THz
if total_energy1 ~= 0
    [~,chirped_pulse] = func.General( tfwhm2_1,omega,ifft(input_field1.fields),1 );
    input_field1.fields = chirped_pulse;
end

tfwhm1_2 = 0.05; % ps
total_energy2 = 5000e3; % nJ
pump_wavelength2 = 2000e-9; % m
freq_shift = c/pump_wavelength2 - sim.f0;
time_delay2 = -time_delay1;
input_field2 = build_MMgaussian(tfwhm1_2,time_window,total_energy2,1,Nt,{'ifft',freq_shift},1,time_delay2);
tfwhm2_2 = 10; % ps
func = calc_chirp;
omega = ifftshift(2*pi*f,1); % 2*pi*THz
[~,chirped_pulse] = func.General( tfwhm2_2,omega,ifft(input_field2.fields),1 );
input_field2.fields = chirped_pulse;

initial_condition = struct('dt',dt,'fields',input_field1.fields + input_field2.fields);
for i = 1:10
    prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

    %%
    save(sprintf('two_color_lossy%u_1.09um%uuJ_%ups_2um%uuJ_%ups_%uumID_%ubar_L%ucm.mat',i,total_energy1/1e3,tfwhm2_1,total_energy2/1e3,tfwhm2_2,gas.core_radius*1e6*2,gas.pressure/1e5,fiber.L0*100),'-v7.3');
end