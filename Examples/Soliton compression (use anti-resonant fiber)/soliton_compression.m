% This code simulates soliton compression in a H2-filled HCF.

close all;  clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.3,3]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;

num_save = 1;
fiber.L0 = 0.3; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % um
gas.temperature = 288; % K
gas.pressure = 1*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.050; % ps
total_energy = 3e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Save
figure;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2,'Color','b');
xlim([-0.1,0.1]);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Power (W)');