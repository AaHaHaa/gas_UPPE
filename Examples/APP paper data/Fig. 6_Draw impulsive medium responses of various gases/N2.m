% This code draws the impulsive response of N2.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.15,10]*1e-6; % m
Nt = 2^17;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'medium response';
sim.pulse_centering = false;
sim.gpu_yes = false;

num_save = 1;
fiber.L0 = 0.00001; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'N2';
gas.num_tubes = 7; % the number of tubes in the anti-resonant fiber
gas.r_tube = 12.5e-6; % m; the tube radius (not core!)
gas.t_tube = 300e-9; % m; the tube's wall thickness of anti-resonant fibers
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagatefigure;
tfwhm1 = 0.01; % ps
total_energy = 100; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm1,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,1);

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
de = imag(prop_output.delta_permittivity(:,1,1,1));
de = de/max(de);
P = abs(prop_output.fields(:,:,1)).^2;
P = P/max(P);
figure;
plot(t,de,'linewidth',2,'Color','k');
ylabel('\Delta\epsilon (norm.)');
ylim([-1.15,1.15]);
xlim([-0.5,15]);
xlabel('Time (ps)');
set(gca,'fontsize',25);
print('N2_response.pdf','-dpdf');

%%
save('N2.mat');