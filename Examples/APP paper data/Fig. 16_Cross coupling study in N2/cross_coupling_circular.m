% This code simulates the depolarization effect when the pump pulse is
% circularly polarized.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.8,3]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'cross coupling';
%sim.gpu_yes = false;
sim.scalar = false;
sim.ellipticity = 1;

num_save = 30;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 150e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 3*1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'N2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.3; % ps
total_energy = 100e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,length(sim.midx),Nt,{'ifft',freq_shift},sqrt([1,1e-3*ones(1,length(sim.midx)-1)]));
initial_condition.fields = initial_condition.fields.*sqrt([1,0.01]);

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = spectrum./max(spectrum(:));

figure;
h = plot(f,spectrum(:,:,end)/max(max(spectrum(:,:,end))),'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('Frequency (THz)');
ylabel('PSD (norm.)');
set(gca,'fontsize',25);
xlim([270,310]);
legend('+','-');
print(gcf,'cross_coupling_circular.pdf','-dpdf');

figure;
h = plot(f,spectrum(:,:,10),'linewidth',10);
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',60,'Color','None','YTick',[]);
xlim([285,295]);
ylim([0,0.12]);
print(gcf,'cross_coupling_circular_inset.pdf','-dpdf');