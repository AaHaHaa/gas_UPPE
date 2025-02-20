close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^9;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SSFS';
sim.pulse_centering = true;
sim.gpuDevice.Index = 2;

num_save = 100;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 10e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'H2';
gas.num_tubes = 7; % the number of tubes in the anti-resonant fiber
gas.r_tube = 12.5e-6; % m; the tube radius (not core!)
gas.t_tube = 300e-9; % m; the tube's wall thickness of anti-resonant fibers
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.1;
total_energy = 1e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,length(sim.midx),Nt,{'ifft',freq_shift},sqrt([1,1e-3*ones(1,length(sim.midx)-1)]));

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = spectrum./max(spectrum(:));

ccc = whitejet_lower(256);
figure;
pcolor(f,prop_output.z,permute(spectrum,[3,1,2]));
shading interp;
colormap(ccc);
set(gca,'fontsize',25);
xlabel('Frequency (THz)');
ylabel('Propagation (m)');
xlim([200,350]);
pcolor_plotbox(gca);
c = colorbar; ylabel(c,'PSD (norm.)'); set(c,'YTick',[0,1]); caxis([0,1]);
print(gcf,'SSFS_H2_transient.pdf','-dpdf');