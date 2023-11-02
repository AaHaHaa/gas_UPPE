% This code duplicates the simulations of N2 under a 100 fs pulse in 
% Fig. 3(e) in "High energy redshifted and enhanced spectral broadening by 
% molecular alignment (2020)"

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.4,2]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
sim.progress_bar_name = 'N2(230fs)';

num_save = 50;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 250e-6; % m
gas.temperature = 288.15; % K
gas.pressure = 1*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 780e-9; % m
gas.gas_material = 'N2';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 0.23; % ps
total_energy = tfwhm*1e7; % nJ
pump_wavelength = 780e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Plot
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,:)),1)).^2*time_window^2/1e6.*3e5./(3e5./f).^2; % uJ/nm
spectrum = squeeze(spectrum(:,1,:)).';

figure;
h = plot(3e5./f,spectrum(end,:).');
xlim([600,950]);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');
print('230fs.jpg','-djpeg');

%%
save('N2_230fs.mat');