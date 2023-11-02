% This code duplicates the simulations of H2 under a 1540-nm pulse in 
% Fig. 14(c) in "Pure rotational stimulated Raman scattering in H2-filled 
% hollow-core photonic crystal fibers (2020)"

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [1.38,1.75]*1e-6; % m
Nt = 2^22;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);
sim.f0 = f0;
sim.progress_bar_name = 'H2 S(1)';

num_save = 1;
fiber.L0 = 10; % propagation length
sim.save_period = 0;%fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 4.5e-6/0.64; % m; their fiber has 9-um MFD, so I need to back-calculate the actual core radius.
gas.temperature = 300; % K
gas.pressure = 16*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1540e-9; % m
gas.gas_material = 'H2';
gas.fiber_type = 'no_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas(lambda<1600) = real(fiber.betas(lambda<1600)) + 1i*log(10^(0.04/10/2));
fiber.betas(lambda>1600) = real(fiber.betas(lambda>1600)) + 1i*log(10^(0.11/10/2));
gas.H2.R.preR = [gas.H2.R.preR([1,2]),zeros(1,length(gas.H2.R.preR)-2)];
gas.H2.V.preR = zeros(size(gas.H2.V.preR)); % no vib

%% Initial condition and Propagate
tfwhm = 12e3; % ps
total_energy = 2.5e3; % nJ
pump_wavelength = 1540e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Plot
filter_lambda = 1600;
probe = struct('dt',prop_output.dt,'fields',prop_output.fields(:,:,end));
residual = edgepass_filter('lowpass', probe, sim.f0, filter_lambda,0.3);
Raman = edgepass_filter('highpass', probe, sim.f0, filter_lambda,0.3);
residual_energy = trapz(t,abs(residual.fields).^2)/1e3; % nJ
Raman_energy = trapz(t,abs(Raman.fields).^2)/1e3; % nJ

spectrum = abs(fftshift(ifft(probe.fields(:,1,:)),1)).^2*time_window^2/1e6.*3e5./(3e5./f).^2; % uJ/nm
spectrum = squeeze(spectrum(:,1,:)).';

final_spectrum = spectrum(end,:)';
final_spectrum = final_spectrum/max(final_spectrum);
figure;
plot(3e5./f,final_spectrum,'linewidth',2);
xlim([1500,1750]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');

%%
save(sprintf('S1_Raman_gain_%u.mat',total_energy),'-v7.3','t','filter_lambda','lambda','probe','pump_wavelength','Raman_energy','residual_energy','tfwhm','time_window','total_energy');