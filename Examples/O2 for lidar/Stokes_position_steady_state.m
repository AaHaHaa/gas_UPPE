% This code simulates the initial Stokes growth in the steady-state Raman
% regime.
% In the steady-state regime, the Stokes growth is the strongest at the
% peak power which is at the pulse center of a Gaussian pump pulse.

close all; clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^20;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);
sim.gpuDevice.Index = 1;
sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'Duplicate O2';
sim.pulse_centering = false;

num_save = 1;
fiber.L0 = 20; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 20e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1064e-9; % m
gas.material = {'O2'};
gas.num_tubes = 7; % the number of tubes in the anti-resonant fiber
gas.r_tube = 40e-6; % m; the tube radius (not core!)
gas.t_tube = 300e-9; % m; the tube's wall thickness of anti-resonant fibers
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%gas.O2.V.preR = gas.O2.V.preR/2;

%% Initial condition and Propagate
tfwhm = 0.5e3; % ps
total_energy = 20e3*0.7; % nJ
pump_wavelength = 1064e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,4);

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
Stokes_pulse = edgepass_spectral_filter('highpass', prop_output, sim.f0, 1080);
Stokes_pulse.fields = Stokes_pulse.fields./sqrt(max(abs(Stokes_pulse.fields).^2,[],1));
pump_pulse = edgepass_spectral_filter('lowpass', prop_output, sim.f0, 1080);
pump_pulse.fields = pump_pulse.fields./sqrt(max(abs(pump_pulse.fields).^2,[],1));

figure;
h = plot(t/1e3,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2,'linewidth',2);
set(gca,'fontsize',25);
set(h(1),'Color','k'); set(h(2),'Color','r');
xlim([-1,1]*tfwhm/1e3*2);
xlabel('Time (ns)');
ylabel('Power (norm.)');
legend('P','S');

factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain
spectrum = abs(fftshift(ifft(prop_output.fields(:,:,end)))).^2*factor_correct_unit.*factor;

figure;
semilogy(lambda,spectrum/max(spectrum),'linewidth',2);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
xlim([1000,1750]);
%ylim([1e-3,1]);
set(gca,'fontsize',25);

output_field = prop_output;
output_field.fields = output_field.fields(:,:,end);

center_lambda = 1064;
bandwidth_lambda = 100;
gaussexpo = 4;
pump = gaussian_spectral_filter(output_field, sim.f0, center_lambda, bandwidth_lambda, gaussexpo);
disp((trapz(t,abs(pump.fields).^2)/1e3)/(trapz(t,abs(input_field.fields).^2)/1e3));

center_lambda = 1280;
bandwidth_lambda = 100;
gaussexpo = 4;
S1 = gaussian_spectral_filter(output_field, sim.f0, center_lambda, bandwidth_lambda, gaussexpo);
disp((trapz(t,abs(S1.fields).^2)/1e3)/(trapz(t,abs(input_field.fields).^2)/1e3));