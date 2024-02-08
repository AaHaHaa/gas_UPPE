close all; clearvars;

addpath('../../../MMTools/gas_UPPE_vector/user_helpers','../../../MMTools/gas_UPPE_vector/broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,3]*1e-6; % m
Nt = 2^15;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SSFS';
sim.pulse_centering = true;
sim.scalar = false;

num_save = 1;
fiber.L0 = 0.1; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 150e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 20*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm = 1; % ps
total_energy = 1.5e6; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,length(sim.midx),Nt,{'ifft',freq_shift},sqrt([1,1e-3*ones(1,length(sim.midx)-1)]));
initial_condition.fields = fft(ifft(initial_condition.fields).*sqrt([1,0.0001]).*exp(1i*2*pi*[zeros(Nt,1),rand(Nt,1)]));

prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%%
spectrum = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2;
spectrum = spectrum./max(spectrum(:));
idx = f>265 & f<280;
max_plot = max(spectrum(idx,2))*1.05;

figure;
h = plot(f,spectrum/max_plot,'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('Frequency (THz)');
ylabel('PSD (norm.)');
ylim([0,1]);
set(gca,'fontsize',25);
xlim(c/pump_wavelength+[-1,1]*150);
legend('x','y');
print(gcf,'capillary_PSD_1d5mJ_1ps.pdf','-dpdf');