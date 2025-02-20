close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,30]*1e-6; % m
Nt = 2^16;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'two-pulse circular';
sim.gpuDevice.Index = 2;
sim.pulse_centering = false;
sim.scalar = false;
sim.ellipticity = 1;

num_save = 20;
fiber.L0 = 0.2; % propagation length
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
gas.material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
tfwhm1_1 = 0.03; % ps
total_energy1 = 200e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
time_delay1 = -20;
input_field1 = build_MMgaussian(tfwhm1_1,time_window,total_energy1,1,Nt,{'ifft',freq_shift},1,time_delay1);
input_field1.fields = input_field1.fields.*sqrt([1,0.01]);

tfwhm1_2 = 10; % ps
total_energy2 = 10e3; % nJ
time_delay2 = -time_delay1;
input_field20 = build_MMgaussian(tfwhm1_2,time_window,total_energy2,1,Nt,{'ifft',freq_shift},1,time_delay2);

prop_output = cell(1,2);
for i = 1:2
    if i == 1
        input_field2.fields = input_field20.fields.*[1,1]/sqrt(2);
    else
        input_field2.fields = input_field20.fields.*[1,-1]/sqrt(2);
    end
    
    initial_condition = struct('dt',dt,'fields',input_field1.fields + input_field2.fields);
    prop_output{i} = UPPE_propagate(fiber,initial_condition,sim,gas);
end
%%
%figure;
%plot(t,abs(prop_output.fields(:,:,end)).^2);
z = 21;

output{1}.fields = (prop_output{1}.fields(:,1,:) + prop_output{1}.fields(:,2,:))/sqrt(2);
output{2}.fields = (prop_output{2}.fields(:,1,:) - prop_output{2}.fields(:,2,:))/sqrt(2);

idx = (f>272 & f<276) | (f>308 & f<310);
spectrum = abs(fftshift(ifft([output{1}.fields(:,:,z),output{2}.fields(:,:,z)]),1)).^2;
norm_factor = max(max(spectrum(idx,:)));

idx = t<0;
pulse2_co = output{1}.fields(:,1,z); pulse2_co(idx) = 0;
pulse2_cross = output{2}.fields(:,1,z); pulse2_cross(idx) = 0;

spectrum_co = abs(fftshift(ifft(pulse2_co),1)).^2/norm_factor;
spectrum_cross = abs(fftshift(ifft(pulse2_cross),1)).^2/norm_factor;

figure;
h = plot(f,[spectrum_co,spectrum_cross]);
set(h(1),'Color','b','linewidth',4); set(h(2),'Color','r','linewidth',2);
set(gca,'fontsize',25);
xlabel('Frequency (THz)');
ylabel('PSD (norm.)');
legend('x','y');
ylim([0,1.5]);
xlim(c/pump_wavelength+[-1,1]*50);
print(gcf,'two_pulse_circular.pdf','-dpdf');