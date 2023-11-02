close all;

addpath('../../../user_helpers');

load('Raman_52fs.mat');

range = {Nt/2:Nt,1:Nt/2};
probe_pulse.fields = prop_output{2}.fields(:,1,end);
probe_pulse.dt = dt;
probe_pulse.fields(range{2},:,:) = 0;

initial.fields = prop_output{2}.fields(:,1,1);
initial.dt = dt;
initial.fields(range{2},:,:) = 0;

%%
lambda1 = 600; % nm
probe_pulse1 = gaussian_spectral_filter(probe_pulse, sim.f0, lambda1, 100);

disp(sum(abs(probe_pulse1.fields).^2)/sum(abs(initial.fields).^2)*100);

%%
lambda2 = 800; % nm
probe_pulse2 = gaussian_spectral_filter(probe_pulse, sim.f0, lambda2, 300);

disp(sum(abs(probe_pulse2.fields).^2)/sum(abs(initial.fields).^2)*100);

%%
lambda3 = 1200; % nm
probe_pulse3 = gaussian_spectral_filter(probe_pulse, sim.f0, lambda3, 300);

disp(sum(abs(probe_pulse3.fields).^2)*dt/1e3);

%% plot
[optimal_value,dechirped_FWHM,dechirped_field,grating_size,mirror_size] = pulse_compressor('Treacy-beta2',pi/6,sim.lambda0*1e9,t,probe_pulse3.fields,1e-3/1000,false,false,-1);
transform_limited_field = calc_transform_limited( probe_pulse3.fields );
figure;
h = plot(t*1e3,abs([transform_limited_field,dechirped_field]).^2/1e9,'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('Time (fs)');
ylabel('Power (GW)');
xlim([-300,300]);
set(gca,'fontsize',20); %set(gcf,'Color','None');
legend('TL','D');
print('compressed_pulse_S.eps','-depsc');

spectrum_wavelength = abs(fftshift(ifft(probe_pulse3.fields),1)).^2.*299.792458./(299.792458./f).^2*time_window^2/1e6;
figure;
plot(lambda/1e3,spectrum_wavelength,'linewidth',2*1.5,'Color','b');
set(gca,'fontsize',20*1.5);
xlabel('Wavelength (μm)');
ylabel('PSD (μJ/μm)');
xlim([1.1,1.3]);
print('Stokes_spectrum.eps','-depsc');

figure;
h = plot(t,abs([probe_pulse1.fields,probe_pulse2.fields,probe_pulse3.fields]).^2/1e6,'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','k'); set(h(3),'Color','r');
set(gca,'fontsize',20);
legend('AS','P','S');
xlim([0,10]);
xlabel('Time (ps)'); ylabel('Power (MW)');
print('each_component.pdf','-dpdf');