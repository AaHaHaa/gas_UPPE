clearvars; close all;

addpath('../../../../MMTools/gas_UPPE/user_helpers');

%%
load('two_color_lossy3_0.70um2000uJ_10ps_2um10uJ_10ps_300umID_20bar_L100cm.mat');

permittivity0 = 8.85418782e-12;
h = 6.62607015e-34;

range = {Nt/2:Nt,1:Nt/2};

% compute spectra
probe = prop_output;
probe.fields(range{2},:,:) = 0;

spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';

% photon number
f_AS = f<299792.458/1.0e3 & f>299792.458/1.2e3;
f_pump = f<299792.458/1.7e3 & f>299792.458/2.3e3;
f_Stokes = f<299792.458/5e3 & f>299792.458/15e3;
photon_number = spectrum*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
total_photon_number = trapz(f*1e12,photon_number.',1)';
photon_number_AS = trapz(f(f_AS)*1e12,photon_number(:,f_AS).',1)';
photon_number_pump = trapz(f(f_pump)*1e12,photon_number(:,f_pump).',1)';
photon_number_S = trapz(f(f_Stokes)*1e12,photon_number(:,f_Stokes).',1)';

fig = figure;
h = plot(probe.z,[total_photon_number,photon_number_pump,photon_number_AS,photon_number_S]/total_photon_number(1),'linewidth',10);
set(h(1),'LineStyle','--','Color','k');
set(h(2),'Color','k'); set(h(3),'Color','b'); set(h(4),'Color','r');
set(gca,'Color','None','XTick',[],'YTick',[]);
ylim([0,1]);
print(gcf,'photon number (700nm).pdf','-dpdf');

%%
range = {Nt/2:Nt,1:Nt/2};
probe_pulse.fields = prop_output.fields(:,1,25);
probe_pulse.dt = dt;
probe_pulse.fields(range{2},:,:) = 0;

cutonoff_lambda = 7e3; % nm
probe_pulse = edgepass_filter('highpass', probe_pulse, sim.f0, cutonoff_lambda);

[optimal_value,dechirped_FWHM,dechirped_field,grating_size,mirror_size] = pulse_compressor('Treacy-beta2',pi/6,sim.lambda0*1e9,t,probe_pulse.fields,1e-3/50,false,false,-1);
transform_limited_field = calc_transform_limited( probe_pulse.fields );
max_P = max(max(abs([transform_limited_field,dechirped_field]).^2/1e6));
figure;
h = plot(t,abs([transform_limited_field,dechirped_field]).^2/1e6,'linewidth',4);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlim([-1,1]);
ylim([0,max_P*1.1]);
set(gca,'fontsize',50,'YTick',[]); %set(gcf,'Color','None');
legend('TL','D');
print('dechirped_pulse_S (700nm).eps','-depsc');

%%
load('two_color_lossy4_1.09um2000uJ_10ps_2um10uJ_10ps_300umID_20bar_L100cm.mat');

permittivity0 = 8.85418782e-12;
h = 6.62607015e-34;

range = {Nt/2:Nt,1:Nt/2};

% compute spectra
probe = prop_output;
probe.fields(range{2},:,:) = 0;

spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';

% photon number
f_AS = f<299792.458/1.0e3 & f>299792.458/1.2e3;
f_pump = f<299792.458/1.7e3 & f>299792.458/2.3e3;
f_Stokes = f<299792.458/5e3 & f>299792.458/15e3;
photon_number = spectrum*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
total_photon_number = trapz(f*1e12,photon_number.',1)';
photon_number_AS = trapz(f(f_AS)*1e12,photon_number(:,f_AS).',1)';
photon_number_pump = trapz(f(f_pump)*1e12,photon_number(:,f_pump).',1)';
photon_number_S = trapz(f(f_Stokes)*1e12,photon_number(:,f_Stokes).',1)';

fig = figure;
h = plot(probe.z,[total_photon_number,photon_number_pump,photon_number_AS,photon_number_S]/total_photon_number(1),'linewidth',10);
set(h(1),'LineStyle','--','Color','k');
set(h(2),'Color','k'); set(h(3),'Color','b'); set(h(4),'Color','r');
set(gca,'Color','None','XTick',[],'YTick',[]);
ylim([0,1]);
print(gcf,'photon number (1090nm).pdf','-dpdf');

%%
range = {Nt/2:Nt,1:Nt/2};

cutonoff_lambda = 1090; % nm

% good pulse
probe_pulse.fields = prop_output.fields(:,1,34);
probe_pulse.dt = dt;
probe_pulse.fields(range{2},:,:) = 0;
probe_pulse = gaussian_spectral_filter(probe_pulse, sim.f0, cutonoff_lambda, 100);

[optimal_value,dechirped_FWHM,dechirped_field,grating_size,mirror_size] = pulse_compressor('Treacy-beta2',pi/6,sim.lambda0*1e9,t,probe_pulse.fields,1e-3/50,false,false,-1);
transform_limited_field = calc_transform_limited( probe_pulse.fields );
max_P = max(max(abs([transform_limited_field,dechirped_field]).^2/1e6));
figure;
h = plot(t,abs([transform_limited_field,dechirped_field]).^2/1e6,'linewidth',4);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlim([-1,1]);
ylim([0,max_P*1.1]);
set(gca,'fontsize',50,'YTick',[]); %set(gcf,'Color','None');
legend('TL','D');
print('dechirped_pulse_AS_good (1090nm).eps','-depsc');

% bad pulse due to back-conversion
probe_pulse.fields = prop_output.fields(:,1,62);
probe_pulse.dt = dt;
probe_pulse.fields(range{2},:,:) = 0;
cutonoff_lambda = 1090; % nm
probe_pulse = gaussian_spectral_filter(probe_pulse, sim.f0, cutonoff_lambda, 100);

[optimal_value,dechirped_FWHM,dechirped_field,grating_size,mirror_size] = pulse_compressor('Treacy-beta2',pi/6,sim.lambda0*1e9,t,probe_pulse.fields,1e-3/50,false,false,-1);
max_P = max(abs(dechirped_field).^2/1e6);
figure;
h = plot(t,abs(dechirped_field).^2/1e6,'linewidth',4);
set(h(1),'Color','r');
xlim([-1,1]);
ylim([0,max_P*1.1]);
set(gca,'fontsize',50,'YTick',[]); %set(gcf,'Color','None');
print('dechirped_pulse_AS_bad (1090nm).eps','-depsc');