clearvars; close all;

addpath('../../../user_helpers/');

load('Raman_enhanced_SPM.mat');

initial_photon = trapz(f,abs(fftshift(ifft(prop_output.fields(:,:,1)),1)).^2./f);
Stokes_pulse = gaussian_spectral_filter(prop_output, sim.f0, c/(c/pump_wavelength-gas.H2.V.omega(2)/2/pi)*1e9, 1200, 5,true);

Stokes_QE = trapz(f,abs(fftshift(ifft(Stokes_pulse.fields(:,:,end)),1)).^2./f)/initial_photon;

%% Dechirped pulse
[~,dechirped_FWHM,dechirped_field] = pulse_compressor('Treacy-t',22*pi/180,sim.lambda0*1e9,t,Stokes_pulse.fields(:,:,end),1e-3/600,false,true);
[transform_limited_field,t_insert,transform_limited_FWHM] = calc_transform_limited( Stokes_pulse.fields(:,:,end),3,t );
figure;
plot(t_insert,abs(transform_limited_field).^2/1e9,'linewidth',2,'Color','b');
hold on;
plot(t,abs(dechirped_field).^2/1e9,'linewidth',2,'Color','r');
hold off;
xlim([-0.8,2]);
ylim([0,2]);
xlabel('Time (ps)');
ylabel('Power (GW)');
set(gca,'fontsize',25);
legend('TL','D');
%print(gcf,'dechirped.pdf','-dpdf');

%% Spectrogram
% Spectrogram of the pump shows that its center disappears and is
% transferred to the Stokes pulse such that Stokes pulse exhibits pump's
% SPM central segment, which further makes the Stokes pulse dechirpable to
% a duration shorter than pump's initial temporal duration.

% Filter
pump_pulse = gaussian_spectral_filter(prop_output, sim.f0, pump_wavelength*1e9, 100, 5);

calc_spectrogram(t,f,pump_pulse.fields(:,:,end)); title('Pump''s spectrogram');
calc_spectrogram(t,f,Stokes_pulse.fields(:,:,end)); title('Stokes''s spectrogram');

%% chirped Stokes pulse
figure;
plot(t,abs(Stokes_pulse.rejected_fields(:,:,end)).^2/1e9,'linewidth',2,'Color','k');
hold on;
plot(t,abs(Stokes_pulse.fields(:,:,end)).^2/1e9,'linewidth',2,'Color','r');
hold off;
xlim([-4,4]);
ylim([0,0.5]);
xlabel('Time (ps)');
ylabel('Power (GW)');
set(gca,'fontsize',25);
legend('P','S');
%print(gcf,'chirped_Stokes.pdf','-dpdf');