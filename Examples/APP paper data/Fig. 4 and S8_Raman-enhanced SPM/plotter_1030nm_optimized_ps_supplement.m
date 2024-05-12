clearvars; close all;

addpath('../../../user_helpers/');

load('optimized_ps_Raman_enhanced_SPM_1030nm_supplement.mat');

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
xlim([-0.4,0.8]);
ylim([0,8]);
xlabel('Time (ps)');
ylabel('Power (GW)');
set(gca,'fontsize',25);
legend('TL','D');
print(gcf,'optimized_ps_dechirped_1030nm_supplement.pdf','-dpdf');

%% chirped Stokes pulse
figure;
plot(t,abs(Stokes_pulse.rejected_fields(:,:,end)).^2/1e9,'linewidth',2,'Color','k');
hold on;
plot(t,abs(Stokes_pulse.fields(:,:,end)).^2/1e9,'linewidth',2,'Color','r');
hold off;
xlim([-2,2]);
ylim([0,1.8]);
xlabel('Time (ps)');
ylabel('Power (GW)');
set(gca,'fontsize',25);
legend('P','S');
print(gcf,'optimized_ps_chirped_Stokes_103nm_supplement.pdf','-dpdf');