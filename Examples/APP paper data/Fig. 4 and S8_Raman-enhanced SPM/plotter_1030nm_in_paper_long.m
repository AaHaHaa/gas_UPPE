clearvars; close all;

addpath('../../../user_helpers/');

load('Raman_enhanced_SPM_1030nm_paper_long.mat');

Stokes_intensity = zeros(Nt,num_repeat);
Stokes_QE = zeros(num_save+1,num_repeat);
for ii = 1
    for j = 1:num_repeat
        initial_photon = trapz(f,abs(fftshift(ifft(prop_output{ii,j}.fields(:,:,1)),1)).^2./f);
        Stokes_pulse = gaussian_spectral_filter(prop_output{ii,j}, sim.f0, c/(c/pump_wavelength-gas.H2.V.omega(2)/2/pi)*1e9, 400, 5);

        Stokes_intensity(:,j) = abs(Stokes_pulse.fields(:,:,end)).^2;
        
        Stokes_QE(:,j) = squeeze(trapz(f,abs(fftshift(ifft(Stokes_pulse.fields),1)).^2./f)/initial_photon);
    end
end
Stokes_intensity_mean = mean(Stokes_intensity,2); Stokes_intensity_std = std(Stokes_intensity,[],2);
norm_factor = max(Stokes_intensity_mean+Stokes_intensity_mean);

Stokes_QE_mean = mean(Stokes_QE,2); Stokes_QE_std = std(Stokes_QE,[],2);

figure;
confplot(t,Stokes_intensity_mean/norm_factor,Stokes_intensity_std/norm_factor,'linewidth',2,'Color','r');
hold on;
plot(t,abs(prop_output{1}.fields(:,:,1)).^2/max(abs(prop_output{1}.fields(:,:,1)).^2),'linewidth',2,'Color','k');
ylabel('Peak power (norm.)');
ylim([0,1]);
set(gca,'fontsize',25);
xlim([-4,4]);
xlabel('Time (ps)');
print(gcf,'pulse_1030nm_long.pdf','-dpdf');

figure;
confplot(prop_output{1}.z,Stokes_QE_mean*100,Stokes_QE_std*100,'linewidth',2,'Color','r');
ylabel('QE (%)');
set(gca,'fontsize',25);
xlim([0,fiber.L0]);
xlabel('Propagation distance (m)');
print(gcf,'QEevolution_1030nm_long.pdf','-dpdf');
this_ylim = get(gca,'YLim');

figure;
confplot(prop_output{1}.z,Stokes_QE_mean*100,Stokes_QE_std*100,'linewidth',2,'Color','r');
set(gca,'fontsize',40);
xlim([0.4,fiber.L0]);
ylim(this_ylim);
set(gca,'XTick',[0.4,0.5])
print(gcf,'QEevolution_1030nm_long_closeview.pdf','-dpdf');