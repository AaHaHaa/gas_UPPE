clearvars; close all;

addpath('../../../user_helpers');

load('MM_35bar.mat');

permittivity0 = 8.85418782e-12;

prop_output_i = prop_output;
prop_output_i.fields = prop_output_i.fields(:,1,:);

%% Maps
spectrum = abs(fftshift(ifft(prop_output_i.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(f,prop_output_i.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-60,0]);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
ylabel('Propagation distance (cm)');
title('Spectrum');
%print(gcf,'spectrumMap.jpg','-djpeg');

figure;
pcolor(3e5./f,prop_output_i.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-60,0]);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlim([600,2000]);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (cm)');
title('Spectrum (uncorrected lambda)');

figure;
fixed_fields = fft(ifft(prop_output_i.fields).*exp(1i*2*pi*ifftshift(f,1).*permute(prop_output_i.t_delay,[2 3 1])));
pcolor(t,prop_output_i.z*100,permute(abs(prop_output_i.fields).^2,[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
c = colorbar; ylabel(c,'Power (W)');
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Propagation distance (cm)');
title('Pulse');
%print(gcf,'pulseMap.jpg','-djpeg');

epsilon = zeros(length(t),length(prop_output_i.z));
for j = 1:length(prop_output_i.z)
    epsilon(:,j) = circshift(prop_output_i.delta_permittivity(:,1,2,j),round(prop_output_i.t_delay(j)/prop_output_i.dt));
end
epsilon_r = epsilon/permittivity0;

figure;
pcolor(t*1e3,prop_output_i.z*100,real(epsilon_r).'); shading interp; 
%cmap = whitejet_centered(2^10); colormap(cmap); %caxis([-40,60]);
colormap(bluewhitered(2^10))
c = colorbar;
set(gca,'fontsize',20);
xlabel('Time (fs)');
ylabel('Propagation distance (cm)');
title('\Delta\epsilon_r^V');
%print(gcf,'coherenceMap.jpg','-djpeg');

%% Final spectrum
ii = size(prop_output_i.fields,3);
z = prop_output_i.z(ii)*100;

figure;
h = plot(t,abs(prop_output_i.fields(:,1,1)).^2);
hold on; h2 = plot(t,abs(prop_output_i.fields(:,1,ii)).^2); hold off;
legend('At 0cm', sprintf('At %3.2fcm',z));
set(h,'linewidth',2); set(h2,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Power (W)');
%print(gcf,'field.jpg','-djpeg');

figure;
h = plot(3e5./f/1e3,log_spectrum(end,:).');
set(h,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('Spectrum (dB)');
%print(gcf,'spectrum.jpg','-djpeg');

%% Raman coherence and population
figure;
h = plot(t,real(epsilon_r(:,1)));
hold on; h2 = plot(t,epsilon_r(:,ii)); hold off;
l = legend('At 0cm', sprintf('At %3.2fcm',z)); set(l,'location','northwest');
set(h,'linewidth',2); set(h2,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('\Delta\epsilon_r^V');
%print(gcf,'p12V.jpg','-djpeg');