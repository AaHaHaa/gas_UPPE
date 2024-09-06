clearvars; close all;

addpath('../../user_helpers');

load('comb_3.8uJ.mat');

%% Maps
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,:)),1)).^2.*3e2./(3e2./f).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(lambda,prop_output.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); colormap(jet);
caxis([-80,0]); xlim([400,1200]); caxis([-60,0]);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (cm)');
title('Spectrum');
%print(gcf,'spectrumMap.jpg','-djpeg');

figure;
fixed_fields = fft(ifft(prop_output.fields).*exp(1i*2*pi*ifftshift(f,1).*permute(prop_output.t_delay,[2 3 1])));
pcolor(t,prop_output.z*100,permute(abs(prop_output.fields).^2,[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
c = colorbar; ylabel(c,'Power (W)');
xlim([-20,20]);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Propagation distance (cm)');
title('Pulse');
%print(gcf,'pulseMap.jpg','-djpeg');

%% Final spectrum
ii = size(prop_output.fields,3);
z = prop_output.z(ii)*100;

figure;
h = plot(t,abs(prop_output.fields(:,1,ii)).^2);
hold on; h2 = plot(t,abs(prop_output.fields(:,1,1)).^2); hold off;
legend(sprintf('At %3.2fcm',z),'At 0cm');
set(h,'linewidth',2); set(h2,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Power (W)');
%print(gcf,'field.jpg','-djpeg');

figure;
h = plot(3e5./f,log_spectrum(end,:).'-max(log_spectrum(end,:)));
xlim([400,1200]); ylim([-60,0]);
set(h,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (dB)');
%print(gcf,'spectrum.jpg','-djpeg');