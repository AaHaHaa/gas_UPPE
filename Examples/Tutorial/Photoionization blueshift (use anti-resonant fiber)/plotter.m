%clearvars; close all;

addpath('../../../user_helpers');

%load('photoionization_blueshift_5.3uJ.mat');

permittivity0 = 8.85e-12;

spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(f,prop_output.z*100,log_spectrum); shading interp; colormap(jet); caxis([-40,0]);
cmap = whitejet_lower; colormap(cmap);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
ylabel('Propagation distance (cm)');
title('Spectrum');
%print(gcf,sprintf('H2-filled_PCF_spectrumMap_%u.jpg',i),'-djpeg');

%fixed_fields = fft(ifft(prop_output_i_i.fields).*exp(1i*2*pi*ifftshift(f,1).*permute(prop_output_i_i.t_delay,[2 3 1])));
[~,max_idx] = max(squeeze(max(abs(prop_output.fields),[],1)));

figure;
pcolor(t*1e3,prop_output.z*100,permute(abs(prop_output.fields).^2,[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
cmap = whitejet_lower; colormap(cmap);
c = colorbar; ylabel(c,'Power (W)');
set(gca,'fontsize',20);
%xlim([-100,100]);
xlabel('Time (fs)');
ylabel('Propagation distance (cm)');
title('Pulse');
%print(gcf,sprintf('H2-filled_PCF_pulseMap_%u.jpg',i),'-djpeg');

figure;
h = plot(lambda,spectrum(end,:));
set(h,'linewidth',2);
set(gca,'fontsize',20);
%xlim([120 220]); %ylim([-40,0]);
ylabel('Intensity (dB)');
xlabel('Wavelength (nm)');
title('Spectrum');
%print(gcf,sprintf('H2-filled_PCF_spectrum_%u.jpg',i),'-djpeg');

ii = max_idx;

figure;
h = plot(t*1e3,abs(prop_output.fields(:,1,1)).^2);
hold on; h2 = plot(t*1e3,abs(prop_output.fields(:,1,ii)).^2); hold off;
legend('At 0cm', ['At ' num2str(prop_output.z(ii)*1e2) 'cm']);
set(h,'linewidth',2); set(h2,'linewidth',2);
%xlim([-100,100]); %ylim([0,4e8]);
set(gca,'fontsize',20);
xlabel('Time (fs)');
ylabel('Power (W)');
title('Pulse');
%print(gcf,sprintf('H2-filled_PCF_pulse_%u.jpg',i),'-djpeg');