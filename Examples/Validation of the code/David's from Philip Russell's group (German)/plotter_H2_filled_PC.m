clearvars; close all;

addpath('../../../user_helpers');

permittivity0 = 8.85e-12;

load('H2-filled_PCF.mat');

spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(3e5./f,prop_output.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-60,0]);
c = colorbar; ylabel(c,'Intensity (dB)');
xlim([400,2000]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (cm)');
title('Spectrum');
print(gcf,'H2-filled_PCF_spectrumMap.jpg','-djpeg');

fixed_fields = fft(ifft(prop_output.fields).*exp(1i*2*pi*ifftshift(f,1).*permute(prop_output.t_delay,[2 3 1])));
[~,max_idx] = max(squeeze(max(abs(prop_output.fields),[],1)));

figure;
pcolor(t,prop_output.z*100,permute(abs(fixed_fields).^2,[3,1,2])); shading interp; colormap(jet(2^10)); %caxis([-40,60]);
colormap(whitejet_lower(2^10))
c = colorbar; ylabel(c,'Power (W)');
set(gca,'fontsize',20);
xlim([-2,2]);
xlabel('Time (ps)');
ylabel('Propagation distance (cm)');
title('Pulse');
print(gcf,'H2-filled_PCF_pulseMap.jpg','-djpeg');

epsilon_R_V = zeros(length(t),length(prop_output.z));
for j = 1:length(prop_output.z)
    epsilon_R_V(:,j) = circshift(prop_output.delta_permittivity(:,1,2,j)/permittivity0,0*round(prop_output.t_delay(j)/prop_output.dt));
end

figure;
pcolor(t,prop_output.z*100,real(epsilon_R_V).'); shading interp; 
%cmap = whitejet_centered(2^10); colormap(cmap); %caxis([-40,60]);
colormap(bluewhitered(2^10))
c = colorbar;
set(gca,'fontsize',20);
xlim([-2,2]);
xlabel('Time (ps)');
ylabel('Propagation distance (cm)');
title('\Delta\epsilon_r^V');
print(gcf,'H2-filled_PCF_epsilonMap.jpg','-djpeg');

ii = size(prop_output.fields,3);%max_idx;

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
%print(gcf,'H2-filled_PCF_pulse.jpg','-djpeg');

figure;
h = plot(3e5./f,spectrum(1,:));
hold on; h2 = plot(3e5./f,spectrum(ii,:)); hold off;
legend('At 0cm', ['At ' num2str(prop_output.z(ii)*1e2) 'cm']);
set(h,'linewidth',2); set(h2,'linewidth',2);
xlim([1750,1850]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('PSD');
title('Spectrum');

figure;
h = plot(t*1e3,abs(prop_output.delta_permittivity(:,1,2,1))/permittivity0);
hold on; h2 = plot(t*1e3,abs(prop_output.delta_permittivity(:,1,2,max_idx)/permittivity0)); hold off;
l = legend('At 0cm', ['At ' num2str(prop_output.z(ii)*1e2) 'cm']); set(l,'location','northwest');
set(h,'linewidth',2); set(h2,'linewidth',2);
set(gca,'fontsize',20);
%xlim([-100,100]); ylim([-0.3,0.3]);
xlabel('Time (fs)');
ylabel('Intensity');
title('\Delta\epsilon_r^V');
%print(gcf,'H2-filled_PCF_epsilonV.jpg','-djpeg');