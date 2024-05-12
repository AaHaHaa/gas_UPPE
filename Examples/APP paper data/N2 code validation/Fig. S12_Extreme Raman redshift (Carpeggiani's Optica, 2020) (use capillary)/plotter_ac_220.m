close all;

addpath('../../../../user_helpers');

if exist('fig','var')
    clear fig;
end

wave_range = [900,1300];
time_range = [-5,5];

spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength = squeeze(sum(spectrum_wavelength,2)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength); log_spectrum_wavelength = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = squeeze(sum(spectrum,2)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(lambda,prop_output.z,log_spectrum_wavelength); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim(wave_range);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (m)');
title('Spectral evolution');
%ylim([0,2]);
%print(gcf,'chirped_spectrumMap.jpg','-djpeg');

figure;
h = plot(lambda/1e3,spectrum_wavelength(end,:).'.*time_window^2/1e6);
xlim(wave_range/1e3);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('Spectrum (\muJ/\mum)');
%print(gcf,'spectrum.jpg','-djpeg');

E = zeros(length(t),length(prop_output.z));
for j = 1:length(prop_output.z)
    E(:,j) = circshift(prop_output.fields(:,1,j),round((prop_output.t_delay(j)-output_calc_pulse_speed.t_delay(j))/prop_output.dt));
end

figure;
pcolor(t,prop_output.z,permute(abs(E).^2,[2,1])); shading interp; colormap(jet); %caxis([-40,60]);
%pcolor(t,prop_output_i.z*100,permute(sum(abs(prop_output_i.fields).^2,2),[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
cmap = whitejet_lower; colormap(cmap); c = colorbar; ylabel(c,'Power (W)');
xlim(time_range);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
title('Pulse evolution');
%print(gcf,'chirped_pulseMap.jpg','-djpeg');

figure;
pcolor(t,output_calc_pulse_speed.z,permute(abs(output_calc_pulse_speed.fields(:,1,:)).^2,[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
%pcolor(t,prop_output_i.z*100,permute(sum(abs(prop_output_i.fields).^2,2),[3,1,2])); shading interp; colormap(jet); %caxis([-40,60]);
cmap = whitejet_lower; colormap(cmap); c = colorbar; ylabel(c,'Power (W)');
xlim(time_range);
set(gca,'fontsize',20);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
title('Pulse evolution without nonlinearity');
%print(gcf,'chirped_pulseMap.jpg','-djpeg');

epsilon_V = zeros(length(t),length(prop_output.z));
for j = 1:length(prop_output.z)
    epsilon_V(:,j) = circshift(prop_output.delta_permittivity(:,1,2,j),round((prop_output.t_delay(j)-output_calc_pulse_speed.t_delay(j))/prop_output.dt));
end

figure;
pcolor(t,prop_output.z,abs(epsilon_V).'); shading interp; 
cmap = whitejet_lower; colormap(cmap); c = colorbar;
xlim(time_range);
set(gca,'fontsize',20);
%xlim([-100,100]);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
title('\Delta\epsilon_V');
%print(gcf,'chirped_p12VMap.jpg','-djpeg');

epsilon_R = zeros(length(t),length(prop_output.z));
for j = 1:length(prop_output.z)
    epsilon_R(:,j) = circshift(prop_output.delta_permittivity(:,1,1,j),round((prop_output.t_delay(j)-output_calc_pulse_speed.t_delay(j))/prop_output.dt));
end

figure;
pcolor(t,prop_output.z,abs(epsilon_R).'); shading interp; 
%cmap = whitejet_centered(2^10); colormap(cmap); %caxis([-40,60]);
cmap = whitejet_lower; colormap(cmap); c = colorbar;
xlim(time_range);
set(gca,'fontsize',20);
%xlim([-100,100]);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
title('\Delta\epsilon_R');
%print(gcf,'chirped_p12R2Map.jpg','-djpeg');