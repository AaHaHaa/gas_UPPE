%% Plotter for the polarized SSFS

%clearvars;
close all;

addpath('../../user_helpers');

freq_range = [166,375]; % [160,330];
wave_range = 3e5./fliplr(freq_range);
time_range = [-5,5];

%% Spectra in two modes
% the first polarization mode
spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength1 = squeeze(spectrum_wavelength(:,1,:)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength1); log_spectrum_wavelength1 = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum1 = log_spectrum - max(log_spectrum(:));
% the second polarization mode
spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength2 = squeeze(spectrum_wavelength(:,2,:)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength2); log_spectrum_wavelength2 = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = squeeze(spectrum(:,2,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum2 = log_spectrum - max(log_spectrum(:));

%% Spectra evolution of the first polarization mode
figure;
pcolor(lambda,prop_output.z,log_spectrum_wavelength1); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim(wave_range);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (m)');
title('Spectral evolution');

%% Spectra evolution of the first polarization mode
figure;
pcolor(lambda,prop_output.z,log_spectrum_wavelength2); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim(wave_range);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Propagation distance (m)');
title('Spectral evolution');

%% Output spectra
figure;
h = plot(lambda/1e3,[spectrum_wavelength1(end,:).',spectrum_wavelength2(end,:).'].*time_window^2/1e6);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('Spectrum (\muJ/\mum)');
xlim(wave_range/1e3);