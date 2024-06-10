close all; clearvars;

load('SPM.mat');

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('t (ps)');
ylabel('Power (W)');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([-0.2,0.2]);

% Spectrum
figure;
c = 299792.458; % nm/ps
wavelength = c./f; % nm
factor_correct_unit = time_window^2/1e6; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e6" is to make pJ into uJ
factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
spectrum_wavelength = abs(fftshift(ifft(prop_output.fields),1)).^2*factor_correct_unit.*factor;
h = plot(lambda,spectrum_wavelength(:,:,end));
xlabel('Wavelength (nm)');
ylabel('PSD (\muJ/nm)');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);