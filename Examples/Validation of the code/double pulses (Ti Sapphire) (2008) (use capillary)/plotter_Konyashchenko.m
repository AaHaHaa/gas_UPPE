% This code plots the figures of Konyschenko's paper
%
% Konyashchenko et al., "Frequency shifting of sub-100 fs laser pulses by 
% stimulated Raman scattering in a capillary filled with pressurized gas," 
% Appl. Phys. B 93, 455-461 (2008)

close all; clearvars;

addpath('../../../user_helpers');

data = readmatrix('Stokes spectrum in two-pulse approach.csv');

[wavenumber2,idx] = sort(data(:,1)); % cm^(-1)
Stokes_spectrum_in_double_pulse = data(idx,2)/max(data(:,2)); % norm.

data = readmatrix('Stokes spectrum in single-pulse approach.csv');

[wavenumber1,idx] = sort(data(:,1)); % cm^(-1)
Stokes_spectrum_in_single_pulse = data(idx,2)/max(data(:,2)); % norm.

figure;
plot(wavenumber1,Stokes_spectrum_in_single_pulse,'linewidth',2,'Color','b');
hold on;
plot(wavenumber2,Stokes_spectrum_in_double_pulse,'linewidth',2,'Color','r');
xlabel('Wavenumber (cm^{-1})');
ylabel('PSD (norm.)');
set(gca,'fontsize',20);
l = legend('single-pulse','2nd Stokes (double-pulse)','location','northwest');
set(l,'fontsize',12);
xlim([7600,9000]); ylim([0,1.3]);
print('Konyashchenko Stokes.pdf','-dpdf');