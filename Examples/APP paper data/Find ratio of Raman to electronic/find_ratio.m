% This code finds the ratio of Raman to the electronic nonlinearities in
% H2.

clearvars; close all;

load('H2.mat');

c = 299792.458;
pump_wavelength = linspace(800,1300,1000)';
pump_omega = c./pump_wavelength*2*pi; % THz

%% contributions
X3 = interp1(f*2*pi,sim.X3,pump_omega); % H2 is from the old code which saves X3 in "sim"
permittivity0 = 8.85e-12;
electronic_contribution = 3/4*permittivity0*X3;

Raman_contribution_real = sum(gas.H2.R.preR*4./(gas.H2.R.omega*1e12))/4 + sum(gas.H2.V.preR./(gas.H2.V.omega*1e12))/4;
Raman_contribution_imag = sum(-gas.H2.R.preR*4/2*(gas.H2.R.T2*1e-12))*0 + sum(-gas.H2.V.preR(2)/2*(gas.H2.V.T2*1e-12));
Raman_contribution = Raman_contribution_real + 1i*Raman_contribution_imag;

total_contribution = electronic_contribution + Raman_contribution;

ratio = abs(imag(total_contribution)./real(total_contribution));

figure;
plot(pump_wavelength,ratio,'linewidth',2,'Color','b');
xlabel('Pump wavelength (nm)');
ylabel('Im[\kappa+R]/Re[\kappa+R]');
set(gca,'fontsize',25);
xlim([pump_wavelength(1),pump_wavelength(end)]);