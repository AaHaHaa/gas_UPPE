% This code finds where there is no Kerr-induced Raman suppression effect
% resulting from differential nonlinear phase accumulation between the pump
% and the Stokes pulses.

clearvars; close all;

load('H2.mat');

c = 299792.458;
pump_wavelength = linspace(800,1300,1000)';
pump_omega = c./pump_wavelength*2*pi; % THz
Stokes_omega = pump_omega - gas.H2.V.omega(2);

%% contributions
X3 = interp1(f*2*pi,sim.X3,pump_omega);
permittivity0 = 8.85e-12;
electronic_contribution = 3/4*permittivity0*X3;

Raman_contribution = sum(gas.H2.R.preR./(gas.H2.R.omega*1e12))*4 + sum(gas.H2.V.preR./(gas.H2.V.omega*1e12));

phase_diff = (pump_omega - 2*Stokes_omega).*electronic_contribution + (pump_omega - 5/4*Stokes_omega)*Raman_contribution;

figure;
plot(pump_wavelength,phase_diff,'linewidth',2,'Color','b');
xlabel('Pump wavelength (nm)');
ylabel('Phase factor');
set(gca,'fontsize',25);
xlim([pump_wavelength(1),pump_wavelength(end)]);