% This code compares the Raman gain of my derived equation to 
% "Measurement of Raman Gain Coefficients of Hydrogen, Deuterium, and
% Methane" by John J. Ottusch and David A. Rockwel

clearvars; close all;

Aeff = pi*(15e-6)^2; % Just assume some effective area; not really important

eta = 115; % gas density in amagat
gas.temperature = 298; % K
gas.pressure = eta*1.01325e5; % Pa
c = 299792458; % m/s
h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (SI unit)
au_polarizability = 1.64878e-41; % F*m^2; atomic unit
permittivity0 = 8.85e-12;
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

wavelength = (250:10000)'*1e-9; % nm
wavenumber = 1./wavelength*1e-2; % cm^-1

%% The calculated Raman gain of my derived equation

% T2
vib.T2 = 1e6/(pi*(8220+384*eta)); % ps

% polaribility
vib.Dalpha = 5.72e-17;

% calibration factor
% To match with the experiments of some papers
vib.polarizability_calibration = 1;
% Apply the calibration
vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;

% vibrational Raman shift
f_vib = 2916; % cm^(-1)
vib.Omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz

% Raman response prefactor
vib.preR = gas.Ng*vib.Dalpha^2/4./(vib.Omega*1e12);

% Raman gain calculation
preX3_V = vib.preR;
hr_max = abs(imag(vib.Omega/(2i*vib.Omega/vib.T2+1/vib.T2^2)*1e-12));
Q = 1/Aeff/(permittivity0*c/2)^2;
omegaR = (wavenumber-2916)*100*c*2*pi;
my_vib_Raman_gain = 2*preX3_V*Q*hr_max*(omegaR/4);

%% Comparison
figure;
h = plot(wavelength*1e9,my_vib_Raman_gain);
wavelength_paper = 532; % nm
paper_vib_Raman_gain = 1.26e-11/Aeff; % 1/m/W
hold on; plot(wavelength_paper,paper_vib_Raman_gain,'Marker','.','MarkerSize',30); hold off;
set(h,'linewidth',2); set(gca,'fontsize',25);
xlabel('Wavelength (nm)');
ylabel('\gamma_g^{sim} (1/m/W)');
legend('\gamma_g^{mine}','\gamma_g^{paper}');