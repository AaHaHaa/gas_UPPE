% This code compares the Raman gain of my derived equation to 
% William K. Bischel and Mark J. Dyer's "Wavelength dependence of the 
% absolute Raman gain coefficient for the Q(1) transition in H2".
%
% They have measured the Raman gain of H2. I used their experimental data 
% to check the validity of my equaiton.

clearvars; close all;

Aeff = pi*(150e-6)^2; % Just assume some effective area; not really important

eta = 20; % gas density in amagat 
gas.temperature = 273.15 + 25; % K
gas.pressure = eta*1.01325e5; % Pa
c = 299792458; % m/s
h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (SI unit)
au_polarizability = 1.64878e-41; % F*m^2; atomic unit
permittivity0 = 8.85e-12;
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

%% Raman gain from the paper
wavelength = (250:2400)'*1e-9; % m; vib Raman shift 4155 cm^(-1) corresponds to 2.4 um
wavenumber = 1./wavelength*1e-2; % cm^-1

dv = 309/eta*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta; % ps; Raman linewidth
paper_vib_Raman_gain = 9.37e6*52*eta/dv*(wavenumber-4155).*1./((8.48e4)^2-wavenumber.^2).^2*1e-2; % m//W
%disp(9.37e6*52*eta/dv*(vp-4155)*1/((8.48e4)^2-vp^2)^2);

%% The calculated Raman gain of my derived equation
% T1
T1 = 1.555e3; % ps
rot.T1 = T1;
vib.T1  = T1;

% T2
rot.T2 = 1e6/(pi*(6.15/eta+114*eta)); % ps
vib.T2 = 1e6/(pi*(309/eta*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta)); % ps

% polaribility
rot.gamma = 2.0239*au_polarizability;
vib.Dalpha = 3.54e-17;
vib.Dgamma = 2.57e-17;

% calibration factor
% To match with the experiments of some papers
rot.polarizability_calibration = 1.1;
vib.polarizability_calibration = 1.05*0.95;
% Apply the calibration
rot.gamma = rot.gamma*rot.polarizability_calibration;
vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;
vib.Dgamma = vib.Dgamma*vib.polarizability_calibration;

% nuclear spin statistical constant
gJ = @(J) mod(J,2)*2 + 1; % 1 if even and 3 if odd
max_J = 7; % only 7 rotational energy levels below the 1st vibrational energy level

% Energy of each rotational state
rot.B0 = 60.8; rot.D0 = 1.6e-2; rot.alpha_e = 3.06; % cm^(-1)

% vibrational Raman shift
f_vib = 4155; % cm^(-1)
vib.Omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz

rot_J = 0:max_J-2;
vib_J = 0:max_J;
EJ = @(J) (rot.B0*J.*(J+1) - rot.D0*J.^2.*(J+1).^2)*1e2*h*c; % in the unit "J"
Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population

% frequency shift of each rotational level
rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
vib.omega = vib.Omega - rot.alpha_e*vib_J.*(vib_J+1)*1e2*c*1e-12*2*pi; % 2*pi*THz

% Raman response prefactor
rot.preR = gas.Ng*1/15/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
vib.preR = gas.Ng*(vib.Dalpha^2+4/45*vib_J.*(vib_J+1)./(2*vib_J-1)./(2*vib_J+3)*vib.Dgamma^2)/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);

% vib Raman gain calculation
preX3_V = sum(vib.preR);
hr_V_max = abs(imag(vib.Omega/(2i*vib.Omega/vib.T2+1/vib.T2^2)*1e-12));
N_1 = 1/(permittivity0*c/2)^2; % normalization constant
omegaR_V = (wavenumber-4155)*100*c*2*pi;
my_vib_Raman_gain = 2*preX3_V*N_1*hr_V_max*(omegaR_V/4);

% rot Raman gain calculation
preX3_R = rot.preR(2);
hr_R_max = abs(imag(rot.omega(2)/(2i*rot.omega(2)/rot.T2+1/rot.T2^2)*1e-12));
omegaR_R = (c./wavelength-rot.omega(2)/2/pi*1e12)*2*pi;
my_rot_Raman_gain = 2*preX3_R*N_1*hr_R_max*(omegaR_R/4);

%{
% alpha from Phil Russell's group
alpha_Belli = 1.55e-41; % F/m^2
preX3_V_Belli = gas.Ng*alpha_Belli^2/2/hbar;
Belli_vib_Raman_gain = 2*preX3_V_Belli*N_1*hr_max*(omegaR/4);
%}
%% Comparison
figure;
yyaxis right;
plot(wavelength*1e9,my_vib_Raman_gain./paper_vib_Raman_gain,'linewidth',2,'Color','r');
xlabel('Wavelength (nm)');
ylabel('\gamma_g^{here}/\gamma_g^{B. and D.''s}');
set(gca,'fontsize',20,'YColor','r');
yyaxis left;
plot(wavelength*1e9,[my_vib_Raman_gain,paper_vib_Raman_gain]*1e11,'linewidth',2,'Color','b');
ylabel('\gamma_g (cm/GW)');
legend('\gamma_g^{here}','\gamma_g^{B. and D.''s}');
set(gca,'fontsize',20,'YColor','b');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3)*1.3,pos(4)]);
print('vib gain comparison.pdf','-dpdf');

figure;
plot(wavelength*1e9,my_rot_Raman_gain*1e11,'linewidth',2);
xlabel('Wavelength (nm)');
ylabel('\gamma_g^{rot} (cm/GW)');
set(gca,'fontsize',20);