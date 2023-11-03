clearvars; close all;

Aeff = pi*(15e-6)^2; % Just assume some effective area; not really important

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
rot.polarizability_calibration = 1;%1/sqrt(2);
vib.polarizability_calibration = 1;%sqrt(2);
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
%{
t = linspace(0,max(rot.T2,vib.T2)*10,2^30)';
R = sum([sum(rot.preR.*exp(-t/rot.T2).*sin(rot.omega.*t),2),...
         sum(vib.preR.*exp(-t/vib.T2).*sin(vib.omega.*t),2)],2);
dt = mean(diff(t));
R = R/trapz(t,R);

Raman_time = trapz(t,t.*R);
%}

t = linspace(0,max(rot.T2,vib.T2)*1,2^20)';
R_rot = sum(rot.preR.*exp(-t/rot.T2).*sin(rot.omega.*t),2);
R_vib = sum(vib.preR.*exp(-t/vib.T2).*sin(vib.omega.*t),2);
figure;
plot(t,[R_vib,R_rot],'linewidth',2);
xlabel('Time (ps)');
ylabel('Raman response');
set(gca,'fontsize',20);
legend('Vib','Rot');
print('Raman response.jpg','-djpeg');