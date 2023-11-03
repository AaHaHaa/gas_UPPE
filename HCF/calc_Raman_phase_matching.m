% Calculate the phase matching relation between the pump, Raman-Stokes, and
% Raman-anti-Stokes wavelengths.
%
% Two pulses are assumed. The first pulse excites phonons which drives the
% Raman process in the second pulse. This code plots out the wave vectors
% for analysis.

%clearvars; close all;

num_disp_orders = 3;
use_gpu = false;%true;

load('info_AR_HC_PCF_H2_30um_1atm_300nm.mat');
gas_material = 'H2';

c = 2.99792458e-4; % speed of ligth; m/ps

Nf = size(beta,1);
num_modes = size(beta,2);

wavelength_min = 0.6; % um
wavelength_max = 20; % um

% Show the Stokes,pump,anti-Stokes phase matching condition 
wavelength_pump1 = 299792.458/(299792.458/1030)/1e3; % um
wavelength_pump2 = 299792.458/(299792.458/2000)/1e3; % um
Stokes_shift_V = 124.38;

tfwhm1 = 10; % ps
total_energy1 = 3e6; % nJ
tfwhm2 = 10; % ps
total_energy2 = 8e6; % nJ
ratio_square_peak_power_to_Gaussian_peak_power = 0.939437278699651;
peak_power1 = total_energy1*1e-9/tfwhm1/1e-12*ratio_square_peak_power_to_Gaussian_peak_power; % W
peak_power2 = total_energy2*1e-9/tfwhm2/1e-12*ratio_square_peak_power_to_Gaussian_peak_power; % W

%% Calculate the propagation constants
wavelength = wavelength*1e6; % um
f_calc = c./wavelength*1e6; % THz; loaded from the file

f = linspace(f_calc(end),f_calc(1),Nf)'; % THz; resample for applying Taylor series expansion
if use_gpu
    beta = mySpline( flipud(f_calc),flipud(beta),f(5:end-4),6,'../cuda/' );
    f = f(5:end-4);
else
    abs_beta = interp1(f_calc,abs(beta),f,'pchip');
    ang_beta = interp1(f_calc,unwrap(angle(beta),[],1),f,'pchip');
    beta = abs_beta.*exp(1i*ang_beta);
end

omega = 2*pi*f; % angular frequencies in 1/ps
df = f(2)-f(1);
domega = 2*pi*df;
beta = real(beta); % beta in 1/m

%% Display the results
coo = distinguishable_colors(num_modes);

% Make beta more visually understandable by subtracting out the zeroth and
% first order terms
min_idx = find(f>c*1e6/wavelength_min,1);
max_idx = find(f<c*1e6/wavelength_max,1,'last');
fit_order = 7;
[betas_Taylor_coeff,~,mu] = polyfit(omega(max_idx:min_idx),real(beta(max_idx:min_idx,1)),fit_order);
beta_ref = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
new_beta = real(beta(max_idx:min_idx,1))-(beta_ref(1)+beta_ref(2)*(omega(max_idx:min_idx)-mu(1))/mu(2));
beta_ref = [beta_ref(1)-beta_ref(2)*mu(1)/mu(2) + (max(new_beta) + min(new_beta))/2;...
            beta_ref(2)/mu(2)];
new_beta = beta - (beta_ref(1)+beta_ref(2)*omega);

f1_pump = c*1e6/wavelength_pump1; % THz
f1_Stokes_V = f1_pump - Stokes_shift_V;
f1_AS_V = f1_pump + Stokes_shift_V;

pump1_idx = find(f>f1_pump,1); f1_pump = f(pump1_idx); beta_pump1 = new_beta(pump1_idx);
Stokes1_V_idx = find(f>f1_Stokes_V,1); f1_Stokes_V = f(Stokes1_V_idx); beta_Stokes1_V = new_beta(Stokes1_V_idx);
AS1_V_idx = find(f>f1_AS_V,1); f1_AS_V = f(AS1_V_idx); beta_AS1_V = new_beta(AS1_V_idx);

f2_pump = c*1e6/wavelength_pump2; % THz
f2_Stokes_V = f2_pump - Stokes_shift_V;
f2_AS_V = f2_pump + Stokes_shift_V;

pump2_idx = find(f>f2_pump,1); f2_pump = f(pump2_idx); beta_pump2 = new_beta(pump2_idx);
Stokes2_V_idx = find(f>f2_Stokes_V,1); f2_Stokes_V = f(Stokes2_V_idx); beta_Stokes2_V = new_beta(Stokes2_V_idx);
AS2_V_idx = find(f>f2_AS_V,1); f2_AS_V = f(AS2_V_idx); beta_AS2_V = new_beta(AS2_V_idx);


%% Nonlinear coefficient
switch gas_material
    case 'H2' % m^2/(W*atm)
              % This value is taken from Wahlstrand, et al., 
              % "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2" (2015)
              % Its value, after converted into X3, is close to the paper by Belli et al.,
              % "Vacuum-ultraviolet to infrared supercontinuum in hydrogen-filled photonic crystal fiber" Optica (2015)
              % with X3 = 2.206e-26 m^2/V^2 at standard conditions
        n2 = 0.65e-23;
    case 'N2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'O2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'air' % Calculate n2 for N2 and O2
               % Add them up with 79% N2 and 21% O2
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2_N2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_N2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2_O2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_O2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2 = 0.79*n2_N2 + 0.21*n2_O2; %clear n2_N2 n2_O2
        n2(isinf(n2)) = 0; % avoid the singularity
    case 'Xe' % m^2/(W*atm)
              % From Shu-Zee Alencious Lo, et al.,
              % "Pulse propagation in hollow-core fiber at high-pressure regime: application to compression of tens of ?J pulses and determination of nonlinear refractive index of xenon at 1.03um" Applied Optics (2018)
        n2 = 50.1e-24;
    case 'Ar' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 7.96e-24;
    case 'Ne' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.85e-24;
    case 'He' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.34e-24;
    case 'Kr' % m^2/(W*atm)
              % From Carsten Br�e, Ayhan Demircan, and G�nter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 18.9e-24;
    case 'CH4'
        n2 = 3.118e-23; % m^2/(W*atm)
end
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
current_folder = current_path(1:sep_pos(end));
dr = mean(diff(r)); dtheta = mean(diff(theta));
SR = calc_SR_tensors(mode_profiles,permute(r,[2,3,1]),dr,dtheta,struct('gpu_yes',use_gpu),current_folder);
gamma2 = n2*(2*pi*f2_pump*1e12)/(c*1e12)*SR;
gamma2_Stokes = n2*(2*pi*f2_Stokes_V*1e12)/(c*1e12)*SR;
gamma1 = n2*(2*pi*f1_pump*1e12)/(c*1e12)*SR;
gamma1_Stokes = n2*(2*pi*f1_Stokes_V*1e12)/(c*1e12)*SR;

%% Plot

coherence_wave1_beta = beta_pump1 + (gamma1-gamma1_Stokes*2)*peak_power1 - beta_Stokes1_V;
coherence_wave2_beta = beta_pump2 + (gamma2-gamma2_Stokes*2)*peak_power2 - beta_Stokes2_V;

figure;
h = plot(f,new_beta,'linewidth',2);
hold on; h2 = plot([f1_Stokes_V,f1_pump],[beta_Stokes1_V,beta_pump1],'linewidth',2); hold off;
hold on; h3 = plot([f2_pump,f2_Stokes_V],[beta_pump2,beta_pump2-coherence_wave1_beta],'linewidth',2); hold off;
hold on; h4 = plot([f2_pump,f2_AS_V],[beta_pump2,beta_pump2+coherence_wave1_beta],'linewidth',2); hold off;
xlim(c*1e6./[wavelength_max,wavelength_min]);
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
%legend('\beta_0','Vib Raman');
ylabel('\beta (1/m)');

fprintf('phase-matching length (pump1): %6.4f(mm*atm)\n',2*pi*1e3/abs(beta_AS1_V+beta_Stokes1_V-beta_pump1*2)*(pressure/1.01325e5));
fprintf('phase-matching length (pump2): %6.4f(mm*atm)\n',2*pi*1e3/abs(beta_AS2_V+beta_Stokes2_V-beta_pump2*2)*(pressure/1.01325e5));
%fprintf('Phonon-amplifying length: %6.4f(mm*atm)\n',2*pi*1e3/abs(beta_Stokes2_V-(beta_pump2-coherence_wave1_beta))*(pressure/1.01325e5));
fprintf('Phonon-amplifying length: %6.4f(mm*atm)\n',2*pi*1e3/abs(coherence_wave1_beta-coherence_wave2_beta)*(pressure/1.01325e5));
fprintf('Anti-Stokes-amplifying length: %6.4f(mm*atm)\n',2*pi*1e3/abs(beta_AS2_V-(beta_pump2+coherence_wave1_beta))*(pressure/1.01325e5));