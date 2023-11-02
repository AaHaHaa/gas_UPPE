% Calculate the Taylor series coefficients of the propagation constant,
% which represents the phase velocity, group velocity, and group
% dispersion, etc.

clearvars; close all;

addpath('../../../user_helpers');

num_disp_orders = 3;
use_gpu = false;%true;

load('info_AR_HC_PCF_H2_60um_30atm_310nm.mat','beta','wavelength');

c = 2.99792458e-4; % speed of ligth; m/ps

Nf = size(beta,1);
num_modes = size(beta,2);

wavelength_min = 0.7; % um
wavelength_max = 2; % um

% Show the Stokes,pump,anti-Stokes phase matching condition 
wavelength_pump = 1.03; % um
Stokes_shift_V = 125;

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

f1_pump = c*1e6/wavelength_pump; % THz
f1_Stokes_V = f1_pump - Stokes_shift_V;
f1_AS_V = f1_pump + Stokes_shift_V;

pump1_idx = find(f>f1_pump,1); f1_pump = f(pump1_idx); beta_pump1 = new_beta(pump1_idx,1);
Stokes1_V_idx = find(f>f1_Stokes_V,1); f1_Stokes_V = f(Stokes1_V_idx); beta_Stokes1_V = new_beta(Stokes1_V_idx,1);
AS1_V_idx = find(f>f1_AS_V,1); f1_AS_V = f(AS1_V_idx); beta_AS1_V = new_beta(AS1_V_idx,1);

beta_pump2 = new_beta(pump1_idx,2);
beta_Stokes2_V = new_beta(Stokes1_V_idx,2);
beta_AS2_V = new_beta(AS1_V_idx,2);


coherence_wave_beta = beta_pump2 - beta_Stokes2_V;

figure;
h = plot(f,new_beta,'linewidth',2);
hold on; plot([f1_Stokes_V,f1_pump,f1_AS_V],[beta_Stokes1_V,beta_pump1,beta_AS1_V],'linewidth',2); hold off;
hold on; plot([f1_Stokes_V,f1_pump,f1_AS_V],[beta_Stokes2_V,beta_pump2,beta_pump2+coherence_wave_beta],'linewidth',2); hold off;
xlim(c*1e6./[wavelength_max,wavelength_min]);
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
%legend('\beta_0','Vib Raman');
ylabel('\beta (1/mm)');