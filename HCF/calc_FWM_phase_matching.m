% Calculate the phase matching relation between the pump, Stokes, and
% anti-Stokes wavelengths.

clearvars; close all;

num_disp_orders = 3;
use_gpu = false;%true;

load('info_MWLW_coating_Ar_300um_25atm.mat','beta','wavelength');

c = 2.99792458e-4; % speed of ligth; m/ps

Nf = size(beta,1);
num_modes = size(beta,2);

wavelength_min = 0.3; % um
wavelength_max = 2; % um

% Show the Stokes,pump,anti-Stokes phase matching condition 
wavelength_pump = 1.03; % um
FWM_Shift = 87.449;

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

f_pump = c*1e6/wavelength_pump; % THz
f_seed = f_pump - FWM_Shift;
f_idler = f_pump + FWM_Shift;

pump_idx = find(f>f_pump,1); f_pump = f(pump_idx); beta_pump = new_beta(pump_idx);
seed_idx = find(f>f_seed,1); f_seed = f(seed_idx); beta_seed = new_beta(seed_idx);
idler_idx = find(f>f_idler,1); f_idler = f(idler_idx); beta_idler = new_beta(idler_idx);

figure;
h = plot(f,new_beta,'linewidth',2);
hold on; h2 = plot([f_seed,f_pump,f_idler],[beta_seed,beta_pump,beta_idler],'linewidth',2); hold off;
xlim(c*1e6./[wavelength_max,wavelength_min]);
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
%legend('\beta_0','Vib Raman');
ylabel('\beta (1/mm)');