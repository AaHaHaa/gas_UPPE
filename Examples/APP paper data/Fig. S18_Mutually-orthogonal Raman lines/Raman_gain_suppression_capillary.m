close all; clearvars;

addpath('../../../broadband UPPE algorithm','../../../user_helpers');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,30]*1e-6; % m
Nt = 2^16;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
sim.gpu_yes = false;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
omega = f*2*pi;
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate([],sim);

gas.core_radius = 150e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 1.25*1.01325e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Include Raman gain suppression
fit_order = 7;
[betas_Taylor_coeff,~,mu] = polyfit(omega,real(fiber.betas),fit_order);
beta_ref = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
new_beta = real(fiber.betas)-(beta_ref(1)+beta_ref(2)*(omega-mu(1))/mu(2));
beta_ref = [beta_ref(1)-beta_ref(2)*mu(1)/mu(2) + (max(new_beta) + min(new_beta))/2;...
            beta_ref(2)/mu(2)];
new_beta = real(fiber.betas - (beta_ref(1)+beta_ref(2)*omega));

Stokes_shift_V = gas.H2.V.omega(2)/2/pi;
Stokes_shift_R = gas.H2.R.omega(2)/2/pi;

new_range = f>min(f)+Stokes_shift_V & f<max(f)-Stokes_shift_V;
f_pump = f(new_range); % THz
f_S_V = f_pump - Stokes_shift_V;
f_AS_V = f_pump + Stokes_shift_V;

f_S_R = f_pump - Stokes_shift_R;
f_AS_R = f_pump + Stokes_shift_R;

beta_pump = new_beta(new_range);
S_V_idx = arrayfun(@(x)find(f>x,1),f_S_V); f_S_V = f(S_V_idx); beta_S_V = new_beta(S_V_idx);
AS_V_idx = arrayfun(@(x)find(f>x,1),f_AS_V); f_AS_V = f(AS_V_idx); beta_AS_V = new_beta(AS_V_idx);

S_R_idx = arrayfun(@(x)find(f>x,1),f_S_R); f_S_R = f(S_R_idx); beta_S_R = new_beta(S_R_idx);
AS_R_idx = arrayfun(@(x)find(f>x,1),f_AS_R); f_AS_R = f(AS_R_idx); beta_AS_R = new_beta(AS_R_idx);

omega_P = f_pump*2*pi*1e12;

X3 = fiber.X3(find(lambda<gas.mode_profile_wavelength*1e9,1));

tfwhm = 1; % ps
total_energy = 1.5e6; % nJ
ratio_square_peak_power_to_Gaussian_peak_power = 0.939437278699651;
peak_power = total_energy*1e-9/tfwhm/1e-12*ratio_square_peak_power_to_Gaussian_peak_power; % W
peak_intensity = peak_power*fiber.SR;

transient_ratio = tfwhm/gas.H2.R.T2;

permittivity0 = 8.8541878176e-12;
I_Rb = sum(6*gas.H2.R.preR.*gas.H2.R.omega./(gas.H2.R.omega.^2+1./gas.H2.R.T2.^2)*1e-12)*peak_intensity;

delta_k_V = 2*beta_pump - beta_AS_V - beta_S_V;
delta_k_R = 2*beta_pump - beta_AS_R - beta_S_R;

R_f_vib_co = sum(4*gas.H2.R.preR.*(gas.H2.R.omega./(gas.H2.R.omega.^2-gas.H2.V.omega(2)^2+2i*gas.H2.V.omega(2)/gas.H2.R.T2+1/gas.H2.R.T2^2))*1e-12) + sum(gas.H2.V.preR.*(gas.H2.V.omega./(gas.H2.V.omega.^2-gas.H2.V.omega(2)^2+2i*gas.H2.V.omega(2)/gas.H2.V.T2+1/gas.H2.V.T2^2))*1e-12);
R_f_rot_co = sum(4*gas.H2.R.preR.*(gas.H2.R.omega./(gas.H2.R.omega.^2-gas.H2.R.omega(2)^2+2i*gas.H2.R.omega(2)/gas.H2.R.T2+1/gas.H2.R.T2^2))*1e-12) + sum(gas.H2.V.preR.*(gas.H2.V.omega./(gas.H2.V.omega.^2-gas.H2.R.omega(2)^2+2i*gas.H2.R.omega(2)/gas.H2.V.T2+1/gas.H2.V.T2^2))*1e-12);

R_f_cross = sum(6*gas.H2.R.preR.*(gas.H2.R.omega./(gas.H2.R.omega.^2-gas.H2.R.omega(2)^2+2i*gas.H2.R.omega(2)/gas.H2.R.T2+1/gas.H2.R.T2^2))*1e-12);

% Scaled from steady-state
R_f_vib_co = real(R_f_vib_co) +1i*imag(R_f_vib_co)*transient_ratio;
R_f_rot_co = real(R_f_rot_co) +1i*imag(R_f_rot_co)*transient_ratio;
R_f_cross = real(R_f_cross) +1i*imag(R_f_cross)*transient_ratio;

% Transient equation with pulse energy
R_f_vib_co = real(R_f_vib_co) -1i*gas.H2.V.preR(2)*total_energy*1e-9/peak_power/2;
R_f_rot_co = real(R_f_rot_co) -1i*4*gas.H2.R.preR(2)*total_energy*1e-9/peak_power/2;
R_f_cross = real(R_f_cross) -1i*6*gas.H2.R.preR(2)*total_energy*1e-9/peak_power/2;

kappa = 1/permittivity0^2/(c*1e12)^2;
kappa_e = 3*permittivity0*X3/4;

g_vib_co = kappa*gas.H2.V.omega(2)*1e12*imag(R_f_vib_co)*peak_intensity + ...
           1/2*abs(real(sqrt((2*kappa*((sqrt(omega_P.^2-(gas.H2.V.omega(2)*1e12)^2)+omega_P).*(kappa_e+R_f_vib_co))*peak_intensity-delta_k_V).*...
                             (2*kappa*((sqrt(omega_P.^2-(gas.H2.V.omega(2)*1e12)^2)-omega_P).*(kappa_e+R_f_vib_co))*peak_intensity+delta_k_V))));
g_rot_co = kappa*gas.H2.R.omega(2)*1e12*imag(R_f_rot_co)*peak_intensity + ...
           1/2*abs(real(sqrt((2*kappa*((sqrt(omega_P.^2-(gas.H2.R.omega(2)*1e12)^2)+omega_P).*(kappa_e+R_f_rot_co))*peak_intensity-delta_k_R).*...
                             (2*kappa*((sqrt(omega_P.^2-(gas.H2.R.omega(2)*1e12)^2)-omega_P).*(kappa_e+R_f_rot_co))*peak_intensity+delta_k_R))));
                         
G_vib_co = 2*real(g_vib_co)/peak_intensity;
G_rot_co = 2*real(g_rot_co)/peak_intensity;

ap = 2*kappa*(sqrt(omega_P.^2-(gas.H2.R.omega(2)*1e12)^2)*(R_f_cross/2+kappa_e/3)*peak_intensity+omega_P*((R_f_cross/2-kappa_e/3)*peak_intensity-I_Rb));
am = 2*kappa*(sqrt(omega_P.^2-(gas.H2.R.omega(2)*1e12)^2)*(R_f_cross/2+kappa_e/3)*peak_intensity-omega_P*((R_f_cross/2-kappa_e/3)*peak_intensity-I_Rb));
g_cross_linear = 1/2*kappa*gas.H2.R.omega(2)*1e12*imag(R_f_cross)*peak_intensity+...
                 1/2*abs(real(sqrt((ap-delta_k_R).*(am+delta_k_R))));

G_cross_linear = 2*real(g_cross_linear)/peak_intensity;

figure;
h = plot(lambda(new_range)/1e3,[G_rot_co,G_vib_co,G_cross_linear]*1e11,'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
xlim([0.7,2.2]);
set(gca,'fontsize',25);
xlabel('Wavelength (μm)');
ylabel('Raman gain (cm/GW)');
legend('co;Rot','co;Vib','cross;Rot');
print(gcf,'Raman gain_for_interleaving_lines.pdf','-dpdf');

% Inset around gain suppression
figure;
h = plot(lambda(new_range)/1e3,[G_rot_co,G_vib_co]*1e11,'linewidth',10);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlim([0.95,1.08]); set(gca,'XTick',[0.95,1,1.05]);
ylim([0,2.08e-7]);
set(gca,'fontsize',60,'YTick',[]);
print(gcf,'Raman_gain_closeview.pdf','-dpdf');