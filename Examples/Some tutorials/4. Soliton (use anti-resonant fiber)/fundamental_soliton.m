% This code generates the fundamental soliton, which maintains its shape
% during propagation, in a Ar-filled HCF.
% Loss of the fiber is ignored to see the maintaining feature of the
% soliton after a long propagation.

close all;  clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.8,3]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;

num_save = 30;
fiber.L0 = 1000; % m; propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
% Load default parameters like 
%
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_deltaZ.threshold = 1e-8; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.Raman_model = 1; Use the Raman model (exist only in Raman-active gases)
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.num_noise_photon_per_bin = 1; Include photon shot noise
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'gas_UPPE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

% Configure gas parameters for the gas_info().
% These parameters vary based on different fiber type to use.
% Please check each example of each fiber for details regarding what
% parameters are required.
gas.core_radius = 30e-6; % m
gas.temperature = 288; % K
gas.pressure = 10*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1030e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'Ar';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 2e-3; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.gas_material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

fiber.betas = real(fiber.betas); % ignore the loss

%% Initial condition
N = 1; % soliton number

pump_wavelength = 1030e-9; % m

tfwhm = 0.5;
beta2 = diff(fiber.betas,2)/(2*pi*(f(2)-f(1)))^2;
beta2_pump = beta2(find(f>3e5/(pump_wavelength*1e9),1)); % ps^2/m
fiiii = fiber.SR;
lambda0 = sim.lambda0;
freq_shift = 3e5/(pump_wavelength*1e9) - sim.f0;
n2_pump = 7.96e-24*gas.pressure/1.013e5; % 7.96e-24 is picked due to Ar
initial_condition = build_MMsoliton(tfwhm, beta2_pump, fiiii, lambda0, time_window, 1, Nt, {'ifft',freq_shift},N,0,n2_pump);

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([-2,2]);

% Spectrum
figure;
h = plot(f,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlabel('Frequency (THz)');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([286,296]);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);
xlim([-2,2]);

% Comparison of spectra
figure;
[x,y] = meshgrid(f,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('Frequency (THz)');
ylabel('z');
title('Spectrum during propagation');
set(gca,'fontsize',14);
xlim([286,296]);