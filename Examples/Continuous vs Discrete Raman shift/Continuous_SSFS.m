% This code simulates soliton self-frequency shift (SSFS) in a 30-um-core
% anti-resonant hollow-core fiber filled with 15-bar H2.
% The input pulse is 50 fs at 1030 nm.
% This is a scalar simulation and doesn't include polarization modes. The
% input pulse is linearly polarized.

close all;  clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.pulse_centering = false;
sim.ellipticity = 0; % linearly polarized for 0 and circularly polarized for 1

num_save = 100;
fiber.L0 = 50; % propagation length
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
% sim.adaptive_dz.threshold = 1e-8; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = true; Consider the Raman (exist only in Raman-active gases)
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
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
gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 5*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1030e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = {'H2'};
gas.num_tubes = 7; % the number of tubes in the anti-resonant fiber
gas.r_tube = 20e-6; % m; the tube radius (not core!)
gas.t_tube = 300e-9; % m; the tube's wall thickness of anti-resonant fibers
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
% 
% e.g.
%    gas.H2.R.T1 - H2's upper population lifetime for the excited rotational level (unused)
%    gas.H2.R.T2 - H2's coherence decay time between the ground state and the excited rotational level
%    gas.H2.R.gamma - H2's polarizability anisotropy
%    gas.H2.R.polarization_calibration - calibration factor to the paper's polarizability anisotropy value (gamma)
%    gas.H2.R.B0 - H2's B0 value for rotational Raman scattering
%    gas.H2.R.D0 - H2's D0 value for rotational Raman scattering
%    gas.H2.R.alpha_e - H2's alpha_e value for rotational Raman scattering
%    gas.H2.R.omega - H2's rotational transition angular frequencies
%    gas.H2.R.preR - H2's prefactors, representing Raman strength, of rotational Raman scattering
%
%    gas.H2.V.T1 - H2's upper population lifetime for the excited vibrational level (unused)
%    gas.H2.V.T2 - H2's coherence decay time between the ground state and the excited vibrational level
%    gas.H2.V.Dalpha - H2's polarizability derivative
%    gas.H2.V.Dgamma - H2's polarizability-anisotropy derivative
%    gas.H2.V.polarization_calibration - calibration factor to the paper's polarizability-derivative value (Dalpha,Dgamma)
%    gas.H2.V.Omega - H2's vibrational transition frequency for the ground rotational level Q(0), v = 0 --> 1, J = 0 --> 0 (unchanged J)
%    gas.H2.V.omega - H2's vibrational transition angular frequencies (calibrated with rotational transitions from Omega)
%    gas.H2.V.preR - H2's prefactors, representing Raman strength, of vibrational Raman scattering
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = real(fiber.betas);

%% Initial condition
N = 1; % soliton number

pump_wavelength = 1030e-9; % m

tfwhm = 0.05;
beta2 = diff(fiber.betas,2)/(2*pi*(f(2)-f(1)))^2;
beta2_pump = beta2(find(f>3e5/(pump_wavelength*1e9),1)); % ps^2/m
fiiii = fiber.SR;
lambda0 = sim.lambda0;
freq_shift = 3e5/(pump_wavelength*1e9) - sim.f0;
n2_pump = 0.65e-23*gas.pressure/1.013e5; % 0.65e-23 is picked due to H2 [see gas_n2.m]
initial_condition = build_MMsoliton(tfwhm, beta2_pump, fiiii, lambda0, time_window, 1, Nt, {'ifft',freq_shift},N,0,n2_pump);

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Animation
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

max_PowerT = max(abs(prop_output.fields(:)).^2);
max_PSD = max(abs(fftshift(ifft(prop_output.fields),1)).^2*factor_correct_unit.*factor,[],'all');

% Temporal and spectral evolutions
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    subplot(1,2,1);
    plot(t,abs(prop_output.fields(:,:,i)).^2/max_PowerT,'linewidth',2,'Color','k');
    xlabel('Time');
    ylabel('Power (norm.)');
    set(gca,'XTick',[0,13.28],'XTickLabel',{'0','\Deltat'});
    set(gca,'YTick',[0,1]);
    xlim([-1,14]);
    ylim([0,1]);
    set(gca,'fontsize',25);
    subplot(1,2,2);
    plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,i)),1)).^2*factor_correct_unit.*factor/max_PSD,'linewidth',2,'Color','k');
    xlabel('Wavelength');
    ylabel('PSD (norm.)');
    set(gca,'XTick',[1030,1188.7],'XTickLabel',{'\lambda_0','\lambda_0+\Delta\lambda'});
    set(gca,'YTick',[0,1]);
    xlim([950,1300]);
    ylim([0,1]);
    set(gca,'fontsize',25);
    set(figs,'Position',[400,290,840,420]);
    drawnow;

    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('Continuous_SSFS');
exportVideo.FrameRate = 10;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

%% Save the results
save('Continuous_SSFS.mat');