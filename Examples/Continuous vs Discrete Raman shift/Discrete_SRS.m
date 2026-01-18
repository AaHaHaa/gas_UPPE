close all; clearvars;

addpath('../../user_helpers','../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^15;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'Discrete SRS';
sim.pulse_centering = false;

num_save = 100;
fiber.L0 = 0.3; % m; propagation length
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
gas.core_radius = 50e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 100*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 800e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = {'H2'};
gas.fiber_type = 'no_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
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

%% Single-pulse pumping
tfwhm = 5; % ps
total_energy = 70e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim(1).f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1);
%tfwhm = 5; % ps
%func = calc_chirp;
%[~,chirped_field] = func.Gaussian( tfwhm,ifftshift(2*pi*f,1),ifft(input_field.fields),1 );
%input_field.fields = chirped_field;

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Animation
center_lambda = 1800;
bandwidth_lambda = 500;
gaussexpo = 4;
Stokes = gaussian_spectral_filter(prop_output, sim.f0, center_lambda, bandwidth_lambda, gaussexpo);
pump = gaussian_spectral_filter(prop_output, sim.f0, pump_wavelength*1e9, bandwidth_lambda, gaussexpo);

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
    h = plot(t,abs([pump.fields(:,:,i),Stokes.fields(:,:,i)]).^2/max_PowerT,'linewidth',2);
    set(h(1),'Color','k','linestyle','-'); set(h(2),'Color','r','linestyle','-');
    xlabel('Time');
    ylabel('Power (norm.)');
    set(gca,'XTick',0,'XTickLabel',{'0'});
    set(gca,'YTick',[0,1]);
    xlim([-10,10]);
    ylim([0,1.5]);
    set(gca,'fontsize',25);
    legend(h,'pump','Stokes');

    subplot(1,2,2);
    h = plot(lambda,abs(fftshift(ifft([pump.fields(:,:,i),Stokes.fields(:,:,i)]),1)).^2*factor_correct_unit.*factor/max_PSD,'linewidth',2);
    set(h(1),'Color','k','linestyle','-'); set(h(2),'Color','r','linestyle','-');
    xlabel('Wavelength');
    ylabel('PSD (norm.)');
    set(gca,'XTick',[pump_wavelength*1e9,center_lambda],'XTickLabel',{'\lambda_P','\lambda_S'});
    set(gca,'YTick',[0,1]);
    xlim([950,1900]);
    ylim([0,1]);
    set(gca,'fontsize',25);
    set(figs,'Position',[400,290,840,420]);
    legend(h,'pump','Stokes');
    drawnow;

    set(figs,'Color',[1,1,1]);
    
    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('Discrete_SRS');
exportVideo.FrameRate = 10;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

%% Save the results
save('Discrete_SRS.mat');