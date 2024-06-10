% This code demonstrates the SPM in a N2-filled capillary.
% In Raman-active gases, the strongly-SPM-broadened pulse features a 
% negative chirp at the leading edge and positive chirp from the center to 
% the trailing edge.
%
% This aims to duplicate Fig. 2(a,b) of the following paper:
% Truong et al., "Spectral broadening and pulse compression in molecular 
% gas-filled hollow-core fibers," in IEEE Journal of Selected Topics in 
% Quantum Electronics (2024).

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.5,3]*1e-6; % m
Nt = 2^13;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SPM';
sim.pulse_centering = false;

num_save = 10;
fiber.L0 = 3.5; % propagation length
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
gas.core_radius = 250e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 10*1e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1025e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.gas_material = 'N2';
gas.fiber_type = 'MWLW_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.gas_material).(Raman_type).(Raman_parameters)
% 
% e.g.
%    gas.N2.R.T2 - N2's coherence decay time between the ground state and the excited rotational level
%    gas.N2.R.gamma - N2's polarizability anisotropy
%    gas.N2.R.polarization_calibration - calibration factor to the paper's polarizability anisotropy value (gamma)
%    gas.N2.R.B0 - N2's B0 value for rotational Raman scattering
%    gas.N2.R.D0 - N2's D0 value for rotational Raman scattering
%    gas.N2.R.alpha_e - N2's alpha_e value for rotational Raman scattering
%    gas.N2.R.omega - N2's rotational transition angular frequencies
%    gas.N2.R.preR - N2's prefactors, representing Raman strength, of rotational Raman scattering
%
%    gas.N2.V.T2 - N2's coherence decay time between the ground state and the excited vibrational level
%    gas.N2.V.Dalpha - N2's polarizability derivative
%    gas.N2.V.Dgamma - N2's polarizability-anisotropy derivative
%    gas.N2.V.polarization_calibration - calibration factor to the paper's polarizability-derivative value (Dalpha,Dgamma)
%    gas.N2.V.Omega - N2's vibrational transition frequency for the ground rotational level Q(0), v = 0 --> 1, J = 0 --> 0 (unchanged J)
%    gas.N2.V.omega - N2's vibrational transition angular frequencies (calibrated with rotational transitions from Omega)
%    gas.N2.V.preR - N2's prefactors, representing Raman strength, of vibrational Raman scattering
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition
tfwhm = 0.28; % ps

total_energy = 1e6; % nJ
pump_wavelength = 1025e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});

%% Propagation
prop_output = UPPE_propagate(fiber,initial_condition,sim,gas);

%% Plot
figure;
h = plot(t,squeeze(abs(prop_output.fields).^2),'linewidth',2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('Temporal profile');
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',14);
legend('Input','Output');
xlim([-0.5,0.8]);

spectrum = squeeze(abs(fftshift(ifft(prop_output.fields),1)).^2);
spectrum = spectrum./max(spectrum); % normalized to one for better visualization
figure;
h = plot(f,spectrum,'linewidth',2);
xlabel('Frequency (THz)');
ylabel('PSD (norm.)');
title('Spectrum');
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',14);
legend('Input','Output');
xlim([150,500]);

%% Spectrogram
% Spectrogram of the output pulse shows that it has negative chirp at the
% leading edge and positive chirp from the center to the trailing edge,
% which is a feature of strong SPM in Raman-active gases.
calc_spectrogram(t,fprop_output.fields(:,:,end));

% Animations of spectrograms
log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [~,~,~,figs,ax] = calc_spectrogram(t,f,prop_output.fields(:,1,i),[-0.5,0.8],[600,2000],400,400,true,true,log_yes);
    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('SPM chirp');
exportVideo.FrameRate = 10;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);