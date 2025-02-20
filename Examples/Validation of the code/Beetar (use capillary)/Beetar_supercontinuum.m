% This code duplicates the simulations of N2 under a 400-uJ and 280-fs
% pulse in Fig. 3 in "Multioctave supercontinuum generation and frequency 
% conversion based on rotational nonlinearity (2020)"
%
% Beetar et al., "Multioctave supercontinuum generation and frequency conversion based on rotational nonlinearity," Sci. Adv. 6, eabb5375 (2020)

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.4,2]*1e-6; % m
Nt = 2^11;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
sim.progress_bar_name = 'N2(Beetar)';

num_save = 50;
fiber.L0 = 3.5; % m; propagation length
sim.save_period = fiber.L0/num_save;
%sim.midx = 1:3;
sim.scalar = false;

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
gas.core_radius = 250e-6; % m
gas.temperature = 300; % K
gas.pressure = 6.2*1.01325e5; % Pa; gas pressure
gas.wavelength_order = 6; % The code recomputes the propagation constant to ensure that it has smooth higher-order derivatives up this order; set this to 6, currently maximum implemented value in mySpline.cu, all the time
gas.mode_profile_wavelength = 1025e-9; % m; the wavelength of the mode profile used to compute SR values and overlap integrals, etc.
gas.material = 'N2';
gas.fiber_type = 'no_coating'; % 'Ag_coating', 'no_coating', 'MWLW_coating' coating types for capillaries
gas.xy_sampling = 101; % spatial sampling number for computing the mode profiles for SR values and overlap integrals, etc.

% Load hollow-core fiber parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
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
total_energy = 400e3; % nJ
pump_wavelength = 1025e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift});
input_field.fields = input_field.fields.*sqrt([1,0.01*ones(1,length(sim.midx)*2-1)]);

%% Propagation
prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%% Plot
spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2*time_window^2/1e6.*3e5./(3e5./f).^2; % uJ/nm

final_spectrum = spectrum(:,:,end);
final_spectrum = final_spectrum/max(final_spectrum(:));
figure;
semilogy(3e5./f,sum(final_spectrum,2),'linewidth',2);
xlim([700,1800]);
ylim([1e-3,1]);
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');

%% Save the results
save('Beetar_supercontinuum.mat');