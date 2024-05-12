% This code simulates the Stokes generation process with a chirped pulse, 
% which has a 200-fs transform-limited duration.
% Many simulations with varying chirped durations are employed to
% demonstrate the effect of chirped duration to the Stokes generation
% efficiency and the Raman-enhanced SPM effect.
% Please see Fig. 4 for details.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,3]*1e-6; % m
Nt = 2^12;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SPM';
sim.gpuDevice.Index = 2;
sim.pulse_centering = false;
sim.adaptive_deltaZ.max_deltaZ = 1e-4;

num_save = 1;
fiber.L0 = 0.5; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 150e-6; % m
gas.temperature = 273.15 + 25; % K
gas.pressure = 20*1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.fiber_type = 'MWLW_coating';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Initial condition and Propagate
log10_tfwhm_lim = log10(0.2);
tfwhm = 10.^linspace(log10_tfwhm_lim,log10_tfwhm_lim+log10(10),30);

total_energy = 500e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_pulse = build_MMgaussian(tfwhm(1),time_window,total_energy,1,Nt,{'ifft',freq_shift});

func = calc_chirp();

num_repeat = 10;
prop_output = cell(length(tfwhm),num_repeat);
for i = 1:length(tfwhm)
    this_tfwhm = tfwhm(i); % ps
    
    [~,chirped_pulse] = func.General( this_tfwhm,ifftshift(f*2*pi),ifft(input_pulse.fields),1 );
    initial_condition.dt = input_pulse.dt;
    initial_condition.fields = chirped_pulse;
    
    for j = 1:num_repeat
        prop_output{i,j} = UPPE_propagate(fiber,initial_condition,sim,gas);
    end
end

%%
save('Raman_enhanced_SPM_1030nm_paper.mat');