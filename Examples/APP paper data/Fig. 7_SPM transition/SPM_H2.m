% This code finds the varying accumulated nonlinear phase w.r.t. pulse
% duration due to varying Raman-induced SPM in H2.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.1,200]*1e-6; % m
Nt = 2^15;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'SPM';
%sim.gpu_yes = false;
sim.gpuDevice.Index = 2;

num_save = 30;
fiber.L0 = 2; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/time_window; % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 1e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 800e-9; % m
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 1e-2; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);
fiber.betas = fiber.betas*0;

%% Initial condition and Propagate
peak_power = 2e3/0.3;

tfwhm_all = [0.015:0.005:0.05,0.06:0.01:0.3]; % ps

prop_output = cell(1,length(tfwhm_all));
for i = 1:length(tfwhm_all)
    tfwhm = tfwhm_all(i);
    total_energy = tfwhm*peak_power; % nJ
    pump_wavelength = 800e-9; % m
    freq_shift = c/pump_wavelength - sim.f0;
    initial_condition = build_MMgaussian(tfwhm,time_window,total_energy,length(sim.midx),Nt,{'ifft',freq_shift},sqrt([1,1e-3*ones(1,length(sim.midx)-1)]));

    prop_output{i} = UPPE_propagate(fiber,initial_condition,sim,gas);
end

%%
nonlinear_phase = zeros(1,length(tfwhm_all));
for i = 1:length(tfwhm_all)
    spectrum0 = abs(fftshift(ifft(prop_output{i}.fields(:,:,1)),1)).^2;
    spectrum = abs(fftshift(ifft(prop_output{i}.fields(:,:,end)),1)).^2;

    spectrum0(spectrum0<max(spectrum0)/100) = 0;
    spectrum(spectrum<max(spectrum)/100) = 0;
    
    domega0 = calc_RMS(f,spectrum0);
    domega = calc_RMS(f,spectrum);

    % Calculate the nonlinear phase based on SPM equation of a Gaussian pulse.
    % Check Eq.(4.1.17) in "Nonlinear Fiber Optics (5ed)" by Agrawal.
    nonlinear_phase(i) = sqrt( ((domega/domega0)^2-1)*3*sqrt(3)/4 );
end

save('SPM_H2.mat','-v7.3');