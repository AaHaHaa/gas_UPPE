% This code simulates the initial Stokes growth in the steady-state Raman
% regime.
% In the steady-state regime, the Stokes growth is the strongest at the
% peak power which is at the pulse center of a Gaussian pump pulse.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.95,1.15]*1e-6; % m
Nt = 2^22;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;
%sim.progress_bar = false;
sim.progress_bar_name = 'Find Stokes position';
sim.pulse_centering = false;

num_save = 1;
fiber.L0 = 1; % propagation length
sim.save_period = fiber.L0/num_save;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Gas info
[fiber,sim] = load_default_UPPE_propagate(fiber,sim);

gas.core_radius = 15e-6; % m
gas.temperature = 288; % K
gas.pressure = 30e5; % Pa
gas.wavelength_order = 6;
gas.mode_profile_wavelength = 1030e-9; % m
gas.gas_material = 'H2';
gas.delta = 300e-9; % m; wall thickness of anti-resonant fibers
gas.f_FEM = 1e-2; % loss factor
gas.fiber_type = 'AR_HC_PCF';
gas.xy_sampling = 101;

[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

gas.H2.R.preR = gas.H2.R.preR.*[0,1,0,0,0,0];
gas.H2.V.preR = gas.H2.V.preR*0;

%% Initial condition and Propagate
tfwhm = 20e3; % ps
total_energy = 70e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,1);

prop_output = UPPE_propagate(fiber,input_field,sim,gas);

%%
Stokes_pulse = edgepass_filter('highpass', prop_output, sim.f0, 1080);
Stokes_pulse.fields = Stokes_pulse.fields./sqrt(max(abs(Stokes_pulse.fields).^2,[],1));
pump_pulse = edgepass_filter('lowpass', prop_output, sim.f0, 1080);
pump_pulse.fields = pump_pulse.fields./sqrt(max(abs(pump_pulse.fields).^2,[],1));

figure;
h = plot(t/1e3,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2,'linewidth',2);
set(gca,'fontsize',25);
set(h(1),'Color','k'); set(h(2),'Color','r');
xlim([-1,1]*tfwhm/1e3*2);
xlabel('Time (ns)');
ylabel('Power (norm.)');
legend('P','S');
print('Stokes_position_steady_state','-dpdf');

%%
save('Stokes_position_steady_state.mat');