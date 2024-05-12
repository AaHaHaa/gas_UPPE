% This code plots the superGaussian pulse shape for the paper to use.

close all; clearvars;

addpath('../../../user_helpers','../../../broadband UPPE algorithm');

%% Setup parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.6,5]*1e-6; % m
Nt = 2^14;
[f0,f_range,time_window,dt] = find_tw_f0(c./wavelength_range,Nt);

sim.f0 = f0;

f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
lambda = c./f*1e9; % nm

%% Pulse
tfwhm1 = 10; % ps
total_energy = 5e3; % nJ
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
input_field = build_MMgaussian(tfwhm1,time_window,total_energy,1,Nt,{'ifft',freq_shift},1,0,10);
%input_field.fields = input_field.fields/max(abs(input_field.fields).^2);

figure;
h = plot(t,abs(input_field.fields).^2,'linewidth',30);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
set(h,'Color','k');
xlim([-15,15]);
print('superGaussian_pulse','-dpdf');