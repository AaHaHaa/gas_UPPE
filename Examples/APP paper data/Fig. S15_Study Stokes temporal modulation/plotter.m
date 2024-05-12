clearvars; close all;

%% TL no seed
load('Stokes_temporal_modulation_TL_without_seed.mat');
Stokes_pulse = edgepass_filter('highpass', prop_output, sim.f0, 1200, 0.5);
pump_pulse = edgepass_filter('lowpass', prop_output, sim.f0, 1200, 0.5);

figure;
h = plot(t,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2/1e9,'linewidth',2);
set(h(1),'Color','k'); set(h(2),'Color','r');
set(gca,'fontsize',25);
xlabel('Time (ps)');
ylabel('Power (GW)');
xlim([-6,6]);
print(gcf,'TL_without_seed.pdf','-dpdf');

%% TL with seed
load('Stokes_temporal_modulation_TL_with_seed.mat');
Stokes_pulse = edgepass_filter('highpass', prop_output, sim.f0, 1200, 0.5);
pump_pulse = edgepass_filter('lowpass', prop_output, sim.f0, 1200, 0.5);

figure;
h = plot(t,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2/1e9,'linewidth',2);
set(h(1),'Color','k'); set(h(2),'Color','r');
set(gca,'fontsize',25);
xlabel('Time (ps)');
ylabel('Power (GW)');
xlim([-6,6]);
print(gcf,'TL_with_seed.pdf','-dpdf');

%% chirp no seed
load('Stokes_temporal_modulation_chirp_without_seed.mat');
Stokes_pulse = edgepass_filter('highpass', prop_output, sim.f0, 1200, 0.5);
pump_pulse = edgepass_filter('lowpass', prop_output, sim.f0, 1200, 0.5);

figure;
h = plot(t,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2/1e9,'linewidth',2);
set(h(1),'Color','k'); set(h(2),'Color','r');
set(gca,'fontsize',25);
xlabel('Time (ps)');
ylabel('Power (GW)');
xlim([-6,6]);
print(gcf,'chirp_without_seed.pdf','-dpdf');

[~,dechirped_FWHM1,dechirped_field] = pulse_compressor('Treacy-t',22*pi/180,sim.lambda0*1e9,t,Stokes_pulse.fields(:,:,end),1e-3/600);
figure;
plot(t,abs(dechirped_field).^2,'linewidth',10,'Color','b');
xlim([-3,3]);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
print(gcf,'chirp_without_seed_dechirped.pdf','-dpdf');

%% chirp with seed
load('Stokes_temporal_modulation_chirp_with_seed.mat');
Stokes_pulse = edgepass_filter('highpass', prop_output, sim.f0, 1200, 0.5);
pump_pulse = edgepass_filter('lowpass', prop_output, sim.f0, 1200, 0.5);

figure;
h = plot(t,abs([pump_pulse.fields(:,:,end),Stokes_pulse.fields(:,:,end)]).^2/1e9,'linewidth',2);
set(h(1),'Color','k'); set(h(2),'Color','r');
set(gca,'fontsize',25);
xlabel('Time (ps)');
ylabel('Power (GW)');
xlim([-6,6]);
print(gcf,'chirp_with_seed.pdf','-dpdf');

[~,dechirped_FWHM2,dechirped_field] = pulse_compressor('Treacy-t',22*pi/180,sim.lambda0*1e9,t,Stokes_pulse.fields(:,:,end),1e-3/600);
figure;
plot(t,abs(dechirped_field).^2,'linewidth',10,'Color','b');
xlim([-3,3]);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
print(gcf,'chirp_with_seed_dechirped.pdf','-dpdf');