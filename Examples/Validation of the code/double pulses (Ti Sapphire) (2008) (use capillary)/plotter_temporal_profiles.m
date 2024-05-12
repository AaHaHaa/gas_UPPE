close all; %clearvars;

addpath('../../../user_helpers');

%load('Raman_52fs.mat');

%% Single-pulse approach: Input
pump_fields = flipud(prop_output{1}.fields(:,:,1));

smooth_range = 50;
pump_phase = unwrap(angle(pump_fields));
pump_phase = smooth(pump_phase,smooth_range);
pump_inst_freq = -(-(pump_phase(3:end,1)-pump_phase(1:end-2,1))/(2*dt)/(2*pi)) + sim.f0; % THz; I use "central difference" to calculate the slope here
pump_I = abs(pump_fields).^2;

figure;
patch([t(2:end-1);flipud(t(2:end-1))],[zeros(Nt-2,1);flipud(pump_I(2:end-1))],[pump_inst_freq;flipud(pump_inst_freq)],'EdgeColor','k');
colormap(flipud(jet)); caxis([360,390]); set(gca,'Visible','Off');
xlim([-2.2,2.2]);

print(gcf,'single_pulse_input.png','-dpng');

%% Single-pulse approach: Output
cutonoff_lambda = 1000;
filtered_output = edgepass_filter('highpass', prop_output{1}, sim.f0, cutonoff_lambda, 5);
Stokes_fields = flipud(filtered_output.fields(:,:,end)); clearvars filtered_output;
pump_fields = flipud(prop_output{1}.fields(:,:,end)) - Stokes_fields;

smooth_range = 50;
pump_phase = unwrap(angle(pump_fields));
pump_phase = smooth(pump_phase,smooth_range);
pump_inst_freq = -(-(pump_phase(3:end,1)-pump_phase(1:end-2,1))/(2*dt)/(2*pi)) + sim.f0; % THz; I use "central difference" to calculate the slope here
pump_I = abs(pump_fields).^2;

figure;
for i = 60000:71000
    patch([t(i:i+1);flipud(t(i:i+1))],[zeros(2,1);flipud(pump_I(i:i+1))],[pump_inst_freq(i-1:i);flipud(pump_inst_freq(i-1:i))],'EdgeColor','None');
    colormap(flipud(jet)); caxis([360,390]); set(gca,'Visible','Off');
    hold on;
end
plot(t,pump_I,'k');
hold off;
%patch([t(2:end-1);t(end-1:-1:2)],[zeros(Nt-2,1);pump_I(2:end-1)],[flipud(pump_inst_freq);pump_inst_freq],'EdgeColor','k');
%colormap(flipud(jet)); caxis([360,390]); set(gca,'Visible','Off');
xlim([-2.2,2.2]);

hold on;
smooth_range = 50;
Stokes_phase = unwrap(angle(Stokes_fields));
Stokes_phase = smooth(Stokes_phase,smooth_range);
Stokes_inst_freq = -(Stokes_phase(3:end,1)-Stokes_phase(1:end-2,1))/(2*dt)/(2*pi) + sim.f0; % THz; I use "central difference" to calculate the slope here
plot(t,abs(Stokes_fields).^2,'linewidth',5,'Color','r');
hold off;

print(gcf,'single_pulse.png','-dpng');

%% Double-pulse approach
cutonoff_lambda = 1000;
filtered_output = edgepass_filter('highpass', prop_output{2}, sim.f0, cutonoff_lambda, 5);
Stokes_fields = flipud(filtered_output.fields(:,:,end)); clearvars filtered_output;
pump_fields = flipud(prop_output{2}.fields(:,:,end)) - Stokes_fields;

smooth_range = 50;
pump_phase = unwrap(angle(pump_fields));
pump_phase = smooth(pump_phase,smooth_range);
pump_inst_freq = -(-(pump_phase(3:end,1)-pump_phase(1:end-2,1))/(2*dt)/(2*pi)) + sim.f0; % THz; I use "central difference" to calculate the slope here
pump_I = abs(pump_fields).^2;

figure;
for i = 57000:64500
    patch([t(i:i+1);flipud(t(i:i+1))],[zeros(2,1);flipud(pump_I(i:i+1))],[pump_inst_freq(i-1:i);flipud(pump_inst_freq(i-1:i))],'EdgeColor','None');
    colormap(flipud(jet)); caxis([360,390]); set(gca,'Visible','Off');
    hold on;
end
plot(t,pump_I,'k');
hold off;
%patch([t(2:end-1);t(end-1:-1:2)],[zeros(Nt-2,1);pump_I(2:end-1)],[flipud(pump_inst_freq);pump_inst_freq],'EdgeColor','k');
%colormap(flipud(jet)); caxis([360,390]); set(gca,'Visible','Off');
xlim([-8,2]);

hold on;
smooth_range = 50;
Stokes_phase = unwrap(angle(Stokes_fields));
Stokes_phase = smooth(Stokes_phase,smooth_range);
Stokes_inst_freq = -(Stokes_phase(3:end,1)-Stokes_phase(1:end-2,1))/(2*dt)/(2*pi) + sim.f0; % THz; I use "central difference" to calculate the slope here
plot(t,abs(Stokes_fields).^2,'linewidth',5,'Color','r');
hold off;

print(gcf,'double_pulse.png','-dpng');