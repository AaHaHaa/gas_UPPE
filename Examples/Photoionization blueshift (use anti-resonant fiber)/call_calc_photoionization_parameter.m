% This code draws several photoionization-related parameters, such as
% Keldysh parameter, electron density, etc.

clearvars; close all;

addpath('../../user_helpers');

load('photoionization_blueshift_5.2uJ.mat');

[Keldysh_parameter,W,ne,g] = calc_photoionization_parameter(prop_output,fiber,sim,gas.gas_material);

[max_peak_power,max_idx] = max(abs(prop_output.fields));
[~,global_max_peak_power_idx] = max(max_peak_power);

figure;
yyaxis left;
plot(t*1e3,Keldysh_parameter(:,global_max_peak_power_idx),'linewidth',2,'Color','b');
ylabel('Keldysh parameter');
%ylim([-20,20]);
yyaxis right;
plot(t*1e3,ne(:,global_max_peak_power_idx),'linewidth',2);
ylabel('n_e (norm.)');
xlabel('Time (fs)');
%xlim([-30,0]);
ax = gca; set(ax,'fontsize',20); set(ax.YAxis(1),'Color','b');
%print(gcf,'Keldysh parameter_n_e.pdf','-dpdf');