% Plot only H2 for the CLEO summary

clearvars; close all;

load('SPM_H2_noRaman.mat');
nonlinear_phase_noRaman = nonlinear_phase;
load('SPM_H2.mat');
figure;
h = plot(tfwhm_all,[nonlinear_phase',nonlinear_phase_noRaman'],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('Pulse duration (ps)');
ylabel('Nonlinear phase (rad)');
set(gca,'fontsize',25);
l = legend('with SRS','without SRS'); legend('boxoff');
pos = get(l,'Position');
set(l,'Position',[pos(1),pos(2)*0.7,pos(3),pos(4)]);
ylim([1.4,4]);
print(gcf,'SPM_transition_H2_CLEO.pdf','-dpdf');