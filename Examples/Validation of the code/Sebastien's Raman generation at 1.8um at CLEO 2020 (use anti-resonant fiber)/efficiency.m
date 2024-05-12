clearvars; close all;

P_SM = [10,12,15,20,25,27,28,29,30,35,40];
E_SM = [0.2,7.68,14.80,26.59,12.14,6.51,1.85,0.05,0.02,0.03,2.93];

P_MM = [10,15,20,25,27,28,29,30,35];
E_MM = [0.06,12.68,23.04,22.17,20.50,24.24,22.10,19.47,10.50];

figure;
h = plot(P_SM,E_SM); hold on;
h2 = plot(P_MM,E_MM); hold off;
set(h,'linewidth',2); set(h2,'linewidth',2); set(gca,'fontsize',20);
xlim([12,35]);
xlabel('Gas pressure (bar)');
ylabel('Raman energy efficiency (%)');
l = legend('SM','MM: 33% LP_{02}');
set(l,'location','southwest');
print('Raman energy efficiency','-djpeg');