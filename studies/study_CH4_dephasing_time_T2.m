clearvars; close all;

temperature = 273.15 + 25;

eta = linspace(0.5,150,10000)';

vib_T2 = 1e6./(pi*(8220+384*eta)); % ps
          
figure;
plot(eta,vib_T2,'linewidth',2,'LineStyle','-','Color','b');
set(gca,'fontsize',20);
xlabel('Gas pressure (bar)');
ylabel('Dephasing time (ps)');
%print('dephasing time.pdf','-dpdf');