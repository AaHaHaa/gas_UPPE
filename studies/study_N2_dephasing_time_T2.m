clearvars; close all;

temperature = 273.15 + 25;

eta = linspace(0.5,150,10000)';

rot_T2 = 1e6./(pi*3570*eta); % ps
vib_T2 = 1e6/(pi*22.5)*ones(size(eta)); % ps
            
figure;
yyaxis left;
plot(eta,vib_T2/1e3,'linewidth',2,'LineStyle','-','Color','b');
set(gca,'YColor','b','fontsize',20);
%ylim([0,2]);
xlabel('Gas pressure (bar)');
ylabel('Dephasing time (ns)');
yyaxis right;
plot(eta,rot_T2,'linewidth',2,'LineStyle','-','Color','r');
set(gca,'YColor','r','fontsize',20);
ylabel('Dephasing time (ps)');
legend('Vib','Rot');
xlim([0,30]);
%print('dephasing time.pdf','-dpdf');