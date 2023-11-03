clearvars; close all;

temperature = 273.15 + 25;

eta = linspace(0.5,150,10000)';

rot_T2 = 1e6./(pi*(6.15./eta+114*eta)); % ps
%vib_T2 = 1e6./(pi*(309./eta+52*eta)); % ps
vib_T2 = 1e6./(pi*(309./eta*(temperature/298)^0.92+(51.8+0.152*(temperature-298)+4.85e-4*(temperature-298)^2)*eta)); % ps

figure;
yyaxis left;
h = plot(eta,[vib_T2,rot_T2]/1e3,'linewidth',2,'LineStyle','-');
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'YColor','b','fontsize',20);
ylim([0,2]);
xlabel('Gas pressure (bar)');
ylabel('Dephasing time (ns)');
yyaxis right;
plot(eta,vib_T2./rot_T2,'linewidth',2,'LineStyle','--');
ylabel('T_{vib}/T_{rot}')
legend('Vib','Rot');
xlim([0,30]);
%print('dephasing time.pdf','-dpdf');