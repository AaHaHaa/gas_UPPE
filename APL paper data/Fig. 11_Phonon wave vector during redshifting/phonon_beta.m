clearvars; close all;

num = 100;
dynamics = (0:num)';

%% transient regime
beta_P = 10;
beta_S = 8;

beta_P_total = beta_P.*(num-dynamics);
beta_S_total = beta_S.*dynamics;

beta_ph_total = beta_P(1)*num - beta_P_total - beta_S_total;
beta_ph = diff(beta_ph_total);

figure;
yyaxis left;
h = plot(dynamics,[beta_P_total,beta_S_total,beta_ph_total],'linewidth',2,'LineStyle','-');
set(h(1),'Color','k','LineStyle','-'); set(h(2),'Color','r'); set(h(3),'Color','b');
set(gca,'YColor','b','YTick',[])
xlabel('Phonon number N^{ph}'); ylabel('Total wave vector \beta_{total}');
set(gca,'fontsize',25);
ylim([0,beta_P(1)*num*1.3]);
l = legend('Pump','Stokes','Phonon');
yyaxis right;
plot(dynamics(2:num+1),beta_ph,'linewidth',4,'LineStyle','--');
ylabel('\beta^{ph}');
set(gca,'YTick',[]','XTick',[]);
l.String = l.String(1:3); set(l,'fontsize',18,'Box','off');
print(gcf,'transient_phonon.pdf','-dpdf');

%% impulsive regime
beta_P = 10*cos(dynamics/2/pi/15).^2;

beta_P_total = beta_P.*(num-dynamics);

beta_ph_total = beta_P(1)*num - beta_P_total;
beta_ph = diff(beta_ph_total);

figure;
yyaxis left;
h = plot(dynamics,[beta_P_total,beta_ph_total],'linewidth',2,'LineStyle','-');
set(h(1),'Color','k','LineStyle','-'); set(h(2),'Color','b');
set(gca,'YColor','b','YTick',[])
xlabel('Phonon number N^{ph}'); ylabel('Total wave vector \beta_{total}');
set(gca,'fontsize',25);
ylim([0,beta_P(1)*num*1.3]);
l = legend('Pump','Phonon');
yyaxis right;
plot(dynamics(2:num+1),beta_ph,'linewidth',4,'LineStyle','--');
ylabel('\beta^{ph}');
set(gca,'YTick',[]','XTick',[]);
l.String = l.String(1:2); set(l,'fontsize',18,'Box','off');
print(gcf,'impulsive_phonon.pdf','-dpdf');