% This code plots the dispersion relation of optical and acoutic branches
% of a diatomic lattice.

clearvars; close all;

C = 1;
M1 = 1;
M2 = 2;
a = 1;

k = linspace(-1.00001,1.00001,100)'*pi/2/a;

omega_optical  = sqrt( C*( (1/M1+1/M2) + sqrt((1/M1+1/M2)^2-4*sin(k*a).^2/M1/M2) ) );
omega_acoustic = sqrt( C*( (1/M1+1/M2) - sqrt((1/M1+1/M2)^2-4*sin(k*a).^2/M1/M2) ) );

figure;
h = plot(k,[omega_acoustic,omega_optical],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
xlabel('k');
ylabel('\omega');
set(gca,'fontsize',25,'YTick',[],'Color','None');
xlim([-1,1]*pi/2/a);
set(gca,'TickLabelInterpreter', 'latex',...
        'XTick',[-1,0,1]*pi/2/a,'XTickLabel',{'$-\frac{\pi}{2a}$','$0$','$\frac{\pi}{2a}$'});
print(gcf,'dispersion.pdf','-dpdf');