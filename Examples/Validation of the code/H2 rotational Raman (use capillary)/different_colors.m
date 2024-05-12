close all;

addpath('../../../user_helpers');

%%
lambda1 = 1540; % nm
f1 = 3e5/lambda1;
pump_pulse = gaussian_spectral_filter(probe, sim.f0, lambda1, 1);

lambda2 = 1698.44; % nm
f2 = 3e5/lambda2;
Stokes_pulse = gaussian_spectral_filter(probe, sim.f0, lambda2, 1);

figure;
h = plot(t,abs([pump_pulse.fields,Stokes_pulse.fields]).^2/1e6);
set(h(1),'Color','k'); set(h(2),'Color','r');
xlim([-20000,20000]);
set(h,'linewidth',2);
set(gca,'fontsize',20);
legend('Pump','Stokes');
xlabel('Time (ps)'); ylabel('Power (MW)');