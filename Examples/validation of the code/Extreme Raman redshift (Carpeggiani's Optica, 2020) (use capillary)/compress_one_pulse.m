close all;

addpath('../../../user_helpers');

clear probe_pulse

j = size(prop_output.fields,3);

initial.fields = prop_output.fields(:,1,1);
initial.dt = dt;

probe_pulse.fields = prop_output.fields(:,1,j);
probe_pulse.dt = dt;

%%
cutonoff_lambda = 1200; % m
Stokes_probe = edgepass_filter('highpass', probe_pulse, sim.f0, cutonoff_lambda,0.15,1,true);
xlim([0.8,1.9]*1e3);
print(gcf,'spectrum.jpg','-djpeg');

%[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power,fig] = analyze_field( t,f,Stokes_probe.fields,'Treacy-t',pi/6,1e-6,true );
[transform_limited_field,t_interp,TL_fwhm,t_fwhm] = calc_transform_limited( Stokes_probe.fields,10,t );
fprintf('t_fwhm = %4.2f(fs)\n',t_fwhm);
figure;
[~,TCenter] = max(abs(Stokes_probe.fields).^2);
h = plot(t*1e3,circshift(abs(Stokes_probe.fields).^2/1e6,-round(t(TCenter)/(t(2)-t(1))))); hold off;
xlim([-200,200]);
set(h,'linewidth',2);
set(gca,'fontsize',20);
xlabel('Time (fs)');
ylabel('Power (MW)');
title('Raman soliton');
print(gcf,'compressed_pulse_S.jpg','-djpeg');

fprintf('Pulse energy = %4.2f(nJ)\n',sum(abs(Stokes_probe.fields).^2*Stokes_probe.dt/1e3));
fprintf('Energy efficiency = %4.2f%%\n',sum(abs(Stokes_probe.fields).^2)/sum(abs(initial.fields(:,:,1)).^2)*100);