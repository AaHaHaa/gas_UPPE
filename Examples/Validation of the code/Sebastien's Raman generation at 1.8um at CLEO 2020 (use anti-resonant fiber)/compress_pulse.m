close all;

addpath('../../../user_helpers');

clear probe_pulse

%%
probe_pulse.fields = prop_output.fields(:,1,end);
probe_pulse.dt = dt;
probe_pulse = edgepass_filter('highpass', probe_pulse, sim.f0, 1300);

fprintf('Stokes energy=%4.2f%%\n',sum(abs(probe_pulse.fields).^2)*dt/1e3/total_energy*100);
fprintf('Quantum efficiency=%4.2f%%\n',sum(abs(probe_pulse.fields).^2)*dt/1e3/total_energy/(1.03/1.8)*100);

%[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power,fig] = analyze_field( t,f,probe_pulse.fields,'Treacy-t',pi/3,1e-3/600,true,true );

%print(gcf,'compressed_pulse.jpg','-djpeg');