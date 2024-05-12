clearvars; close all;

addpath('../../../user_helpers/');

load('Raman_enhanced_SPM_1030nm_supplement.mat');

Stokes_QE = zeros(length(tfwhm),num_repeat);
dechirped_peak_power = zeros(length(tfwhm),num_repeat);
for ii = 1:length(tfwhm)
    for j = 1:num_repeat
        initial_photon = trapz(f,abs(fftshift(ifft(prop_output{ii,j}.fields(:,:,1)),1)).^2./f);
        Stokes_pulse = gaussian_spectral_filter(prop_output{ii,j}, sim.f0, c/(c/pump_wavelength-gas.H2.V.omega(2)/2/pi)*1e9, 500, 5);

        [~,~,dechirped_field] = pulse_compressor('Treacy-t',22*pi/180,sim.lambda0*1e9,t,Stokes_pulse.fields(:,:,end),1e-3/600);


        Stokes_QE(ii,j) = trapz(f,abs(fftshift(ifft(Stokes_pulse.fields(:,:,end)),1)).^2./f)/initial_photon;
        dechirped_peak_power(ii,j) = max(abs(dechirped_field).^2);
    end
end
Stokes_QE_mean = mean(Stokes_QE,2); Stokes_QE_std = std(Stokes_QE,[],2);
dechirped_peak_power_mean = mean(dechirped_peak_power,2); dechirped_peak_power_std = std(dechirped_peak_power,[],2);

plot_tfwhm = log(tfwhm);

figure;
yyaxis left;
confplot(plot_tfwhm,dechirped_peak_power_mean/1e9,dechirped_peak_power_std/1e9,'linewidth',2,'Color','b','LineStyle','-');
set(gca,'YColor','b');
ylabel('Peak power (GW)');
ylim([0,6]);
yyaxis right;
confplot(plot_tfwhm,Stokes_QE_mean*100,Stokes_QE_std*100,'linewidth',2,'LineStyle','-');
ylabel('QE (%)');
ylim([0,15]);
set(gca,'fontsize',18,'XTick',log([0.2,0.4,0.7,1,2]),'XTickLabel',{'0.2','0.4','0.7','1','2'});
xlim([plot_tfwhm(1),plot_tfwhm(end)]);
xlabel('Pulse duration (ps)');
set(gcf,'Units','Pixels');
fig_pos = get(gcf,'Position'); set(gcf,'Position',[round(fig_pos(1)*0.5),fig_pos(2)*0.5,fig_pos(3)*2,fig_pos(4)*1]);
ax_pos = get(gca,'Position'); set(gca,'Position',[0.1,ax_pos(2),ax_pos(3),ax_pos(4)]); set(gca,'Units','Pixels');
print(gcf,'SPM_1030nm_supplement.pdf','-dpdf','-bestfit');

Stokes_pulse = gaussian_spectral_filter(prop_output{7,1}, sim.f0, c/(c/pump_wavelength-gas.H2.V.omega(2)/2/pi)*1e9, 500, 5);

[~,dechirped_FWHM,dechirped_field] = pulse_compressor('Treacy-t',22*pi/180,sim.lambda0*1e9,t,Stokes_pulse.fields(:,:,end),1e-3/600);
[transform_limited_field,t_insert,transform_limited_FWHM] = calc_transform_limited( Stokes_pulse.fields(:,:,end),3,t );
figure;
plot(t_insert,abs(transform_limited_field).^2/1e9,'linewidth',2,'Color','b');
hold on;
plot(t,abs(dechirped_field).^2/1e9,'linewidth',2,'Color','r');
hold off;
xlim([-0.4,0.4]);
ylim([0,20]);
xlabel('Time (ps)');
ylabel('Power (GW)');
set(gca,'fontsize',25);
legend('TL','D');
print(gcf,'dechirped_1030nm_supplement.pdf','-dpdf');