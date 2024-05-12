clearvars; close all;

addpath('../../../user_helpers');

load('Raman_52fs.mat');

permittivity0 = 8.85418782e-12;
h = 6.62607015e-34;

range = {1:Nt/2, Nt/2:Nt};

%% compute spectra
spectrum = cell(1,2);
j = 1;
for i = 1:2
    
    range_i = range{i};
    
    idx = true(Nt,1); idx(range_i) = false;

    single_pulse.fields = prop_output{1}.fields;
    single_pulse.dt = dt;
    two_pulses.fields = prop_output{2}.fields;
    two_pulses.dt = dt;
    two_pulses.fields(idx,:,:) = 0;

    spectrum_wavelength1 = abs(fftshift(ifft(single_pulse.fields),1)).^2.*299.792458./(299.792458./f).^2;
    spectrum_wavelength1 = squeeze(spectrum_wavelength1(:,1,:)).';
    spectrum1 = abs(fftshift(ifft(single_pulse.fields),1)).^2;
    
    spectrum_wavelength2 = abs(fftshift(ifft(two_pulses.fields),1)).^2.*299.792458./(299.792458./f).^2;
    spectrum_wavelength2 = squeeze(spectrum_wavelength2(:,1,:)).';
    log_spectrum_wavelength = 10*log10(spectrum_wavelength2); log_spectrum_wavelength = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
    spectrum{i} = abs(fftshift(ifft(two_pulses.fields),1)).^2;
    spectrum{i} = squeeze(spectrum{i}(:,1,:)).';
    log_spectrum = 10*log10(spectrum{i}); log_spectrum = log_spectrum - max(log_spectrum(:));

    %{
    figure;
    pcolor(f,prop_output{2}.z*100,log_spectrum); shading interp;
    cmap = whitejet_lower; colormap(cmap); caxis([-60,0]);
    c = colorbar; ylabel(c,'Intensity (dB)');
    set(gca,'fontsize',20);
    xlim([100,600]);
    xlabel('Frequency (THz)');
    ylabel('Propagation distance (cm)');
    title('Spectral evolution');
    %print(gcf,sprintf('Raman_52fs_spectrumMap_%u.jpg',i),'-djpeg');
    %}
    
    ii = size(two_pulses.fields,3);
    z = prop_output{2}.z(ii)*100;

    if ~exist('fig','var')
        fig = figure;
    else
        figure(fig); hold on;
    end
    if i == 1
        plot(1./(299792.458./f*1e-9)/100,spectrum1(:,end).*time_window^2/1e6/(1e12/299792458/100),'linewidth',2,'Color','b');
        hold on;
    end
    hp(j) = plot(1./(299792.458./f*1e-9)/100,spectrum{i}(end,:).'.*time_window^2/1e6/(1e12/299792458/100),'linewidth',2); j = j+1;
    xlim([7600,9000]); ylim([0,0.07]);
    set(gca,'fontsize',20);
    xlabel('Wavenumber (cm^{-1})');
    ylabel('PSD (μJ/cm^{-1})');
    if i == 2 % finishing it
        set(hp(1),'Color','b','LineStyle','--'); set(hp(2),'Color','r','LineStyle','-');
        l = legend('single-pulse','1st Stokes (double-pulse)','2nd Stokes (double-pulse)');
        set(l,'fontsize',12,'location','northwest');
        print(gcf,'Raman_52fs_spectrum.pdf','-dpdf');
    end
end
%{
figure;
plot(t,abs(prop_output{2}.delta_permittivity(:,1,2,1)),'linewidth',2);
hold on; plot(t,abs(prop_output{2}.delta_permittivity(:,1,2,ii)),'linewidth',2); hold off;
l = legend('At 0cm', sprintf('At %3.2fcm',z)); set(l,'location','northwest');
set(gca,'fontsize',20);
xlim([-20,20]);
xlabel('Time (ps)');
ylabel('\Delta\epsilon_V');
%print(gcf,'Raman_52fs_epsilonV.jpg','-djpeg');
%}
%% photon number
f_trigger_AS = f<299792.458/0.5e3 & f>299792.458/0.7e3;
f_trigger_pump = f<299792.458/0.7e3 & f>299792.458/1e3;
f_trigger_S = f<299792.458/1e3 & f>299792.458/1.3e3;
f_AS = f<299792.458/0.5e3 & f>299792.458/0.7e3;
f_pump = f<299792.458/0.7e3 & f>299792.458/1e3;
f_Stokes = f<299792.458/1e3 & f>299792.458/1.3e3;
photon_number = spectrum{1}*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
total_photon_number1 = trapz(f*1e12,photon_number.',1)';
photon_number1_AS = trapz(f(f_trigger_AS)*1e12,photon_number(:,f_trigger_AS).',1)';
photon_number1_pump = trapz(f(f_trigger_pump)*1e12,photon_number(:,f_trigger_pump).',1)';
photon_number1_S = trapz(f(f_trigger_S)*1e12,photon_number(:,f_trigger_S).',1)';
photon_number = spectrum{2}*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
total_photon_number2 = trapz(f*1e12,photon_number.',1)';
photon_number2_AS = trapz(f(f_AS)*1e12,photon_number(:,f_AS).',1)';
photon_number2_pump = trapz(f(f_pump)*1e12,photon_number(:,f_pump).',1)';
photon_number2_S = trapz(f(f_Stokes)*1e12,photon_number(:,f_Stokes).',1)';

fig = figure;
h1 = plot(prop_output{2}.z*100,[total_photon_number1,photon_number1_AS,photon_number1_pump,photon_number1_S]/total_photon_number1(1)+1.2,'linewidth',2);
set(h1(1),'LineStyle','--','Color','k');
set(h1(2),'Color','b'); set(h1(3),'Color','k'); set(h1(4),'Color','r');
hold on;
h2 = plot(prop_output{2}.z*100,[total_photon_number2,photon_number2_AS,photon_number2_pump,photon_number2_S]/total_photon_number2(1),'linewidth',2);
hold off;
set(h2(1),'LineStyle','--','Color','k');
set(h2(2),'Color','b'); set(h2(3),'Color','k'); set(h2(4),'Color','r');
set(gca,'fontsize',20);
ylim([0,2.2]);
set(gca,'YTick',[0,1,1.2,2.2],'YTickLabel',{'0','1','0','1'});
lh1 = legend(h1,'T','AS','P','S','location','northwest');
xlabel('Propagation distance (cm)');
ylabel('Photon number (norm.)');
pos = get(fig,'Position'); set(fig,'Position',[pos(1:2),pos(3)*1.3,pos(4)]);
ax_pos_photon = get(gca,'Position');
print('photon number.eps','-depsc');

%% temporal and phonon wave evolution
fields = permute(sum(abs(prop_output{2}.fields).^2,2),[3,1,2]);
max_fields = max(fields(:));
fields = fields/max_fields;
epsilon_V = abs(squeeze(prop_output{2}.delta_permittivity(:,1,2,:))/permittivity0).';
max_epsilon_V = max(epsilon_V(:));
epsilon_V = epsilon_V/max_epsilon_V;
total_data = [fields; epsilon_V];

fig = figure;
pcolor(t,[prop_output{2}.z,prop_output{2}.z+prop_output{2}.z(end)]*100,total_data);
shading interp;
cmap = whitejet_lower; colormap(cmap);
xlim([-10,10]);
xlabel('Time (ps)');
ylabel('Propagation distance (cm)');
pos = get(fig,'Position'); set(fig,'Position',[pos(1:2),pos(3)*1.5,pos(4)*1.05]);
colorbar; % tmp colorbar for setting sizes below
ax = gca; ax_pos = get(ax,'Position');
set(ax,'YTick',[0,5,10,15,20],'YTickLabel',{'0','5','0/10','5','10'},'Position',[ax_pos_photon(1),ax_pos_photon(2),ax_pos_photon(3)/1.5*1.4,ax_pos_photon(4)/1.05]);
ax_pos = get(ax,'Position');
set(ax,'fontsize',20);
ax_high = axes('Position',[ax_pos(1),ax_pos(2)+ax_pos(4)/2,ax_pos(3),ax_pos(4)/2],'Visible','Off'); colormap(cmap); axes(ax);
set(ax_high,'fontsize',20);
c1 = colorbar(ax); ylabel(c1,'Power');
c1_pos = get(c1,'Position'); set(c1,'Position',[c1_pos(1:3),c1_pos(4)/2*0.8]); set(ax,'Position',ax_pos);
set(c1,'Limits',[0,1],'Ticks',[0,1]);
c2 = colorbar(ax_high); ylabel(c2,'|\Delta\epsilon_r^V|');
set(c2,'Position',[c1_pos(1),c1_pos(2)+0.37,c1_pos(3),c1_pos(4)/2*0.8]); set(ax,'Position',ax_pos);
set(c2,'Limits',[0,1],'Ticks',[0,1]);
hold on;
plot([t(1),t(end)],prop_output{2}.z(end)*100*ones(1,2),'linewidth',2,'Color','k');
hold off;
print('pulseMap.png','-dpng');