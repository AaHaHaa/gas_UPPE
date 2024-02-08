clearvars; close all;

addpath('../../../../MMTools/gas_UPPE/user_helpers');

%%
total_photon_number_all = zeros(101,10);
photon_number_AS_all    = zeros(101,10);
photon_number_pump_all  = zeros(101,10);
photon_number_S_all     = zeros(101,10);
for ii = 1:10
    load(sprintf('two_color_lossy%u_0.70um2000uJ_10ps_2um5000uJ_10ps_300umID_20bar_L100cm.mat',ii));

    permittivity0 = 8.85418782e-12;
    h = 6.62607015e-34;

    range = {Nt/2:Nt,1:Nt/2};

    % compute spectra
    probe = prop_output;
    probe.fields(range{2},:,:) = 0;

    spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
    spectrum = squeeze(spectrum(:,1,:)).';

    % photon number
    f_AS = f<299792.458/1.0e3 & f>299792.458/1.2e3;
    f_pump = f<299792.458/1.7e3 & f>299792.458/2.3e3;
    f_Stokes = f<299792.458/5e3 & f>299792.458/15e3;
    photon_number = spectrum*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
    total_photon_number_all(:,ii) = trapz(f'*1e12,photon_number,2);
    photon_number_AS_all(:,ii) = trapz(f(f_AS)'*1e12,photon_number(:,f_AS),2);
    photon_number_pump_all(:,ii) = trapz(f(f_pump)'*1e12,photon_number(:,f_pump),2);
    photon_number_S_all(:,ii) = trapz(f(f_Stokes)'*1e12,photon_number(:,f_Stokes),2);
end
total_photon_number_mean = mean(total_photon_number_all,2); total_photon_number_std = std(total_photon_number_all,[],2);
photon_number_AS_mean = mean(photon_number_AS_all,2);       photon_number_AS_std = std(photon_number_AS_all,[],2);
photon_number_pump_mean = mean(photon_number_pump_all,2);   photon_number_pump_std = std(photon_number_pump_all,[],2);
photon_number_S_mean = mean(photon_number_S_all,2);         photon_number_S_std = std(photon_number_S_all,[],2);

figure;
confplot(probe.z,photon_number_pump_mean/total_photon_number_mean(1),photon_number_pump_std/total_photon_number_mean(1),'linewidth',10,'Color','k');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_pump (700nm).pdf','-dpdf');
figure;
confplot(probe.z,photon_number_S_mean/total_photon_number_mean(1),photon_number_S_std/total_photon_number_mean(1),'linewidth',10,'Color','r');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_S (700nm).pdf','-dpdf');
figure;
confplot(probe.z,photon_number_AS_mean/total_photon_number_mean(1),photon_number_AS_std/total_photon_number_mean(1),'linewidth',10,'Color','b');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_AS (700nm).pdf','-dpdf');

%%
total_photon_number_all = zeros(101,10);
photon_number_AS_all    = zeros(101,10);
photon_number_pump_all  = zeros(101,10);
photon_number_S_all     = zeros(101,10);
for ii = 1:10
    load(sprintf('two_color_lossy%u_1.09um2000uJ_10ps_2um5000uJ_10ps_300umID_20bar_L100cm.mat',ii));

    permittivity0 = 8.85418782e-12;
    h = 6.62607015e-34;

    range = {Nt/2:Nt,1:Nt/2};

    % compute spectra
    probe = prop_output;
    probe.fields(range{2},:,:) = 0;

    spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
    spectrum = squeeze(spectrum(:,1,:)).';

    % photon number
    f_AS = f<299792.458/1.0e3 & f>299792.458/1.2e3;
    f_pump = f<299792.458/1.7e3 & f>299792.458/2.3e3;
    f_Stokes = f<299792.458/5e3 & f>299792.458/15e3;
    photon_number = spectrum*time_window^2*1e-24./(h*f'*1e12); % 1/Hz
    total_photon_number_all(:,ii) = trapz(f'*1e12,photon_number,2);
    photon_number_AS_all(:,ii) = trapz(f(f_AS)'*1e12,photon_number(:,f_AS),2);
    photon_number_pump_all(:,ii) = trapz(f(f_pump)'*1e12,photon_number(:,f_pump),2);
    photon_number_S_all(:,ii) = trapz(f(f_Stokes)'*1e12,photon_number(:,f_Stokes),2);
end
total_photon_number_mean = mean(total_photon_number_all,2); total_photon_number_std = std(total_photon_number_all,[],2);
photon_number_AS_mean = mean(photon_number_AS_all,2);       photon_number_AS_std = std(photon_number_AS_all,[],2);
photon_number_pump_mean = mean(photon_number_pump_all,2);   photon_number_pump_std = std(photon_number_pump_all,[],2);
photon_number_S_mean = mean(photon_number_S_all,2);         photon_number_S_std = std(photon_number_S_all,[],2);

figure;
confplot(probe.z,photon_number_pump_mean/total_photon_number_mean(1),photon_number_pump_std/total_photon_number_mean(1),'linewidth',10,'Color','k');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_pump (1090nm).pdf','-dpdf');
figure;
confplot(probe.z,photon_number_S_mean/total_photon_number_mean(1),photon_number_S_std/total_photon_number_mean(1),'linewidth',10,'Color','r');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_S (1090nm).pdf','-dpdf');
figure;
confplot(probe.z,photon_number_AS_mean/total_photon_number_mean(1),photon_number_AS_std/total_photon_number_mean(1),'linewidth',10,'Color','b');
ylim([0,1]);
set(gca,'Color','None','XTick',[],'YTick',[]);
print(gcf,'photon number_AS (1090nm).pdf','-dpdf');