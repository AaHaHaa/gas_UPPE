% This code plots the results of SSFS in H2 or N2, of Fig. S17 in the
% supplement.
% Please change the loading file name below to plot each of their results.

%clearvars;
close all;

addpath('../../../user_helpers');

idx = 101;

load('SSFS_H2.mat');

de = imag(prop_output.delta_permittivity(:,1,1,idx));
de = de/max(de);
P = abs(prop_output.fields(:,:,idx)).^2;
P = P/max(P);
figure;
yyaxis right;
plot(t,de,'linewidth',2);
ylabel('\Delta\epsilon (norm.)');
ylim([-0.5,1.15]);
yyaxis left;
plot(t,P,'linewidth',2,'Color','b');
set(gca,'YColor','b');
ylabel('Power (norm.)');
ylim([0,1.3]);
xlabel('Time (ps)');
set(gca,'fontsize',25);
xlim([26.5,27.5]);
print(gcf,'SSSFS (H2).pdf','-dpdf');

% Spectral evolution
probe.fields = prop_output.fields;
probe.dt = dt;
spectrum_wavelength = abs(fftshift(ifft(probe.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength = squeeze(spectrum_wavelength(:,1,:)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength); log_spectrum_wavelength = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(lambda/1e3,prop_output.z,log_spectrum_wavelength); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim([0.8,1.8]);
c = colorbar; ylabel(c,'PSD (dB)');
xlabel('Wavelength (µm)');
ylabel('Propagation (m)');
set(gca,'fontsize',25,'Units','pixels');
pos = get(gca,'Position');
set(gca,'Position',[pos(1:2),pos(3)*0.95,pos(4)]);
pcolor_plotbox(gca);
print(gcf,'spectral_evolution (H2).pdf','-dpdf');

%%
load('SSFS_N2.mat');

de = imag(prop_output.delta_permittivity(:,1,1,idx));
de = de/max(de);
P = abs(prop_output.fields(:,:,idx)).^2;
P = P/max(P);
figure;
yyaxis right;
plot(t,de,'linewidth',2);
ylabel('\Delta\epsilon (norm.)');
ylim([-0.5,1.15]);
yyaxis left;
plot(t,P,'linewidth',2,'Color','b');
set(gca,'YColor','b');
ylabel('Power (norm.)');
ylim([0,1.3]);
xlabel('Time (ps)');
set(gca,'fontsize',25);
xlim([-60,60]);
print(gcf,'SSSFS (N2).pdf','-dpdf');

% Inset
figure;
yyaxis right;
plot(t,de,'linewidth',10);
ylim([-0.5,1.15]);
set(gca,'YTick',[]);
yyaxis left;
plot(t,P,'linewidth',10,'Color','b');
set(gca,'YColor','b');
ylim([0,1.3]);
set(gca,'fontsize',60);
xlim([34,36]);
set(gca,'Color','None','YTick',[],'XTick',[34,35,36]);
print(gcf,'SSSFS (N2)_supplement_closeview.pdf','-dpdf');

% Spectral evolution
probe.fields = prop_output.fields;
probe.dt = dt;
spectrum_wavelength = abs(fftshift(ifft(probe.fields),1)).^2.*3e2./(3e2./f).^2;
spectrum_wavelength = squeeze(spectrum_wavelength(:,1,:)).';
log_spectrum_wavelength = 10*log10(spectrum_wavelength); log_spectrum_wavelength = log_spectrum_wavelength - max(log_spectrum_wavelength(:));
spectrum = abs(fftshift(ifft(probe.fields),1)).^2;
spectrum = squeeze(spectrum(:,1,:)).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(lambda/1e3,prop_output.z,log_spectrum_wavelength); shading interp;
cmap = whitejet_lower; colormap(cmap); caxis([-30,0]);
xlim([0.8,1.6]);
c = colorbar; ylabel(c,'PSD (dB)');
xlabel('Wavelength (µm)');
ylabel('Propagation (m)');
set(gca,'fontsize',25,'Units','pixels');
pos = get(gca,'Position');
set(gca,'Position',[pos(1:2),pos(3)*0.95,pos(4)]);
pcolor_plotbox(gca);
print(gcf,'spectral_evolution (N2).pdf','-dpdf');