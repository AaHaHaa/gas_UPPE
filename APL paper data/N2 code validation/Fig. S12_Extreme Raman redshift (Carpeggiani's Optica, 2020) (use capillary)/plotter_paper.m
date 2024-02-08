close all; clearvars;

%% Fig. 1d
data = readmatrix('Extreme red shift figure 1d.csv');
data = sortrows(data);

load('redshift_N2_6.0m_530um_4.0bar.mat');
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,end)),1)).^2./lambda.^2;
spectrum = spectrum/max(spectrum);
data_PSD = data(:,2)*abs(trapz(lambda,spectrum)/trapz(data(:,1),data(:,2)));
figure;
plot(lambda,spectrum,'linewidth',2,'Color','b');
hold on;
plot(data(:,1),data_PSD,'linewidth',2,'Color','r');
hold off;
xlim([700,1800]);
set(gca,'fontsize',25);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
l = legend('Ours','Prior work'); set(l,'location','northwest');
print(gcf,'figure_bd.pdf','-dpdf');

%% Fig. 1c
data = readmatrix('Extreme red shift figure 1c.csv');
data = sortrows(data);

load('redshift_N2_5.5m_1000um_0.9bar.mat');
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,end)),1)).^2./lambda.^2;
spectrum = spectrum/max(spectrum);
data_PSD = data(:,2)*max(spectrum)/max(data(:,2));
figure;
plot(lambda,spectrum,'linewidth',2,'Color','b');
hold on;
plot(data(:,1),data_PSD,'linewidth',2,'Color','r');
hold off;
xlim([900,1300]);
set(gca,'fontsize',25);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
l = legend('Ours','Prior work'); set(l,'location','northwest');
print(gcf,'figure_ac.pdf','-dpdf');