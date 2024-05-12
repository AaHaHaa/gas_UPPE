close all; clearvars;

%% Fig. 3e
data = readmatrix('Fan redshifted figure 3e.csv');
data = sortrows(data);

load('N2_40fs.mat');
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,end)),1)).^2./lambda.^2;
spectrum = spectrum/max(spectrum);
data_PSD = data(:,2)*max(spectrum)/max(data(:,2));
figure;
plot(lambda,spectrum,'linewidth',2,'Color','b');
hold on;
plot(data(:,1),data_PSD,'linewidth',2,'Color','r');
hold off;
xlim([600,1000]);
set(gca,'fontsize',25);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
l = legend('Ours','Prior work'); set(l,'location','northwest');
print(gcf,'figure_3e.pdf','-dpdf');

%% Fig. 3f
data = readmatrix('Fan redshifted figure 3f.csv');
data = sortrows(data);

load('N2_100fs.mat');
spectrum = abs(fftshift(ifft(prop_output.fields(:,1,end)),1)).^2./lambda.^2;
spectrum = spectrum/max(spectrum);
data_PSD = data(:,2)*max(spectrum)/max(data(:,2));
figure;
plot(lambda,spectrum,'linewidth',2,'Color','b');
hold on;
plot(data(:,1),data_PSD,'linewidth',2,'Color','r');
hold off;
xlim([600,1000]);
set(gca,'fontsize',25);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
l = legend('Ours','Prior work'); set(l,'location','northwest');
print(gcf,'figure_3f.pdf','-dpdf');