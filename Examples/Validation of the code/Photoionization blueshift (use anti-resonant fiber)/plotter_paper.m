clearvars; close all;

addpath('../../../user_helpers');

spectrum_wavelength = cell(1,4);

load('photoionization_blueshift_3.2uJ.mat');

spectrum_wavelength{1} = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2.*299.792458./(299.792458./f).^2;
spectrum_wavelength{1} = spectrum_wavelength{1}/max(spectrum_wavelength{1});

load('photoionization_blueshift_4.9uJ.mat');

spectrum_wavelength{2} = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2.*299.792458./(299.792458./f).^2;
spectrum_wavelength{2} = spectrum_wavelength{2}/max(spectrum_wavelength{2});

load('photoionization_blueshift_5.1uJ.mat');

spectrum_wavelength{3} = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2.*299.792458./(299.792458./f).^2;
spectrum_wavelength{3} = spectrum_wavelength{3}/max(spectrum_wavelength{3});

load('photoionization_blueshift_5.6uJ.mat');

spectrum_wavelength{4} = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2.*299.792458./(299.792458./f).^2;
spectrum_wavelength{4} = spectrum_wavelength{4}/max(spectrum_wavelength{4});

spectrum_wavelength = cell2mat(spectrum_wavelength) + (0:3);

figure;
h = plot(lambda,spectrum_wavelength,'linewidth',2,'Color','b');
hold on;
plot(lambda,ones(Nt,1).*[1,2,3],'linewidth',2,'LineStyle','--','Color','k');
set(gca,'fontsize',20);
ylim([0,4]);
ylabel('PSD (norm.)');
xlabel('Wavelength (nm)');
xlim([400,900]);
print(gcf,'all_spectra.pdf','-dpdf');

%%
[Keldysh_parameter,~,relative_ne] = calc_photoionization_parameter(prop_output,fiber,sim,gas.gas_material);
fprintf('min Keldysh parameter: %6.4f\n',min(Keldysh_parameter(:)));

[max_peak_power,max_idx] = max(abs(prop_output.fields));
[~,global_max_peak_power_idx] = max(max_peak_power);

figure;
yyaxis left;
plot(t*1e3,Keldysh_parameter(:,1),'linewidth',2,'Color','b');
ylabel('Keldysh parameter');
ylim([0,20]);
yyaxis right;
plot(t*1e3,abs(prop_output.fields(:,1)).^2/1e6,'linewidth',2);
ylabel('Power (MW)');
xlabel('Time (fs)');
xlim([-30,30]);
ax = gca; set(ax,'fontsize',20); set(ax.YAxis(1),'Color','b');
print(gcf,'Keldysh parameter (initial).pdf','-dpdf');

figure;
yyaxis left;
plot(t*1e3,Keldysh_parameter(:,global_max_peak_power_idx),'linewidth',2,'Color','b');
ylabel('Keldysh parameter');
ylim([0,20]);
yyaxis right;
plot(t*1e3,abs(prop_output.fields(:,global_max_peak_power_idx)).^2/1e6,'linewidth',2);
ylabel('Power (MW)');
xlabel('Time (fs)');
xlim([-30,0]);
ax = gca; set(ax,'fontsize',20); set(ax.YAxis(1),'Color','b');
print(gcf,'Keldysh parameter (max peak power).pdf','-dpdf');

figure;
yyaxis left;
plot(t*1e3,relative_ne(:,1),'linewidth',2,'Color','b');
ylabel('n_e (norm.)');
%ylim([0,20]);
yyaxis right;
plot(t*1e3,abs(prop_output.fields(:,1)).^2/1e6,'linewidth',2);
ylabel('Power (MW)');
xlabel('Time (fs)');
xlim([-30,30]);
ax = gca; set(ax,'fontsize',20); set(ax.YAxis(1),'Color','b');
print(gcf,'relative ne (initial).pdf','-dpdf');

figure;
yyaxis left;
plot(t*1e3,relative_ne(:,global_max_peak_power_idx),'linewidth',2,'Color','b');
ylabel('n_e (norm.)');
%ylim([0,20]);
yyaxis right;
plot(t*1e3,abs(prop_output.fields(:,global_max_peak_power_idx)).^2/1e6,'linewidth',2);
ylabel('Power (MW)');
xlabel('Time (fs)');
xlim([-30,0]);
ax = gca; set(ax,'fontsize',20); set(ax.YAxis(1),'Color','b');
print(gcf,'relative ne (max peak power).pdf','-dpdf');