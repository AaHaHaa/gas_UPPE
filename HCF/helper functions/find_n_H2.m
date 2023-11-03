clearvars; close all;

addpath('SLMtools');

% refractive index
data = dlmread('n_H2 (Peck).csv',',',1,0);
data_wl = data(:,1); % um
data_n = data(:,2);

%% Sellmeier formula
wavenumber_R = 1./data_wl;
refractivity = @(wavenumber) 14895.6./(180.7-wavenumber.^2) + 4903.7./(92-wavenumber.^2); % 10^6(n-1)
n_Sellmeier = refractivity(wavenumber_R)/1e6 + 1;

%% Extended range of wavelengths
wavelength_L = linspace(0.01,0.167,1000)';
wavelength_R = [linspace(1.7,5,10),linspace(5.1,15,30)]'; % um
wavenumber_L = 1./wavelength_L;
wavenumber_R = 1./wavelength_R;

n_Sellmeier_calc_L = refractivity(wavenumber_L)/1e6 + 1;
n_Sellmeier_calc_R = refractivity(wavenumber_R)/1e6 + 1;

slm = slmengine([data_wl;wavelength_R],[data_n;n_Sellmeier_calc_R],'knots',[data_wl;wavelength_R],... 
   'decreasing','on','plot','on','extrapolation','cubic','jerk','negative');

n_data_calc_R = slmeval(wavelength_R,slm,0); % refractive index of H2
n_data_calc_L = slmeval(wavelength_L,slm,0); % refractive index of H2

%% Plot
figure;
h = plot(data_wl,[data_n,n_Sellmeier]);
legend('Data','Fitted Sellmeier');
set(h,'linewidth',2);
xlabel('Wavelength (\mum)'); ylabel('n');
title('Loaded data');

figure;
h = plot(wavelength_R,[n_data_calc_R,n_Sellmeier_calc_R]);
legend('Fitted Data','Fitted Sellmeier');
set(h,'linewidth',2);
xlabel('Wavelength (\mum)'); ylabel('n');
title('Extended wavelength');

figure;
h = plot(wavelength_L,[n_data_calc_L,n_Sellmeier_calc_L]);
legend('Fitted Data','Fitted Sellmeier');
set(h,'linewidth',2);
xlabel('Wavelength (\mum)'); ylabel('n');
title('Extended wavelength');

%% Save
save('n_H2.mat','slm');