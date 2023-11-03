clearvars; close all;

addpath('SLMtools');

% refractive index
data = dlmread('n_CH4 (Rollefson).csv',',',1,0);
data_wl = data(:,1); % um
data_n = data(:,2);

%% Sellmeier formula
wavenumber_R = 1./data_wl;
refractivity = @(wavenumber) 0.873e3 + 2.07e-2*(16.50./((3019.6e-4)^2-wavenumber.^2) + 7.5./((1304e-4)^2-wavenumber.^2)); % 10^6(n^2-1)
n_Sellmeier = sqrt(refractivity(wavenumber_R)/1e6 + 1);

%% Extended range of wavelengths (um)
wavelength1 = linspace(0.01,3,1000)';     wavenumber1 = 1./wavelength1;
wavelength2 = linspace(3,4,1000)';        wavenumber2 = 1./wavelength2;
wavelength3 = linspace(4,7,1000)';        wavenumber3 = 1./wavelength3;
wavelength4 = linspace(7,8,1000)';        wavenumber4 = 1./wavelength4;
wavelength5 = linspace(8,15,5000)';       wavenumber5 = 1./wavelength5;

wavelength = [wavelength1;wavelength2;wavelength3;wavelength4;wavelength5];
wavenumber = [wavenumber1;wavenumber2;wavenumber3;wavenumber4;wavenumber5];
n_Sellmeier_calc = sqrt(refractivity(wavenumber)/1e6 + 1);

%% Plot
figure;
h = plot(data_wl,[data_n,n_Sellmeier]);
legend('Data','Fitted Sellmeier');
set(h,'linewidth',2);
xlabel('Wavelength (\mum)'); ylabel('n');
title('Loaded data');
figure;
h = plot(wavelength,n_Sellmeier_calc);
set(h,'linewidth',2);
xlabel('Wavelength (\mum)'); ylabel('n');
title('Extended wavelength');
ylim([1.0003,1.0006]);

%% Save
%save('n_CH4.mat','slm');