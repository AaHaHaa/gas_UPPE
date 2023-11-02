clearvars; close all;

addpath('../../../user_helpers');

spectra = cell(1,3);
load('N2_40fs.mat');
spectra{1} = spectrum(end,:).';
load('N2_40fs_noRot.mat');
spectra{2} = spectrum(end,:).';
load('N2_40fs_noVib.mat');
spectra{3} = spectrum(end,:).';

figure;
h = plot(3e5./f,[spectra{1},spectra{2},spectra{3}]);
xlim([600,1000]);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');
l = legend('Vib+Rot','Vib','Rot');
set(l,'fontsize',20,'location','northwest');

%%
spectra = cell(1,3);
load('N2_100fs.mat');
spectra{1} = spectrum(end,:).';
load('N2_100fs_noRot.mat');
spectra{2} = spectrum(end,:).';
load('N2_100fs_noVib.mat');
spectra{3} = spectrum(end,:).';

figure;
h = plot(3e5./f,[spectra{1},spectra{2},spectra{3}]);
xlim([600,1000]);
set(h,'linewidth',2); set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Spectrum (\muJ/nm)');
l = legend('Vib+Rot','Vib','Rot');
set(l,'fontsize',20,'location','northwest');